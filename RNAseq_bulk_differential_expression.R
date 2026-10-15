#0. clear workspace and Set working directory
rm(list = ls())
setwd('/Users/ibaidiez/Desktop/projects/ROSMAP/analysis/NatComm')

#Load libraries
library(DESeq2)
library(dplyr)
library(org.Hs.eg.db)
library(tibble)

rnaseq_data_path <- 'data/Gene_Expression/RNAseq_bulk/brain'
metadata_file <- "metadata/batch12_dlpfc_metadata.csv"


################################
#1. Load  data
################################

#1.1 Load raw counts data
cts_b1 = read.table(paste0(rnaseq_data_path,"/ROSMAP_batch1_gene_all_counts_matrix_clean.txt"), header=T,sep='\t' )
colnames(cts_b1) = gsub("X","",colnames(cts_b1)) #Remove Leading X from column names

#row.names(cts_b1)<- cts_b1$feature
cts_b2 = read.table(paste0(rnaseq_data_path,"/ROSMAP_batch2_gene_all_counts_matrix_clean.txt"), header=T,sep='\t' )
colnames(cts_b2) = gsub("X","",colnames(cts_b2)) #Remove Leading X from column names
colnames(cts_b2) = gsub("\\.","-",colnames(cts_b2))


#1.2 Merge data
cts <- dplyr::left_join(cts_b1, cts_b2, by = 'feature') %>%
  tibble::column_to_rownames(var = "feature")
rm(cts_b1,cts_b2)

#1.3 Remove genes with N_unmapped, N_multimapping, N_noFeature, N_ambiguous names
cts <- cts[ row.names(cts)[ !(row.names(cts) %in%  c('N_unmapped','N_multimapping','N_noFeature','N_ambiguous') )],]

#row.names(cts) <- cts$feature # change row names
#cts <- as.matrix(cts)

#1.4 Load metadata
metadata = read.csv(metadata_file, header=T, sep='\t', stringsAsFactors = F )
row.names(metadata)<- metadata$sample

metadata$Group <- factor(metadata$Group)
metadata$sequencingBatch <- factor(metadata$sequencingBatch)
metadata$msex <- factor(metadata$msex)
metadata$apoe_genotype <- factor(metadata$apoe_genotype)
metadata$ceradsc <- factor(metadata$ceradsc)
metadata$cogdx <- factor(metadata$cogdx)
metadata$dcfdx_lv <- factor(metadata$dcfdx_lv)


#1.5 Reorder counts based on metadata
cts <- cts[, rownames(metadata)]
#setdiff(rownames(metadata),colnames(cts))

all(rownames(metadata) == colnames(cts)) # check metadata and counts data matchs
all(colnames(cts) == rownames(metadata)) # are they in the same order?


row.names(metadata)<- metadata$projid
colnames(cts)<- metadata$projid



################################
#2. Deseq2
################################

#2.1 construct a DESeqDataSet
#Experimental design to model the samples
# - Allows controlling additional factors (other than the variable of interest)
# - If you have paired samples (if the same subject receives two treatments e.g. before and after treatment), then you need to include the subject (sample) and treatment information in the design formula for estimating the treatment effect while considering differences in subjects.
#   If sample and treatments are represented as subjects and condition in coldata table, then the design formula should be design = ~ subjects + condition. The factor of interest such as condition should go at the end of the formula.


#Option 1: set outlier as another group
cts<-cts[,-which(names(cts) %in% c('9841821','55820335','20964295','10253148','20254588','92393245','31726180'))]
metadata<-metadata[-which(rownames(metadata) %in%  c('9841821','55820335','20964295','10253148','20254588','92393245','31726180')),]


dds <- DESeqDataSetFromMatrix(countData = cts, colData = metadata, design = ~ msex + pmi + age_death + sequencingBatch + Group )
#metadata$tdp_st4[is.na(metadata$tdp_st4)] <- 0
#dds <- DESeqDataSetFromMatrix(countData = cts, colData = metadata, design = ~ msex + pmi + age_death + sequencingBatch + tdp_st4 + Group )

dds


#2.2 select the reference level for condition comparisons.
#    The comparisons of other conditions will be compared against this reference
dds$Group <- relevel(dds$Group, ref = "AD")

#2.3 Pre-filtering
# While it is not necessary to pre-filter low count genes before running the DESeq2 functions,
# there are two reasons which make pre-filtering useful: by removing rows in which there are very few reads
#  1. we reduce the memory size of the dds data object,
#  2. we increase the speed of the transformation and testing functions within DESeq2.

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#2.4 Run DESeq
dds <- DESeq(dds)



################################
#3. Differential gee expression RES vs AD
################################

#Check computed contrasts
resultsNames(dds)


#3.1 Select contrast of interest
res <- results(dds, name="Group_RESIL_vs_AD") 
resOrdered <- res[order(res$pvalue),] #reorder results based on pvalues

#3.2 Convert esembl id to gene ids
df=as.data.frame(resOrdered)
genes_symbol <- mapIds(org.Hs.eg.db, keys=gsub("\\..*","",rownames(df)), keytype = "ENSEMBL", column = "SYMBOL")

df<-add_column(df,genes_symbol , .after = 0)
write.csv(df, file="batch12_ResilvsAD_dlpf_2026.csv")
#write.csv(df, file="batch12_ResilvsAD_dlpf_2026_correctedagesexpmitdp.csv")



#3.3 Save normalized counts
vsd <- vst(dds, blind=FALSE)
nomr_cts <- assay(vsd)

df=as.data.frame(nomr_cts)
genes_symbol <- mapIds(org.Hs.eg.db, keys=gsub("\\..*","",rownames(df)), keytype = "ENSEMBL", column = "SYMBOL")
df<-add_column(df,genes_symbol , .after = 0)
write.csv(df, file="dlpf_batch12_nomalizedcounts.csv")


