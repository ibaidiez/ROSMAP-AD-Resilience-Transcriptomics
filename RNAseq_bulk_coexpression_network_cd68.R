#0. clear workspace and Set working directory
rm(list = ls())
setwd('/Users/ibaidiez/Desktop/projects/ROSMAP/analysis/NatComm')

#Load libraries
library(DESeq2)

metadata_file <- 'metadata/metadata_dlpfc.csv'
rnaseq_data_path <- 'data/Gene_Expression/RNAseq_bulk/brain'

RIN_thr <- 5



################################
#0. Load  metadata
################################

#0.1 Load metadata
metadata = read.csv(metadata_file, header=T, sep='\t', stringsAsFactors = F )
row.names(metadata)<- metadata$sample

metadata$sequencingBatch <- factor(metadata$sequencingBatch)
metadata$msex <- factor(metadata$msex)
metadata$apoe_genotype <- factor(metadata$apoe_genotype)
metadata$ceradsc <- factor(metadata$ceradsc)
metadata$cogdx <- factor(metadata$cogdx)
metadata$dcfdx_lv <- factor(metadata$dcfdx_lv)

#0.2 Filter samples by RIN quality
metadata<-metadata[is.na(metadata$batch)==FALSE,]
metadata<-metadata[is.na(metadata$RIN)==FALSE,]
metadata<-metadata[is.na(metadata$pmi)==FALSE,]
metadata<-metadata[metadata$RIN>=RIN_thr,]


num_samples <- nrow(metadata)


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

cts_b3 = read.table(paste0(rnaseq_data_path,"/ROSMAP_batch3_gene_all_counts_matrix_clean.txt"), header=T,sep='\t' )
colnames(cts_b3) = gsub("X","",colnames(cts_b3)) #Remove Leading X from column names

#1.2 Merge data
cts <- dplyr::left_join(cts_b1, cts_b2, by = 'feature') %>%
  dplyr::left_join(cts_b3, by = 'feature') %>%
  tibble::column_to_rownames(var = "feature")
rm(cts_b1,cts_b2,cts_b3)

#1.3 Remove genes with N_unmapped, N_multimapping, N_noFeature, N_ambiguous names
cts <- cts[ row.names(cts)[ !(row.names(cts) %in%  c('N_unmapped','N_multimapping','N_noFeature','N_ambiguous') )],]
#colnames(cts)



#1.4 Filter and reorder counts based on metadata
cts <- cts[, rownames(metadata)]
#setdiff(rownames(metadata),colnames(cts))


all(rownames(metadata) == colnames(cts)) # check metadata and counts data matchs
all(colnames(cts) == rownames(metadata)) # are they in the same order?


row.names(metadata) <- metadata$projid
colnames(cts)<- metadata$projid




################################
#2. Deseq2
################################

#2.1 construct a DESeqDataSet
#Experimental design to model the samples
# - Allows controlling additional factors (other than the variable of interest)
# - If you have paired samples (if the same subject receives two treatments e.g. before and after treatment), then you need to include the subject (sample) and treatment information in the design formula for estimating the treatment effect while considering differences in subjects.
#   If sample and treatments are represented as subjects and condition in coldata table, then the design formula should be design = ~ subjects + condition. The factor of interest such as condition should go at the end of the formula.


dds <- DESeqDataSetFromMatrix(countData = cts, colData = metadata, design = ~ msex + pmi + age_death + sequencingBatch  )
#dds <- DESeqDataSetFromMatrix(countData = cts, colData = metadata, design = ~ msex + pmi + age_death + sequencingBatch + tdp_st4 + Group )

dds

#2.3 Pre-filtering
# While it is not necessary to pre-filter low count genes before running the DESeq2 functions,
# there are two reasons which make pre-filtering useful: by removing rows in which there are very few reads
#  1. we reduce the memory size of the dds data object,
#  2. we increase the speed of the transformation and testing functions within DESeq2.

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#2.4 Run DESeq
dds <- DESeq(dds)



#2.5 Export normalized expression
vsd <- vst(dds, blind = FALSE) 
nomr_cts <- assay(vsd)

# Save both 'norm_cts' and 'metadata' into the specified file
#save(nomr_cts, metadata, file = 'rnaseq_normdata.RData')



#2.6 Remove Batch/Covariate Effects (The "Adjustment")
library(limma)
design <- model.matrix(~msex + pmi + age_death + sequencingBatch, data = colData(dds))

# We keep the intercept but remove the influence of the specific columns
# If you have a specific condition (e.g. AD vs Control) you want to KEEP, 
# you should include that in the 'design' but exclude it from the removal.
nomr_cts_adj <- removeBatchEffect(nomr_cts, 
                             covariates = metadata[, c("pmi", "age_death")],
                             batch = metadata$sequencingBatch,
                             batch2 = metadata$msex)




#update gene names
genes_symbol <- mapIds(org.Hs.eg.db, keys=gsub("\\..*","",rownames(nomr_cts_adj)), keytype = "ENSEMBL", column = "SYMBOL")
keep_idx <- !is.na(genes_symbol)

nomr_cts_adj <- nomr_cts_adj[keep_idx, ]
genes_symbol <- genes_symbol[keep_idx]
genes_symbol_final <- make.unique(as.character(genes_symbol))

rownames(nomr_cts_adj) <- genes_symbol_final


save(nomr_cts_adj, metadata, file = 'rnaseq_normdata.RData')

rm(cts,design,norm_cts,vsd,genes_symbol,genes_symbol_final,keep,keep_idx)

################################
#3. Co-expression analysis of CD68
################################
#load('rnaseq_normdata.RData')
#cd68_pos <- which(rownames(nomr_cts_adj) == "CD68")

# 3.1 Generate the Co-expression (Correlation) Matrix
coexp_matrix <- cor(t(nomr_cts_adj), method = "pearson")


#3.3 Hierarchical clustering with dynamic cut 
dist_matrix <- as.dist(1 - coexp_matrix)
gene_tree <- hclust(dist_matrix, method = "average") # "average" or "complete" are standard
#plot(gene_tree, labels = FALSE, main = "Gene Clustering Dendrogram", xlab = "", sub = "")


#3.4 Perform dynamic tree cut of the dendrogram rather than fixed k
library(dynamicTreeCut)

dynamic_modules <- cutreeDynamic(dendro = gene_tree, distM = as.matrix(dist_matrix),
                                 deepSplit = 2, pamRespectsDendro = FALSE,
                                 minClusterSize = 30) # Minimum 30 genes per module

# How many modules did it find?
table(dynamic_modules)
num_modules <- length(unique(dynamic_modules))

#3.5 Find wich module has CD68
gene_module_df <- data.frame(
  Gene = rownames(nomr_cts_adj), 
  Module = dynamic_modules
)

# Find the specific module ID for CD68
cd68_module_id <- gene_module_df$Module[gene_module_df$Gene == "CD68"]

cat("CD68 is in Module:", cd68_module_id)

# Extract all other genes in that same module and their correlation values
cd68_module_genes <- gene_module_df$Gene[gene_module_df$Module == cd68_module_id]

cd68_correlations <- coexp_matrix[cd68_module_genes, "CD68"]
cd68_results <- data.frame(
  Gene = names(cd68_correlations),
  Correlation_with_CD68 = as.numeric(cd68_correlations)
)

#ort by highest correlation
cd68_results <- cd68_results[order(-cd68_results$Correlation_with_CD68), ]
write.csv(cd68_results, file = "CD68_Module_Correlations.csv", row.names = FALSE)


#3.6 Perform enrichment of the cluster
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

go_enrich <- enrichGO(gene = cd68_module_genes,
                      OrgDb = org.Hs.eg.db,
                      keyType = "SYMBOL",
                      ont = "BP", # Biological Process
                      pAdjustMethod = "BH",
                      qvalueCutoff = 0.05)

dotplot(go_enrich, showCategory = 15) + ggtitle("GO Enrichment for CD68 Module")






















