#0. clear workspace and Set working directory
rm(list = ls())
setwd('/Users/ibaidiez/Desktop/projects/ROSMAP/analysis/NatComm')

#Load libraries
#library(BPCells)
library(Seurat)
#library(SeuratData)
#library(SeuratObject)
#library(SeuratDisk)
library(data.table)
library(Matrix)
#library(scCustomize)
library(DoubletFinder)
library(ggplot2)
library(patchwork)
library(tidyverse)

exp1_sparse_matrix_file <- '../../data/Gene_Expression/snRNAseq/snRNAseq_exp1/ROSMAP_Brain.snRNAseq_counts_sparse_format_20230420.csv';
exp1_features_file <- '../../data/Gene_Expression/snRNAseq/snRNAseq_exp1/ROSMAP_Brain.snRNAseq_metadata_genes_20230420.csv';
exp1_cells_file <- '../../data/Gene_Expression/snRNAseq/snRNAseq_exp1/ROSMAP_Brain.snRNAseq_metadata_cells_20230420.csv';
biospecimen_file <- '../../data/Metadata/ROSMAP_biospecimen_metadata.csv';

exp2_data_path <- "../../data/Gene_Expression/snRNAseq/snRNAseq_exp2"
out_folder<- 'snRAseq_exp12'
dir.create(out_folder, showWarnings = FALSE)

libraries_to_discard <- c( "200309-B17-A" , "200805-B43-B" )

# needs to be set for large dataset analysis
options(future.globals.maxSize = 1e12)


log_folder<-paste0(out_folder,'/logs')

#Processing parameters
cell_type <- 12 # Exc, Inh, VN, OPCs, Astrc, Oligo, Micr  (KNN cell type)
pca_npcs <- 50
num_features <- 3000
num_dimensions_umap <- 50

#QC
percent_mt_thr <- 10
nCount_RNA_thr <- 400
nFeature_RNA_thr <-253

#outliers removal
iqr_nCount_RNA <- 3
iqr_nFeature_RNA <- 3
iqr_percent.mt <- 3


dir.create(paste0(log_folder,'/',cell_type), showWarnings = FALSE)
dir.create(paste0(log_folder,'/',cell_type,'/qc'), showWarnings = FALSE)


#several
features_qc <- c("P2RY12","CX3CR1",#Surveiling
                 "SPP1","TMEM163","DDIT4","SLC2A3", #antiinflamatory
                 "FTH1","C1QA","APOE", # pro-inflamatory
                 "IFI44L","MX2","GLDN","TPRG1", "PTPRG", # interferon
                 "TOP2A" # profilerative
                 )
features_qc_sep <- c(2.5,6.5,9.5,14.5)

########################## INTEGRATE DATA ################################# 
seurat_list <- c()
layers <- c()
num_layers <- 0

files <- list.files(path = exp2_data_path, pattern = "*.mtx.gz")

#################################################################################
#0.  Cells metadata & remove outliers
#################################################################################
cells_metadata <- read.csv(paste0(log_folder, "/",cell_type,"/qc.csv" ), header=T,stringsAsFactors = F)
rownames(cells_metadata) <-  cells_metadata$X
cells_metadata$X <- NULL

#################################################################################
#1. Load  cells from reference dataset
#################################################################################

#1.1 Load dataset metadata
features <- read.csv(exp1_features_file, header=T,stringsAsFactors = F)
cells <- read.csv(exp1_cells_file, header=T,stringsAsFactors = F)
biospecimen_metadata <- read.csv(biospecimen_file, header=T,stringsAsFactors = F)
biospecimen_metadata <- biospecimen_metadata[biospecimen_metadata$assay == "scrnaSeq",]
rownames(biospecimen_metadata)<- biospecimen_metadata$specimenID
rownames(cells) <- cells$cell_name

#1.2 Load counts 
counts <- fread(exp1_sparse_matrix_file,skip=1)
i=as.vector(counts$V1)
j=as.vector(counts$V2)
x=as.vector(counts$V3)
rm(counts)

counts <- sparseMatrix(i, j, x = x)
rownames(counts)=features$x
colnames(counts)=cells$cell_name
rm(i,j,x)

#1.3 Filter counts to keep cells of interes
batch_cells <-rownames(cells_metadata[ cells_metadata$libraryBatch=="exp1" & cells_metadata$qc=="Pass" ,])
counts<-counts[,batch_cells]
cells<- cells[batch_cells,]


#1.4 Add seurat object
num_layers <- num_layers +1
layers[[num_layers]] <- "exp1"
seurat_list[[num_layers]] <- CreateSeuratObject(counts = counts, meta.data = cells_metadata[batch_cells,])
seurat_list[[num_layers]]$dataset <- "exp1"

rm(batch_cells,counts,features,cells, exp1_sparse_matrix_file,exp1_features_file,exp1_cells_file)
rm(biospecimen_metadata)


#################################################################################
#2. Load  cells from the rest of the datasets
#################################################################################

#2.1 Read counts
for (i in 1:length(files)){
  name <- gsub('.matrix.mtx.gz','', files[i]) #Basename of file
  str<-sprintf("[%d/%d] %s", i, length(files), name)
  print(str)
  
  if ( name %in% libraries_to_discard )
  {
    print(" Skip")
    next;
  }
  
  #2.1.1 Read data
  counts_file <- paste0(exp2_data_path,'/', name ,'.matrix.mtx.gz')
  features_file <- paste0(exp2_data_path,'/', name ,'.features.tsv.gz')
  cells_file <- paste0(exp2_data_path,'/', name ,'.barcodes.tsv.gz')
  
  cts <- ReadMtx(mtx = counts_file,
                 features = features_file,
                 cells = cells_file,
                 cell.column = 1,
                 feature.column = 2,
                 cell.sep = "\t",
                 feature.sep = "\t",
                 skip.cell = 0,
                 skip.feature = 0)
  
  
  # Rename colnames to avoid overlapping names between files
  colnames(cts) <-  paste0(name,'_',colnames(cts))
  
  
  #2.1.2 Filter cells to cell type
  #batch_cells <-rownames(cells_metadata[ cells_metadata$libraryBatch==name & cells_metadata$qc=="Pass" ,])
  batch_cells <-rownames(filter(cells_metadata, libraryBatch==name & qc == "Pass" & doublets == "Singlet" & dist_outliers_nCount_RNA <= iqr_nCount_RNA & dist_outliers_nFeature_RNA <= iqr_nFeature_RNA & dist_outliers_percent.mt  <= iqr_percent.mt   ))
   
  if ( length(batch_cells)<2 )
  {
    print(" Data not found")
    rm(cts)
    next
  }
  
  cts<-cts[,batch_cells]
  

  
  #2.1.3 Create seurat object
  num_layers <- num_layers + 1
  layers[[num_layers]] <-  name
  seurat_list[[num_layers]]  <- CreateSeuratObject(counts = cts, meta.data = cells_metadata[batch_cells,])
  seurat_list[[num_layers]]$dataset <- name
  
  
  rm(batch_cells,cts,str,name)
  

}
  


#################################################################################
#3. Integrate data
#################################################################################

#3.1 Mege data and split layers base on dataset
seurat_merged <-Reduce(function(x,y) merge(x,y,merge.data = FALSE) ,seurat_list)
rm(seurat_list)
seurat_merged[["RNA"]] <- JoinLayers(seurat_merged[["RNA"]])
seurat_merged[["RNA"]] <- split(seurat_merged[["RNA"]], f = seurat_merged$libraryBatch)

#dir.create(paste0(out_folder, "/seurat/",cell_type ) , recursive = TRUE, showWarnings = FALSE)
#SaveSeuratRds(object = seurat_merged,file = paste0(out_folder, "/seurat/",cell_type,"/seurat_object.Rds" ) )

#VlnPlot(seurat_merged, features = "nFeature_RNA")
#VlnPlot(seurat_merged, features = "nCount_RNA")

#3.2 Basic preprocessing
seurat_merged <- NormalizeData(seurat_merged)
seurat_merged <- FindVariableFeatures(seurat_merged,  nfeatures = num_features)
seurat_merged <- ScaleData(seurat_merged)
seurat_merged <- RunPCA(seurat_merged,npcs=pca_npcs)


#3.3 Analisis without integration
#3.3.1 Find clusters and plot data
seurat_merged <- FindNeighbors(seurat_merged, dims = 1:num_dimensions_umap, reduction = "pca")
seurat_merged <- FindClusters(seurat_merged, resolution = 2, cluster.name = "unintegrated_clusters")
seurat_merged <- RunUMAP(seurat_merged, dims = 1:num_dimensions_umap, reduction = "pca", reduction.name = "umap.unintegrated")
#DimPlot(seurat_merged, reduction = "umap.unintegrated", group.by = c( "libraryBatch"))

png(file=paste0(log_folder,'/',cell_type,'/batch_beforecorrection.png'), width=1250, height=400)
p <- DimPlot(seurat_merged, reduction = "umap.unintegrated", group.by = c( "libraryBatch"))
print(p)
dev.off()

#FeaturePlot(seurat_merged, features = c("pANN","dist_outliers_nCount_RNA","dist_outliers_nFeature_RNA","dist_outliers_percent.mt"))

#Check cluster and print histogram of how many cells in each
#DimPlot(seurat_merged, reduction = "umap.unintegrated", group.by = c( "unintegrated_clusters"),label = TRUE, pt.size = 0.5)  + NoLegend()
#tapply(seurat_merged$unintegrated_clusters,seurat_merged$unintegrated_clusters,length)

#Remove small clusters
#cells_rm <- colnames(seurat_merged[,seurat_merged$unintegrated_clusters==53 | seurat_merged$unintegrated_clusters==52 | seurat_merged$unintegrated_clusters==51  ])
#seurat_merged <- seurat_merged[,!colnames(seurat_merged) %in% cells_rm]
#rm(cells_rm)


#3.4 Analisis with integration
#3.4.1 Integrate layersusing CCA Integration
#seurat_merged <- IntegrateLayers(object = seurat_merged, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",k.weight = FALSE)
seurat_merged <- IntegrateLayers(object = seurat_merged, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca")

#SaveSeuratRds(object = seurat_merged,file = paste0(out_folder, "/seurat/",cell_type,"/seurat_integrated.Rds" ) )

#3.4.2 Re-join layers after integration
seurat_merged[["RNA"]] <- JoinLayers(seurat_merged[["RNA"]])

#3.4.3 Perform clustering
seurat_merged <- FindNeighbors(seurat_merged, reduction = "integrated.cca", dims = 1:num_dimensions_umap)
seurat_merged <- FindClusters(seurat_merged, resolution = 1)
#seurat_merged <- FindClusters(seurat_merged, resolution = 0.2)

# Visualization
seurat_merged <- RunUMAP(seurat_merged, dims = 1:num_dimensions_umap, reduction = "integrated.cca")

png(file=paste0(log_folder,'/',cell_type,'/batch_corrected.png'), width=1250, height=400)
p <- DimPlot(seurat_merged, reduction = "umap", group.by = c( "libraryBatch"))
print(p)
dev.off()


png(file=paste0(log_folder,'/',cell_type,'/clusters.png'), width=600, height=500)
p <- DimPlot(seurat_merged, reduction = "umap", group.by = c("seurat_clusters") ,label = TRUE, pt.size = 0.5)
print(p)
dev.off()


png(file=paste0(log_folder,'/',cell_type,'/feature_plot_qc.png'), width=600, height=500)
p <-FeaturePlot(seurat_merged, reduction = "umap", features = c("pANN","dist_outliers_nCount_RNA","dist_outliers_nFeature_RNA","dist_outliers_percent.mt"))
print(p)
dev.off()




#Number of counts and genes
png(file=paste0(log_folder,'/',cell_type,'/qc_all.png'), width=2000, height=800)
seurat_merged <- PercentageFeatureSet(seurat_merged, pattern = "^MT-", col.name = "percent.mt")

df<- seurat_merged@meta.data[,c("nCount_RNA","nFeature_RNA","percent.mt","seurat_clusters")]
p1<- ggplot(df, aes( x=seurat_clusters, y = nCount_RNA)) +
  ggdist::stat_halfeye(adjust = .5,  width = .6, .width = 0,  justification = -.3,  point_colour = NA) + 
  geom_boxplot( width = .25,  outlier.shape = NA ) +
  geom_point( size = 1.5, alpha = .2, position = position_jitter(seed = 1, width = .1 ) ) +
  geom_hline(yintercept=nCount_RNA_thr, linetype="dashed", color = "red") +
  theme_bw()
p2<- ggplot(df, aes( x=seurat_clusters, y = nFeature_RNA)) +
  ggdist::stat_halfeye(adjust = .5,  width = .6, .width = 0,  justification = -.3,  point_colour = NA) + 
  geom_boxplot( width = .25,  outlier.shape = NA ) +
  geom_point( size = 1.5, alpha = .2, position = position_jitter(seed = 1, width = .1 ) ) +
  geom_hline(yintercept=nFeature_RNA_thr, linetype="dashed", color = "red") +
  theme_bw()
p3<- ggplot(df, aes( x=seurat_clusters, y = percent.mt)) +
  ggdist::stat_halfeye(adjust = .5,  width = .6, .width = 0,  justification = -.3,  point_colour = NA) + 
  geom_boxplot( width = .25,  outlier.shape = NA ) +
  geom_point( size = 1.5, alpha = .2, position = position_jitter(seed = 1, width = .1 ) ) +
  theme_bw()

p <- p1 + p2 + p3 + plot_layout(ncol = 3)
print(p)
dev.off()



#Plot markers short
png(file=paste0(log_folder,'/',cell_type,'/markers.png'), width=800, height=600)
p<-DotPlot(seurat_merged, features = features_qc) +
  geom_vline(xintercept = features_qc_sep, linetype="dashed") +
  RotatedAxis()
print(p)
dev.off()

#Plot markers long
png(file=paste0(log_folder,'/',cell_type,'/markers_long.png'), width=2000, height=800)
features<-c('TOP2A', #proliferative
            'CYFIP1','ARHGAP22','P2RY12','CX3CR1', #surveiling
            'CALM2','PRNP','CLU','CALM3','CALM1',
            'CST3', 'PTPRC','TNFRSF1B',
            'IL17RA','IL4R', 'CD14','IL1B',
            'CEBPD', 'TMEM163','SPP1','VSIG4','TLR2', #Reactive
            'FTH1', 'FTL', #enhanced redox
            'CASP1', 'CASP5', 'DNAJB6','DNAJB1', 'CCND3', 'NLRP1', 'TGFBR1', #stress response
            'SPIDR', 'PRKCE', 'PPARG', 'MSR1', 'CXCR4', 'ATG7', 'GPNMB', 'ABCA1', 'TREM2','APOE',  #lipid associated
            'ITPR2', 'TGFBR2','ITPR1',
            'OAS2', 'PNPT1', 'OAS3', 'IFIH1', 'OAS1', 'DDX58','MX2', 'STAT1', 'IFIT3', 'IFIT2',  #interferon
            'ITGAX', 'CCL3L1', 'CCL4', 'CCL3', 'NFKB1', 'EGR3', 'NFKBID', 'FOS', 'JUNB','NLRP3', #inflammatory
            'LYN', 'ELMO1', 'AKT3', 'SERPINE1', 'RUNX1', 'PIK3R1'
)
p<-DotPlot(seurat_merged, features = features) +
  geom_vline(xintercept = c(1.5,5.5,10.5,13.5,17.5,22.5,24.5,31.5,41.5,44.5,54.5,64.5), linetype="dashed") +
  RotatedAxis()
print(p)
dev.off()


#Save data
dir.create(paste0(out_folder,'/seurat_integrated'), showWarnings = FALSE)
SaveSeuratRds(object = seurat_merged,file = paste0(out_folder, "/seurat_integrated/",cell_type,".Rds" ) )


#3.4.4 cell proportions
df<-t(table(seurat_merged@meta.data$seurat_clusters, seurat_merged@meta.data$individualID))
dir.create(paste0(out_folder,'/individual_stats'), showWarnings = FALSE)
dir.create(paste0(out_folder,'/individual_stats/cell_proportions'), showWarnings = FALSE)
write.csv(df, file= paste0(out_folder, "/individual_stats/cell_proportions/",cell_type,".csv" ))
rm(df)




#Generate lower resolution clusters
seurat_merged <- LoadSeuratRds(paste0(out_folder, "/seurat_integrated/",cell_type,".Rds" ))
seurat_merged <- FindClusters(seurat_merged, resolution = 0.2)

# Visualization
seurat_merged <- RunUMAP(seurat_merged, dims = 1:num_dimensions_umap, reduction = "integrated.cca")
png(file=paste0(log_folder,'/',cell_type,'/clusters_small.png'), width=600, height=500)
p <- DimPlot(seurat_merged, reduction = "umap", group.by = c("seurat_clusters") ,label = TRUE, pt.size = 0.5)
print(p)
dev.off()

#Plot markers short
png(file=paste0(log_folder,'/',cell_type,'/clusters_small_markers.png'), width=800, height=600)
p<-DotPlot(seurat_merged, features = features_qc) +
  geom_vline(xintercept = features_qc_sep, linetype="dashed") +
  RotatedAxis()
print(p)
dev.off()

#Plot markers long
png(file=paste0(log_folder,'/',cell_type,'/clusters_small_markers_long.png'), width=2000, height=800)
features<-c('TOP2A', #proliferative
            'CYFIP1','ARHGAP22','P2RY12','CX3CR1', #surveiling
            'CALM2','PRNP','CLU','CALM3','CALM1',
            'CST3', 'PTPRC','TNFRSF1B',
            'IL17RA','IL4R', 'CD14','IL1B',
            'CEBPD', 'TMEM163','SPP1','VSIG4','TLR2', #Reactive
            'FTH1', 'FTL', #enhanced redox
            'CASP1', 'CASP5', 'DNAJB6','DNAJB1', 'CCND3', 'NLRP1', 'TGFBR1', #stress response
            'SPIDR', 'PRKCE', 'PPARG', 'MSR1', 'CXCR4', 'ATG7', 'GPNMB', 'ABCA1', 'TREM2','APOE',  #lipid associated
            'ITPR2', 'TGFBR2','ITPR1',
            'OAS2', 'PNPT1', 'OAS3', 'IFIH1', 'OAS1', 'DDX58','MX2', 'STAT1', 'IFIT3', 'IFIT2',  #interferon
            'ITGAX', 'CCL3L1', 'CCL4', 'CCL3', 'NFKB1', 'EGR3', 'NFKBID', 'FOS', 'JUNB','NLRP3', #inflammatory
            'LYN', 'ELMO1', 'AKT3', 'SERPINE1', 'RUNX1', 'PIK3R1'
)
p<-DotPlot(seurat_merged, features = features) +
  geom_vline(xintercept = c(1.5,5.5,10.5,13.5,17.5,22.5,24.5,31.5,41.5,44.5,54.5,64.5), linetype="dashed") +
  RotatedAxis()
print(p)
dev.off()

#proportions small
df<-t(table(seurat_merged@meta.data$seurat_clusters, seurat_merged@meta.data$individualID))
dir.create(paste0(out_folder,'/individual_stats'), showWarnings = FALSE)
dir.create(paste0(out_folder,'/individual_stats/cell_proportions'), showWarnings = FALSE)
write.csv(df, file= paste0(out_folder, "/individual_stats/cell_proportions/",cell_type,"_small.csv" ))
rm(df)



