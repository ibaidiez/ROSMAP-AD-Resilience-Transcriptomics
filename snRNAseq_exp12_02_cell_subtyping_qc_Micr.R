
#0. clear workspace and Set working directory
rm(list = ls())
setwd('/Users/ibaidiez/Desktop/projects/ROSMAP/analysis/NatComm')


#Load libraries
library(Seurat)
library(data.table)
library(Matrix)
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
percent.mt_up_thr <- 10



multiplet_rate<- 0.076 # For douublets. estimate based on next table https://github.com/chris-mcginnis-ucsf/DoubletFinder/issues/76


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
#seurat_list <- c()
#layers <- c()
#num_layers <- 0

#################################################################################
#0.  Cells metadata & remove outliers
#################################################################################
files <- list.files(path = exp2_data_path, pattern = "*.mtx.gz")

#0.1 Load metadata of cells
cells_metadata<-c()
for (i in 1:length(files)){
  name <- gsub('.matrix.mtx.gz','', files[i]) #Basename of file
  str<-sprintf("[%d/%d] %s", i, length(files), name)
  print(str)
  if  (file.exists(paste0(out_folder, "/cells_types/",name,".csv"))==FALSE ) {
    next
  }
  
  celltps <- read.csv( paste0(out_folder, "/cells_types/",name,".csv") , header=T,stringsAsFactors = F)
  cells_metadata<-rbind(cells_metadata,celltps)
  rm(celltps)
}
rownames(cells_metadata) <- cells_metadata$X

cells_metadata <- cells_metadata[cells_metadata$knn_clust==cell_type,c("individualID","libraryBatch","sequencingBatch","refUMAP_1","refUMAP_2","clustSeur.predicted.id")]
cells_metadata$qc <- 'Pass'
colnames(cells_metadata) <- c("individualID","libraryBatch","sequencingBatch","refUMAP_1","refUMAP_2","seurat_clusters","qc")


#0.2 Include reference cells
cells_metadata_ref <- read.csv(paste0(out_folder,'/cells_types/exp1.csv'), header=T,stringsAsFactors = F)
rownames(cells_metadata_ref) <- cells_metadata_ref$X

#cells_metadata_ref <- cells_metadata_ref[cells_metadata_ref$clusters7cells==cell_type,c("individualIDs","libraryBatch","sequencingBatch","umap_1","umap_2")]
cells_metadata_ref <- cells_metadata_ref[cells_metadata_ref$seurat_clusters==cell_type,c("individualIDs","libraryBatch","sequencingBatch","umap_1","umap_2","seurat_clusters")]
cells_metadata_ref$qc <- 'Pass'
colnames(cells_metadata_ref) <- c("individualID","libraryBatch","sequencingBatch","refUMAP_1","refUMAP_2","seurat_clusters","qc")

cells_metadata <-rbind(cells_metadata_ref,cells_metadata)
rm(cells_metadata_ref)


#0.3 Identify outliers
library("dbscan")
k<-100
ds_thr=4.5

cellt_embd<-as.matrix(cells_metadata[,c("refUMAP_1","refUMAP_2")])
d<-rowMeans(kNNdist(cellt_embd, k = k, all = TRUE))
ds<-scale(log(d))
cells_metadata[ds>=ds_thr,]$qc <- "Outlier"
rm(cellt_embd,d,ds,k,name,str,i)

#0.4 Check mismatch between knn and seurat cluster
cells_metadata[cells_metadata$seurat_clusters!=cell_type,]$qc <- "CellMismatch"
print(sprintf("Microglia: %d [%d outliers %.2f %%]",nrow(cells_metadata),sum(cells_metadata$qc=="Outlier"),100*sum(cells_metadata$qc=="Outlier")/nrow(cells_metadata)))
print(sprintf("Microglia: %d [%d mismatch %.2f %%]",nrow(cells_metadata),sum(cells_metadata$qc=="CellMismatch"),100*sum(cells_metadata$qc=="CellMismatch")/nrow(cells_metadata)))
cells_metadata$doublets <- ""
cells_metadata$pANN <- ""
cells_metadata$nCount_RNA <- ""
cells_metadata$nFeature_RNA <- ""
#cells_metadata$kdist_nCount_RNA <- 0
#cells_metadata$kdist_nFeature_RNA <- 0
#cells_metadata$kdist_percent.mt <- 0
cells_metadata$dist_outliers_nCount_RNA <- 0
cells_metadata$dist_outliers_nFeature_RNA <- 0
cells_metadata$dist_outliers_percent.mt <- 0



png(file=paste0(log_folder,'/',cell_type,'/cluster_outliers_mismatch.png'), width=500, height=400)

p<- ggplot() +
  geom_point(aes(x = refUMAP_1, y = refUMAP_2), data = data.frame(cells_metadata[cells_metadata$qc=="Pass",]), size=1,alpha=0.5, color= rgb(0.75,0.75,0.75) ) + # correctly classified
  geom_point(aes(x = refUMAP_1, y = refUMAP_2), data = data.frame(cells_metadata[cells_metadata$qc=="Outlier",]), size=1,alpha=0.5, color= rgb(228,26,28,maxColorValue = 255) ) + # Outliers 
  geom_point(aes(x = refUMAP_1, y = refUMAP_2), data = data.frame(cells_metadata[cells_metadata$qc=="CellMismatch",]), size=1,alpha=0.5, color= rgb(33,102,172,maxColorValue = 255) ) + # Mismatch 
  xlab('UMAP1') + ylab('UMAP2') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

print(p)
dev.off()


cells_metadata$seurat_clusters <- NULL


#################################################################################
#1. Load  cells from reference dataset
#################################################################################

##1.1 Read cell type info and filter  only the cell type of interest
#cells_metadata <- read.csv(paste0(out_folder,'/cells_types/exp1.csv'), header=T,stringsAsFactors = F)
#cells_metadata <- cells_metadata[cells_metadata$clusters7cells==cell_type,]

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

#1.4 Create seurat object 
seurat_object  <- CreateSeuratObject(counts = counts, meta.data =cells_metadata[batch_cells,c("individualID","libraryBatch","sequencingBatch","qc")])
rm(batch_cells,counts)

#1.5 Compute percentage MT genes and remove low quality genes
seurat_object <- PercentageFeatureSet(seurat_object, pattern = "^MT-", col.name = "percent.mt")

#p1 <- VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#p2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = 'lm')
#p <- p1 + p2
#print(p)


#df<-seurat_object@meta.data[,c("nFeature_RNA", "nCount_RNA", "percent.mt")]
#df$x <- 1
#png(file=paste0(log_folder,'/',cell_type,'/qc_all/exp1.png'), width=600, height=600)
#p1<- ggplot(df, aes( x=x, y = nFeature_RNA)) +
#  ggdist::stat_halfeye(adjust = .5,  width = .6, .width = 0,  justification = -.3,  point_colour = NA) + 
#  geom_boxplot( width = .25,  outlier.shape = NA ) +
#  geom_point( size = 1.5, alpha = .2, position = position_jitter(seed = 1, width = .1 ) ) +
#  theme_bw()
#p2<- ggplot(df, aes( x=x, y = nCount_RNA)) +
#  ggdist::stat_halfeye(adjust = .5,  width = .6, .width = 0,  justification = -.3,  point_colour = NA) + 
#  geom_boxplot( width = .25,  outlier.shape = NA ) +
#  geom_point( size = 1.5, alpha = .2, position = position_jitter(seed = 1, width = .1 ) ) +
#  theme_bw()
#p3<- ggplot(df, aes( x=x, y = percent.mt)) +
#  ggdist::stat_halfeye(adjust = .5,  width = .6, .width = 0,  justification = -.3,  point_colour = NA) + 
#  geom_boxplot( width = .25,  outlier.shape = NA ) +
#  geom_point( size = 1.5, alpha = .2, position = position_jitter(seed = 1, width = .1 ) ) +
#  theme_bw()
#p4<-ggplot(df, aes(x = nCount_RNA, y = nFeature_RNA)) + 
#  geom_point( size = 1.5, alpha = .2) + 
#  geom_smooth(method = 'lm') +
#  theme_bw()
#p<- p1 + p2 + p3 + p4
#print(p)
#dev.off()
#rm(df,p,p1,p2,p3,p4)



#seurat_object <- subset(seurat_object, subset = nCount_RNA > 400 & nFeature_RNA > 253  & percent.mt < 2)
########seurat_object <- subset(seurat_object,  percent.mt < 10)

#1.6 . Remove doubles (DoubletFinder)
#1.6.1 Basic preprocessing
seurat_object <- NormalizeData(seurat_object)
seurat_object <- FindVariableFeatures(seurat_object, nfeatures = num_features) 
seurat_object <- ScaleData(seurat_object) 
seurat_object <- RunPCA(seurat_object, npcs=pca_npcs)
seurat_object <- FindNeighbors(seurat_object, dims = 1:num_dimensions_umap) #construction of k-NN neighbor graph
seurat_object <- FindClusters(seurat_object,resolution=1.5)
seurat_object <- RunUMAP(seurat_object, dims = 1:num_dimensions_umap)


#1.6.2 pK Identification (no ground-truth data available to compute it)
sweep.res.list_cells <- paramSweep(seurat_object, PCs = 1:pca_npcs,sct = FALSE) # Performs steps 1,2,3 4. Provides list of the proportion of artificial neighbours for combination of pN and pK
sweep.stats <- summarizeSweep(sweep.res.list_cells, GT = FALSE) # Summarize de results
bcmvn <- find.pK(sweep.stats) # metric to compute best pk (highest value)
#ggplot(bcmvn, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line() # plot results

pK <- bcmvn %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
  filter(BCmetric == max(BCmetric)) %>%
  dplyr::select(pK)
pK <- as.numeric(as.character(pK[[1]]))
rm(bcmvn,sweep.res.list_cells,sweep.stats)

#1.6.3 Homotypic Doublet Proportion Estimate
annotations <- seurat_object@meta.data$seurat_clusters # Get the cluster of each cell
homotypic.prop <- modelHomotypic(annotations)         #compute perfect of homotopic doublets
nExp_poi <- round(multiplet_rate*nrow(seurat_object@meta.data))  #Expected number of doublets. Assuming 7.5% doublet formation rate - tailor for your dataset
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

#1.6.4 Run doublet pfinder with estimated parameters
#seurat_object <- doubletFinder(seurat_object,
#                               PCs = 1:pca_npcs,
#                              pN = 0.25,
#                              pK = pK,
#                              nExp = nExp_poi.adj,
#                              reuse.pANN = FALSE, sct = FALSE)

seurat_object <- doubletFinder(seu = seurat_object, 
                               PCs = 1:pca_npcs,
                               pN = 0.25,
                               pK = pK, # the neighborhood size used to compute the number of artificial nearest neighbours
                               nExp = nExp_poi.adj) # number of expected real doublets



##1.6.5 Compute ncount, nfeature, andpercent.mt distance of each cell in each cluster

clusts <- levels(seurat_object@meta.data$seurat_clusters)
dist_outliers_nCount_RNA <- cells_metadata[rownames(seurat_object@meta.data),]$dist_outliers_nCount_RNA
dist_outliers_nFeature_RNA <- cells_metadata[rownames(seurat_object@meta.data),]$dist_outliers_nFeature_RNA
dist_outliers_percent.mt  <- cells_metadata[rownames(seurat_object@meta.data),]$dist_outliers_percent.mt
for(j in 1:length(clusts) )
{
  indxc <- seurat_object@meta.data$seurat_clusters==clusts[j]
  
  #nCountRNA
  x <- seurat_object@meta.data[indxc,]$nCount_RNA
  dist_outliers_nCount_RNA[indxc] <- (  x - quantile(x,.75) ) / IQR(x)
  
  #nFeatureRNA
  x <- seurat_object@meta.data[indxc,]$nFeature_RNA
  dist_outliers_nFeature_RNA[indxc] <- (  x - quantile(x,.75) ) / IQR(x)
  
  
  #percent.mt
  x <- seurat_object@meta.data[indxc,]$percent.mt
  iqr<-IQR(x)
  if (iqr < 2 )
  {
    iqr <-2
  }
  dist_outliers_percent.mt[indxc] <- (  x - quantile(x,.75) ) / iqr
  
  
}

cells_metadata[rownames(seurat_object@meta.data),]$dist_outliers_nCount_RNA <- dist_outliers_nCount_RNA 
cells_metadata[rownames(seurat_object@meta.data),]$dist_outliers_nFeature_RNA <- dist_outliers_nFeature_RNA
cells_metadata[rownames(seurat_object@meta.data),]$dist_outliers_percent.mt <- dist_outliers_percent.mt

rm(clusts,indxc,x,dist_outliers_nCount_RNA,dist_outliers_nFeature_RNA,dist_outliers_percent.mt)




# visualize doublets
#png(file=paste0(log_folder,'/',cell_type,'/doublets_exp1.png'), width=600, height=500)
#p<-DimPlot(seurat_object, reduction = 'umap', group.by = colnames(seurat_object@meta.data)[ncol(seurat_object@meta.data)]) + plot_annotation(title = 'exp1')
#print(p)
#dev.off()

##table(seurat_object@meta.data[,ncol(seurat_object@meta.data)]) # number of singlets and doublets
library(ggdist)
png(file=paste0(log_folder,'/',cell_type,'/qc/exp1.png'), width=2000, height=800)
df<- seurat_object@meta.data[,c("nCount_RNA","nFeature_RNA","percent.mt","seurat_clusters")]
df$outliers_nCount_RNA<- cells_metadata[rownames(seurat_object@meta.data),]$dist_outliers_nCount_RNA> 3
df$outliers_nFeature_RNA<- cells_metadata[rownames(seurat_object@meta.data),]$dist_outliers_nFeature_RNA> 3
df$outliers_percent.mt<- cells_metadata[rownames(seurat_object@meta.data),]$dist_outliers_percent.mt > 3

p1<-DimPlot(seurat_object, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
p2<-DimPlot(seurat_object, reduction = 'umap', group.by = colnames(seurat_object@meta.data)[ncol(seurat_object@meta.data)]) + plot_annotation(title = "exp1")
p3<-FeaturePlot(seurat_object, features = colnames(seurat_object@meta.data)[ncol(seurat_object@meta.data)-1])
p4<- ggplot(df, aes( x=seurat_clusters, y = nCount_RNA)) +
  ggdist::stat_halfeye(adjust = .5,  width = .6, .width = 0,  justification = -.3,  point_colour = NA) + 
  geom_boxplot( width = .25,  outlier.shape = NA ) +
  geom_point( data =df, aes( x=seurat_clusters, y = nCount_RNA, color = outliers_nCount_RNA),size = 1.5, alpha = .2, position = position_jitter(seed = 1, width = .1 ) ) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  geom_hline(yintercept=nCount_RNA_thr, linetype="dashed", color = "red") +
  theme_bw()
p5<-  ggplot(df, aes( x=seurat_clusters, y = nFeature_RNA)) +
  ggdist::stat_halfeye(adjust = .5,  width = .6, .width = 0,  justification = -.3,  point_colour = NA) + 
  geom_boxplot( width = .25,  outlier.shape = NA ) +
  geom_point( data =df, aes( x=seurat_clusters, y = nFeature_RNA, color = outliers_nFeature_RNA),size = 1.5, alpha = .2, position = position_jitter(seed = 1, width = .1 ) ) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  geom_hline(yintercept=nFeature_RNA_thr, linetype="dashed", color = "red") +
  theme_bw()
p6<- ggplot(df, aes( x=seurat_clusters, y = percent.mt)) +
  ggdist::stat_halfeye(adjust = .5,  width = .6, .width = 0,  justification = -.3,  point_colour = NA) + 
  geom_boxplot( width = .25,  outlier.shape = NA ) +
  geom_point( data =df, aes( x=seurat_clusters, y = percent.mt, color = outliers_percent.mt),size = 1.5, alpha = .2, position = position_jitter(seed = 1, width = .1 ) ) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
  geom_hline(yintercept=percent.mt_up_thr, linetype="dashed", color = "red") +
  theme_bw()

p7<-DotPlot(seurat_object, features = features_qc) +
  geom_vline(xintercept = features_qc_sep, linetype="dashed") +
  RotatedAxis()

p8<-ggplot(df, aes(x = nCount_RNA, y = nFeature_RNA)) + 
  geom_point( size = 1.5, alpha = .2) + 
  geom_smooth(method = 'lm') +
  theme_bw()
p <- p1 + p3 + p2 +p7 + p4 +p5 + p6 +p8  + plot_layout(ncol = 4)
print(p)
dev.off()
rm(df,p,p1,p2,p3,p4,p5,p6,p7,p8)


##1.6.6 Remove doublets
cells_metadata[rownames(seurat_object@meta.data),]$doublets=seurat_object@meta.data[,ncol(seurat_object@meta.data)]
cells_metadata[rownames(seurat_object@meta.data),]$pANN=seurat_object@meta.data[,ncol(seurat_object@meta.data)-1]
cells_metadata[rownames(seurat_object@meta.data),]$nCount_RNA=seurat_object@meta.data$nCount_RNA
cells_metadata[rownames(seurat_object@meta.data),]$nFeature_RNA=seurat_object@meta.data$nFeature_RNA

#seurat_object <-seurat_object[,seurat_object@meta.data[,ncol(seurat_object@meta.data)]=="Singlet"]
#rm(annotations,homotypic.prop,nExp_poi,nExp_poi.adj,pK )

#1.7 Add seurat object
#num_layers <- num_layers +1
#layers[[num_layers]] <- "exp1"
#seurat_list[[num_layers]] <- CreateSeuratObject(counts = seurat_object@assays$RNA$counts, meta.data = seurat_object@meta.data[,c("individualID","libraryBatch","sequencingBatch")])
#seurat_list[[num_layers]]$dataset <- "exp1"
rm(seurat_object,features,cells, exp1_sparse_matrix_file,exp1_features_file,exp1_cells_file)
rm(biospecimen_metadata)




#################################################################################
#2. Load  cells from the rest of the datasets
#################################################################################

#2.1 Read counts
for (i in 1:length(files)){
  name <- gsub('.matrix.mtx.gz','', files[i]) #Basename of file
  str<-sprintf("[%d/%d] %s", i, length(files), name)
  print(str)
  
  
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
  batch_cells <-rownames(filter(cells_metadata, libraryBatch==name & qc == "Pass"))
  if ( length(batch_cells)<2 )
  {
    print(" Data not found")
    rm(cts)
    next
  }
  
  cts<-cts[,batch_cells]
  

  
  #2.1.3 Create seurat object 
  seurat_object  <- CreateSeuratObject(counts = cts, meta.data = cbind(cells_metadata[batch_cells,c("individualID","libraryBatch","sequencingBatch")]))
  rm(batch_cells)
  
  #2.1.4 Compute percentage MT genes and remove low quality genes
  seurat_object <- PercentageFeatureSet(seurat_object, pattern = "^MT-", col.name = "percent.mt")
  #VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  
  #df<-seurat_object@meta.data[,c("nFeature_RNA", "nCount_RNA", "percent.mt")]
  #df$x <- 1
  #png(file=paste0(log_folder,'/',cell_type,'/qc_all/',name,'.png'), width=600, height=600)
  #p1<- ggplot(df, aes( x=x, y = nFeature_RNA)) +
  #  ggdist::stat_halfeye(adjust = .5,  width = .6, .width = 0,  justification = -.3,  point_colour = NA) + 
  #  geom_boxplot( width = .25,  outlier.shape = NA ) +
  #  geom_point( size = 1.5, alpha = .2, position = position_jitter(seed = 1, width = .1 ) ) +
  #  theme_bw()
  #p2<- ggplot(df, aes( x=x, y = nCount_RNA)) +
  #  ggdist::stat_halfeye(adjust = .5,  width = .6, .width = 0,  justification = -.3,  point_colour = NA) + 
  #  geom_boxplot( width = .25,  outlier.shape = NA ) +
  #  geom_point( size = 1.5, alpha = .2, position = position_jitter(seed = 1, width = .1 ) ) +
  #  theme_bw()
  #p3<- ggplot(df, aes( x=x, y = percent.mt)) +
  #  ggdist::stat_halfeye(adjust = .5,  width = .6, .width = 0,  justification = -.3,  point_colour = NA) + 
  #  geom_boxplot( width = .25,  outlier.shape = NA ) +
  #  geom_point( size = 1.5, alpha = .2, position = position_jitter(seed = 1, width = .1 ) ) +
  #  theme_bw()
  #p4<-ggplot(df, aes(x = nCount_RNA, y = nFeature_RNA)) + 
  #  geom_point( size = 1.5, alpha = .2) + 
  #  geom_smooth(method = 'lm') +
  #  theme_bw()
  #p<- p1 + p2 + p3 + p4
  #print(p)
  #dev.off()
  #rm(df,p,p1,p2,p3,p4,p5,p6)
  
  
  #seurat_object <- subset(seurat_object, subset = nCount_RNA > 400 & nFeature_RNA > 253  & percent.mt < 2)
  #seurat_object <- subset(seurat_object,  percent.mt < percent_mt_thr) ####################
  
  #cells_metadata[colnames(seurat_object$percent.mt>=percent_mt_thr),]$qc <- "MT"
  
  
  if ( ncol(seurat_object)<pca_npcs )
  {
    print(" Number of cells < pca_npcs")
    rm(seurat_object)
    next
  }
  
  
  #2.1.6 . Remove doubles (DoubletFinder)
  #2.1.6.1 Basic preprocessing
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object, nfeatures = num_features) 
  seurat_object <- ScaleData(seurat_object) 
  seurat_object <- RunPCA(seurat_object, npcs=pca_npcs)
  seurat_object <- FindNeighbors(seurat_object, dims = 1:num_dimensions_umap) #construction of k-NN neighbor graph
  seurat_object <- FindClusters(seurat_object,resolution=1.5)
  seurat_object <- RunUMAP(seurat_object, dims = 1:num_dimensions_umap)
  
  #2.1.6.2 pK Identification (no ground-truth data available to compute it)
  sweep.res.list_cells <- paramSweep(seurat_object, PCs = 1:pca_npcs,sct = FALSE) # Performs steps 1,2,3 4. Provides list of the proportion of artificial neighbours for combination of pN and pK
  sweep.stats <- summarizeSweep(sweep.res.list_cells, GT = FALSE) # Summarize de results
  bcmvn <- find.pK(sweep.stats) # metric to compute best pk (highest value)
  #ggplot(bcmvn, aes(pK, BCmetric, group = 1)) + geom_point() + geom_line() # plot results
  
  pK <- bcmvn %>% # select the pK that corresponds to max bcmvn to optimize doublet detection
    filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  pK <- as.numeric(as.character(pK[[1]]))
  rm(bcmvn,sweep.res.list_cells,sweep.stats)
  
  #2.1.6.3 Homotypic Doublet Proportion Estimate
  annotations <- seurat_object@meta.data$seurat_clusters # Get the cluster of each cell
  homotypic.prop <- modelHomotypic(annotations)         #compute perfect of homotopic doublets
  nExp_poi <- round(multiplet_rate*nrow(seurat_object@meta.data))  #Expected number of doublets. Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  #2.1.6.4 Run doublet pfinder with estimated parameters
  #seurat_object <- doubletFinder(seurat_object,
  #                               PCs = 1:pca_npcs,
  #                               pN = 0.25,
  #                               pK = pK,
  #                               nExp = nExp_poi.adj,
  #                               reuse.pANN = FALSE, sct = FALSE)
  
  seurat_object <- doubletFinder(seu = seurat_object, 
                                 PCs = 1:pca_npcs,
                                 pN = 0.25,
                                 pK = pK, # the neighborhood size used to compute the number of artificial nearest neighbours
                                 nExp = nExp_poi.adj) # number of expected real doublets
  
  
  
  
  #2.1.6.5 If a cluster > 70% doublets all doubles
  df<-seurat_object@meta.data[c(which(colnames(seurat_object@meta.data)=="seurat_clusters"),ncol(seurat_object@meta.data))]
  df<-as.data.frame.matrix(table(df))
  df_relative <- sweep(x = df, MARGIN = 1, STATS = rowSums(df), FUN = "/")
  seurat_object@meta.data[seurat_object$seurat_clusters==rownames(df_relative[df_relative[,"Doublet"]>=0.7,]),ncol(seurat_object@meta.data)] <- "Doublet" 
  rm(df,df_relative)
  
  
  #2.1.6.6 Compute ncount, nfeature, andpercent.mt distance of each cell in each cluster
  clusts <- levels(seurat_object@meta.data$seurat_clusters)
  dist_outliers_nCount_RNA <- cells_metadata[rownames(seurat_object@meta.data),]$dist_outliers_nCount_RNA
  dist_outliers_nFeature_RNA <- cells_metadata[rownames(seurat_object@meta.data),]$dist_outliers_nFeature_RNA
  dist_outliers_percent.mt  <- cells_metadata[rownames(seurat_object@meta.data),]$dist_outliers_percent.mt
  for(j in 1:length(clusts) )
  {
    indxc <- seurat_object@meta.data$seurat_clusters==clusts[j]
    
    #nCountRNA
    x <- seurat_object@meta.data[indxc,]$nCount_RNA
    dist_outliers_nCount_RNA[indxc] <- (  x - quantile(x,.75) ) / IQR(x)
    
    #nFeatureRNA
    x <- seurat_object@meta.data[indxc,]$nFeature_RNA
    dist_outliers_nFeature_RNA[indxc] <- (  x - quantile(x,.75) ) / IQR(x)
    
    
    #percent.mt
    x <- seurat_object@meta.data[indxc,]$percent.mt
    iqr<-IQR(x)
    if (iqr < 2 )
    {
      iqr <-2
    }
    dist_outliers_percent.mt[indxc] <- (  x - quantile(x,.75) ) / iqr
    
    
  }
  
  cells_metadata[rownames(seurat_object@meta.data),]$dist_outliers_nCount_RNA <- dist_outliers_nCount_RNA 
  cells_metadata[rownames(seurat_object@meta.data),]$dist_outliers_nFeature_RNA <- dist_outliers_nFeature_RNA
  cells_metadata[rownames(seurat_object@meta.data),]$dist_outliers_percent.mt <- dist_outliers_percent.mt

  rm(clusts,indxc,x,dist_outliers_nCount_RNA,dist_outliers_nFeature_RNA,dist_outliers_percent.mt,i)
  
  
  
  # visualize doublets
  png(file=paste0(log_folder,'/',cell_type,'/qc/',name,'.png'), width=2000, height=800)
  df<- seurat_object@meta.data[,c("nCount_RNA","nFeature_RNA","percent.mt","seurat_clusters")]
  df$outliers_nCount_RNA<- cells_metadata[rownames(seurat_object@meta.data),]$dist_outliers_nCount_RNA> 5
  df$outliers_nFeature_RNA<- cells_metadata[rownames(seurat_object@meta.data),]$dist_outliers_nFeature_RNA> 5
  df$outliers_percent.mt<- cells_metadata[rownames(seurat_object@meta.data),]$dist_outliers_percent.mt > 3
  
  
  p1<-DimPlot(seurat_object, reduction = 'umap', group.by = 'seurat_clusters', label = TRUE)
  p2<-DimPlot(seurat_object, reduction = 'umap', group.by = colnames(seurat_object@meta.data)[ncol(seurat_object@meta.data)]) + plot_annotation(title = name)
  p3<-FeaturePlot(seurat_object, features = colnames(seurat_object@meta.data)[ncol(seurat_object@meta.data)-1])

  p4<- ggplot(df, aes( x=seurat_clusters, y = nCount_RNA)) +
    ggdist::stat_halfeye(adjust = .5,  width = .6, .width = 0,  justification = -.3,  point_colour = NA) + 
    geom_boxplot( width = .25,  outlier.shape = NA ) +
    geom_point( data =df, aes( x=seurat_clusters, y = nCount_RNA, color = outliers_nCount_RNA),size = 1.5, alpha = .2, position = position_jitter(seed = 1, width = .1 ) ) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
    geom_hline(yintercept=nCount_RNA_thr, linetype="dashed", color = "red") +
    theme_bw()
  p5<-  ggplot(df, aes( x=seurat_clusters, y = nFeature_RNA)) +
    ggdist::stat_halfeye(adjust = .5,  width = .6, .width = 0,  justification = -.3,  point_colour = NA) + 
    geom_boxplot( width = .25,  outlier.shape = NA ) +
    geom_point( data =df, aes( x=seurat_clusters, y = nFeature_RNA, color = outliers_nFeature_RNA),size = 1.5, alpha = .2, position = position_jitter(seed = 1, width = .1 ) ) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
    geom_hline(yintercept=nFeature_RNA_thr, linetype="dashed", color = "red") +
    theme_bw()
  p6<- ggplot(df, aes( x=seurat_clusters, y = percent.mt)) +
    ggdist::stat_halfeye(adjust = .5,  width = .6, .width = 0,  justification = -.3,  point_colour = NA) + 
    geom_boxplot( width = .25,  outlier.shape = NA ) +
    geom_point( data =df, aes( x=seurat_clusters, y = percent.mt, color = outliers_percent.mt),size = 1.5, alpha = .2, position = position_jitter(seed = 1, width = .1 ) ) +
    scale_color_manual(values = c("FALSE" = "black", "TRUE" = "red")) +
    geom_hline(yintercept=percent.mt_up_thr, linetype="dashed", color = "red") +
    theme_bw()
  
  p7<-DotPlot(seurat_object, features = features_qc) +
    geom_vline(xintercept = features_qc_sep, linetype="dashed") +
    RotatedAxis()
  
  p8<-ggplot(df, aes(x = nCount_RNA, y = nFeature_RNA)) + 
    geom_point( size = 1.5, alpha = .2) + 
    geom_smooth(method = 'lm') +
    theme_bw()
  p <- p1 + p3 + p2 +p7 + p4 +p5 + p6 +p8  + plot_layout(ncol = 4)
  print(p)
  dev.off()
  rm(df,p,p1,p2,p3,p4,p5,p6,p7,p8)
  

  
  #table(seurat_object@meta.data[,ncol(seurat_object@meta.data)]) # number of singlets and doublets
  
  #2.1.6.7 Remove doublets
  cells_metadata[rownames(seurat_object@meta.data),]$doublets=seurat_object@meta.data[,ncol(seurat_object@meta.data)]
  cells_metadata[rownames(seurat_object@meta.data),]$pANN=seurat_object@meta.data[,ncol(seurat_object@meta.data)-1]
  cells_metadata[rownames(seurat_object@meta.data),]$nCount_RNA=seurat_object@meta.data$nCount_RNA
  cells_metadata[rownames(seurat_object@meta.data),]$nFeature_RNA=seurat_object@meta.data$nFeature_RNA
  

  #seurat_object <-seurat_object[,seurat_object@meta.data[,ncol(seurat_object@meta.data)]=="Singlet"]
  rm(annotations,homotypic.prop,nExp_poi,nExp_poi.adj,pK )
  
  #2.1.8 Create seurat object
  #num_layers <- num_layers + 1
  #seurat_list[[num_layers]]  <- CreateSeuratObject(counts = seurat_object@assays$RNA$counts, meta.data = seurat_object@meta.data[,c("individualID","libraryBatch","sequencingBatch")])
  #seurat_list[[num_layers]]$dataset <- name
  #layers[[num_layers]] <-  name
  
  rm(cts,seurat_object,str,name)
  
}

write.csv(cells_metadata, file= paste0(log_folder, "/",cell_type,"/qc.csv" ))



  



