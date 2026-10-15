
#0. clear workspace and Set working directory
rm(list = ls())
setwd('/Users/ibaidiez/Desktop/projects/ROSMAP/analysis/NatComm')



#Load libraries
library(Seurat)
library(tidyverse)
library(data.table)
library(Matrix)


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
dir.create(log_folder, showWarnings = FALSE)



########################## REFERENCE DATASET ################################# 

#################################################################################
#1. Load reference dataset
#################################################################################

#1.1 Load Metadata
features <- read.csv(exp1_features_file, header=T,stringsAsFactors = F)
cells <- read.csv(exp1_cells_file, header=T,stringsAsFactors = F)
rownames(cells)<-cells$cell_name
biospecimen_metadata <- read.csv(biospecimen_file, header=T,stringsAsFactors = F)
biospecimen_metadata <- biospecimen_metadata[biospecimen_metadata$assay == "scrnaSeq",]
rownames(biospecimen_metadata)<- biospecimen_metadata$specimenID


#1.2 Read sparce matrix and convert to count matrix
counts <- fread(exp1_sparse_matrix_file,skip=1)
i=as.vector(counts$V1)
j=as.vector(counts$V2)
x=as.vector(counts$V3)
rm(counts)

counts <- sparseMatrix(i, j, x = x)
rownames(counts)=features$x
colnames(counts)=cells$cell_name
rm(i,j,x)

# #1.3 Remove cells based on posterior processins steps
#counts <- counts[,!colnames(counts) %in% cells_rm]

#1.4 Create seurat object 
cells <- cells[colnames(counts),]
individualIDs<-biospecimen_metadata[cells$specimenID,"individualID"]
seurat_ref <- CreateSeuratObject(counts = counts, meta.data = cbind(cells[,c("specimenID","Cell.Type","sub.cluster")],individualIDs))
rm(counts,features,cells, individualIDs, exp1_sparse_matrix_file,exp1_features_file,exp1_cells_file)



################################
#2. Standard preprocessing
################################

#2.1 Data normalization and scaling.
seurat_ref <- NormalizeData(seurat_ref) #Normalize feture expression
seurat_ref <- FindVariableFeatures(seurat_ref,selection.method = "vst", nfeatures = 3000) #calculate a subset of features that exhibit high cell-to-cell variation in the dataset
seurat_ref <- ScaleData(seurat_ref) #Shifts the expression of each gene, so that the mean expression across cells is 0 and the variance across cells is 1.

#Remove counts
#str(seurat_cts@assays)
seurat_ref<-DietSeurat(seurat_ref,counts = FALSE,data = TRUE, scale.data = TRUE)

#2.2 Dimensionality reduction
pca_npcs<-30
seurat_ref <- RunPCA(seurat_ref, npcs=pca_npcs)

# png(file=paste0(log_folder,'/Num_components.png'), width=800, height=500)
# ElbowPlot(seurat_ref,ndims=pca_npcs)
# dev.off()

#2.3 Clustering of cells
num_components <- pca_npcs
seurat_ref <- FindNeighbors(seurat_ref, dims = 1:num_components) #construction of k-NN neighbor graph
seurat_ref <- FindClusters(seurat_ref,resolution=0.2,algorithm = 1)
seurat_ref <- RunUMAP(seurat_ref, dims = 1:num_components, return.model = TRUE)


#Plot UMAP
#DimPlot(seurat_ref, reduction = "umap")
#DimPlot(seurat_ref, reduction = "umap",label = TRUE, pt.size = 0.5)

png(file=paste0(log_folder,'/clusters.png'), width=600, height=500)
p <- DimPlot(seurat_ref, reduction = "umap",label = TRUE, pt.size = 0.5)
print(p)
dev.off()



#2.4 Use known makers to identify clusters
# Dot plots - the size of the dot corresponds to the percentage of cells expressing the
# feature in each cluster. The color represents the average expression level

#write.csv(rownames(seurat_ref@assays$RNA@data), file="genes.csv")

png(file=paste0(log_folder,'/clusters_markers.png'), width=800, height=500)
features <-c("SLC17A7","GAD1","GRIK1","PLP1","MBP","VCAN","PDGFRA","GFAP","SLC1A2","CD74","P2RY12","FLT1","CLDN5" )
#features <-c("SLC17A7","GAD1","GRIK1","PLP1","MBP","VCAN","PDGFRA","GFAP","SLC1A2","CD74","C3","P2RY12","FLT1","CLDN5" )
p <- DotPlot(seurat_ref, features = features) + RotatedAxis() +
  geom_vline(xintercept = c(1.5,3.5,5.5,7.5,9.5,11.5), linetype="dashed") 

print(p)
dev.off()


#Plor reordered
#my_levels <- c(0,5,14,8,4,9,13,15,3,6,11,1,7,2,12,10)
#factor(Idents(seurat_ref), levels= my_levels)
#Idents(seurat_ref) <- factor(Idents(seurat_ref), levels= my_levels)
#DotPlot(seurat_ref, features = features, ) + RotatedAxis() +
#  geom_vline(xintercept = c(1.5,3.5,5.5,7.5,9.5,11.5), linetype="dashed") 

#my_levels <- c(0:15)
#factor(Idents(seurat_ref), levels= my_levels)
#Idents(seurat_ref) <- factor(Idents(seurat_ref), levels= my_levels)


#Plot genes in the UMAP
#FeaturePlot(seurat_ref, features = features, raster=FALSE)
#FeaturePlot(seurat_ref, features = c("SLC17A7","GAD1","SLC1A2","C3"), raster=FALSE)

#Plot UMAP with name of cells
#new.cluster.ids <- c( "ExcNeur0","Oligodendrocytes", "Astrocytes", "InhNeur3", "ExcNeur4", "ExcNeur5", "InhNeur6", "OPCs",
#                      "ExcNeur8", "ExcNeur9", "Vascular Nitch", "InhNeur11", "Microglia", "ExcNeur13", "ExcNeur14", "ExcNeur15" )
#names(new.cluster.ids) <- levels(seurat_ref)
#seurat_ref <- RenameIdents(seurat_ref, new.cluster.ids)
#DimPlot(seurat_ref, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
#rm(new.cluster.ids )

#2.5 Finding differentially expressed features (cluster biomarkers)
#markers_cluster <- FindMarkers(seurat_cts, ident.1 = "ExtNeur4", min.pct = 0.5)
#head(markers_cluster)


# find markers for every cluster compared to all remaining cells, report only the positive  ones
# markers <- FindAllMarkers(seurat_cts, logfc.threshold = 0.25,  min.pct = 0.1, only.pos = TRUE, test.use = 'DESeq2', slot = 'counts')

#2.6 Define 7 main cell types
#Create Clusters7
clusters7cells <- vector(length = length(seurat_ref@meta.data$seurat_clusters))
clusters7cells[seurat_ref@meta.data$seurat_clusters==0]= "Exc" #Excitatory neurons
clusters7cells[seurat_ref@meta.data$seurat_clusters==5]= "Exc" #Excitatory neurons
clusters7cells[seurat_ref@meta.data$seurat_clusters==8]= "Exc" #Excitatory neurons
clusters7cells[seurat_ref@meta.data$seurat_clusters==14]= "Exc" #Excitatory neurons
clusters7cells[seurat_ref@meta.data$seurat_clusters==4]= "Exc" #Excitatory neurons
clusters7cells[seurat_ref@meta.data$seurat_clusters==9]= "Exc" #Excitatory neurons
clusters7cells[seurat_ref@meta.data$seurat_clusters==13]= "Exc" #Excitatory neurons
clusters7cells[seurat_ref@meta.data$seurat_clusters==15]= "Exc" #Excitatory neurons
clusters7cells[seurat_ref@meta.data$seurat_clusters==3]= "Inh" #Inhibitory neurons
clusters7cells[seurat_ref@meta.data$seurat_clusters==6]= "Inh" #Inhibitory neurons
clusters7cells[seurat_ref@meta.data$seurat_clusters==11]= "Inh" #Inhibitory neurons
clusters7cells[seurat_ref@meta.data$seurat_clusters==1]= "Oligo" #Oligodendrocytes
clusters7cells[seurat_ref@meta.data$seurat_clusters==7]= "OPCs" #OPCs
clusters7cells[seurat_ref@meta.data$seurat_clusters==2]= "Astrc" #Astrocytes
clusters7cells[seurat_ref@meta.data$seurat_clusters==12]= "Micr" #Microglia
clusters7cells[seurat_ref@meta.data$seurat_clusters==10]= "VN" #Vascular Nitch
clusters7cells <-factor(clusters7cells)
seurat_ref@meta.data <- cbind(seurat_ref@meta.data,clusters7cells)
rm(clusters7cells)

png(file=paste0(log_folder,'/7celltypes.png'), width=600, height=500)
p <- DimPlot(seurat_ref, reduction = "umap", group.by = "clusters7cells", label = TRUE, pt.size = 0.5) + NoLegend()
print(p)
dev.off()

#Compare obtained cell types and provided cell types
cell_types<- unique(seurat_ref$clusters7cells)
lbls_cell <- integer(ncol(seurat_ref))
names(lbls_cell) <- colnames(seurat_ref)
for (i in 1:length(cell_types)){
  celltype<- levels(cell_types)[i]
  print(celltype)
  ctps<-seurat_ref[,seurat_ref$clusters7cells==celltype]$Cell.Type
  t<-table(ctps)
  print(t)
  
  max_name<-names(which(t==max(t)))
  lbls_cell[names(ctps[ctps==max_name])]=1 #correct
  lbls_cell[names(ctps[ctps!=max_name])]=2 #incorrect
  rm(t,max_name)
}
rm(celltype,i,my_levels,ctps)


df <- data.frame(x = seurat_ref@reductions$umap@cell.embeddings[,c("umap_1")], y = seurat_ref@reductions$umap@cell.embeddings[,c("umap_2")],
                 classification = lbls_cell )
ggplot() + 
  geom_point(aes(x = x, y = y), data = df[df$classification==1,], size=0.25,alpha=0.1, color= rgb(0.75,0.75,0.75) ) + # correctly classified
  geom_point(aes(x = x, y = y), data = df[df$classification==2,], size=0.25,alpha=0.9, color= rgb(228,26,28,maxColorValue = 255) ) + # incorrectly classified
  xlab('UMAP1') + ylab('UMAP2') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
rm(df,lbls_cell)


# Estimation: 172,641 correct [99.99% correct]
# Astrc: 29,486 astrocytes [100% correct]
# Exc: 74,992 Exc, 6 Inh [99.99% correct]
# Inh: 24,929 Inh, 5 Exc [99.99% correct]
# Micr: 2,836 Micr, 1 VLMC [99.96% correct]
# Oligo: 29,543 Oligo, 1 Micr [99.99% correct]
# OPCs: 7,493 OPCs [100% correct]
# VN: 2,296 Endothelial 923 pericytes 143 VLMC  2 Exc 3 Inh [99.85% correct]



#2.7 Store seurat object for subsequent analyses

#Write cell types info
dir.create(paste0(out_folder, "/cells_types" ), showWarnings = FALSE)
dir.create(paste0(out_folder, "/seurat" ), showWarnings = FALSE)
seurat_ref@meta.data$libraryBatch <- "exp1"
seurat_ref@meta.data$sequencingBatch <- "exp1"
write.csv(cbind(seurat_ref@reductions$umap@cell.embeddings,seurat_ref@meta.data[,c("individualIDs","libraryBatch","sequencingBatch","clusters7cells","seurat_clusters")]), file= paste0(out_folder, "/cells_types/exp1.csv" ))
SaveSeuratRds(object = seurat_ref,file = paste0(out_folder, "/seurat/ref_exp1.Rds" ) )



#Write cell type proportions
thr_prediction_score<-0.75
df_comp <- as.data.frame.matrix(table(seurat_ref$individualIDs, seurat_ref$clusters7cells))
df_comp_relative <- sweep(x = df_comp, MARGIN = 1, STATS = rowSums(df_comp), FUN = "/")
colnames(df_comp_relative)<-paste0("prop_",colnames(df_comp_relative))
df<-cbind(df_comp,df_comp_relative)
df$library <- "exp1"
df$individualID <- rownames(df)
write.csv(df, file= paste0(out_folder, "/individual_cell_composition_ref.csv" ))
rm(df_comp,df_comp_relative,df)




########################## PROJECT REST OF DATA  ################################# 
dir.create(paste0(out_folder,"/cells_types"), showWarnings = FALSE)


#################################################################################
#0. Prepare metadata experiment 2
#################################################################################
demultiplex_file <- '../../data/Metadata/ROSMAP_snRNAseq_demultiplexed_ID_mapping.csv';
scrna_metadata_file <- '../../data/Metadata/ROSMAP_assay_scrnaSeq_metadata.csv';

#1.1 Load  cell specific metadata
if  (file.exists(paste0(out_folder, "/data/cell_metadata.csv" ))==FALSE )
{
  metad.demultiplexID <- read.csv(demultiplex_file, header=T,stringsAsFactors = F)
  metad.demultiplexID$barcode <- apply( metad.demultiplexID[ , 1:2 ] , 1 , paste , collapse = "_" )
  metad.demultiplexID<-metad.demultiplexID[, c("barcode","libraryBatch","individualID")]
  specimenID<-apply( metad.demultiplexID[ ,c("libraryBatch","individualID")] , 1 , paste , collapse = "_" )
  #rownames(metad.demultiplexID)<-metad.demultiplexID$barcode
  
  metad.biospecimen <- read.csv(biospecimen_file, header=T,stringsAsFactors = F)
  metad.biospecimen  <- metad.biospecimen [metad.biospecimen $assay == "scrnaSeq",]
  rownames(metad.biospecimen)<-metad.biospecimen$specimenID
  metad.demultiplexID$exclude <- metad.biospecimen[specimenID,"exclude"]
  rm(metad.biospecimen)
  
  metad.scrna <- read.csv(scrna_metadata_file, header=T,stringsAsFactors = F)
  metad.scrna <- metad.scrna[metad.scrna[,"platformLocation"]=="Broad",]
  rownames(metad.scrna)<-metad.scrna$specimenID
  metad.demultiplexID$sequencingBatch <- metad.scrna[specimenID,"sequencingBatch"]
  rm(metad.scrna,specimenID)
  
  #write cell metadata
  dir.create(paste0(out_folder, "/data" ) , showWarnings = FALSE)
  write.csv(metad.demultiplexID, file= paste0(out_folder, "/data/cell_metadata.csv" ),row.names=FALSE )
}else
{
  metad.demultiplexID <- read.csv(paste0(out_folder, "/data/cell_metadata.csv" ), header=T,stringsAsFactors = F)
  rownames(metad.demultiplexID)<-metad.demultiplexID$barcode
}

#1.2 Split the data into libraries
dir.create(paste0(out_folder, "/cell_metadata" ) , showWarnings = FALSE)
metad.demultiplexID$libraryBatch <- as.factor(metad.demultiplexID$libraryBatch)
library_batches<-unique(metad.demultiplexID$libraryBatch)
for (i in 1:length(library_batches)){
  name <- library_batches[i]
  str<-sprintf("[%d/%d] %s", i, length(library_batches), name)
  print(str)
  metads<-subset(metad.demultiplexID,  libraryBatch == name)
  write.csv(metads, file= paste0(out_folder, "/cell_metadata/",name,".csv" ),row.names=FALSE )
}

rm(i,demultiplex_file,scrna_metadata_file,metad.demultiplexID,metads,name,str,library_batches)



#################################################################################
#2. Project data
#################################################################################
#seurat_ref <- LoadSeuratRds( paste0(out_folder, "/seurat/ref_exp1.Rds" ))

#Each file is a different batch
files <- list.files(path = exp2_data_path, pattern = "*.mtx.gz")

dir.create(paste0(log_folder,'/mapping7cell'), showWarnings = FALSE)
dir.create(paste0(log_folder,'/mapping'), showWarnings = FALSE)


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
  
  
  #2.1.2  Load metadata and remove cells not assigned to participants
  metadata_library <- read.csv(paste0(out_folder, "/cell_metadata/",name,".csv" ), header=T,stringsAsFactors = F)
  rownames(metadata_library)<-metadata_library$barcode
  
  #Remove cells not assigned to subjects
  cts<-cts[,is.na(metadata_library[colnames(cts),"libraryBatch"])==FALSE]
  
  #Remove cells marked as exclude
  cts<-cts[,metadata_library[colnames(cts),"exclude"]==FALSE]
  
 if  (ncol(cts)==0)
 {
   rm(cts,metadata_library)
   next
 }
  

  #2.1.4 Create seurat object
  seurat_cts <- CreateSeuratObject(counts = cts, meta.data = metadata_library[colnames(cts),c("libraryBatch","individualID","sequencingBatch")])
  seurat_cts$dataset <- name
  rm(cts,metadata_library)
  
  #2.1.5  Basic preprocessing
  seurat_cts <- NormalizeData(seurat_cts)
  
  #2.1.6  Map dataset to reference
  # Find a set of anchors between a reference and query object. 
  anchors <- FindTransferAnchors(reference = seurat_ref, query = seurat_cts, reduction= "pcaproject", dims = 1:30, reference.reduction = "pca")

  #Transfer cell type from reference to dataset 
  predictions <- TransferData(anchorset = anchors, refdata = seurat_ref$clusters7cells, dims = 1:30)
  colnames(predictions)<-paste0("clust7.",colnames(predictions))
  seurat_cts <- AddMetaData(seurat_cts, metadata = predictions)
  rm(predictions)
  predictions <- TransferData(anchorset = anchors, refdata = seurat_ref$seurat_clusters, dims = 1:30)
  colnames(predictions)<-paste0("clustSeur.",colnames(predictions))
  seurat_cts <- AddMetaData(seurat_cts, metadata = predictions)
  rm(predictions)
  
  
  #2.1.7 UMAP projection
  # Perform dataset integration using a pre-computed Anchorset of specified low dimensional representations.
  seurat_cts <- IntegrateEmbeddings(anchorset = anchors, reference = seurat_ref, query = seurat_cts, new.reduction.name = "ref.pca")
  seurat_cts <- ProjectUMAP(query = seurat_cts, query.reduction = "ref.pca", reference = seurat_ref, reference.reduction = "pca", reduction.model = "umap")
  
 #2.1.8 Use KNN to classify cell types and compute a distance score
  ref_embd<- seurat_ref@reductions$umap@cell.embeddings
  colnames(ref_embd) <- c("x","y")
  cts_embd<- seurat_cts@reductions$ref.umap@cell.embeddings
  colnames(cts_embd) <- c("x","y")
  lbls_clust7<-factor(seurat_ref$clusters7cells)
  lbls_clustSeurat<-factor(seurat_ref$seurat_clusters)
  
  k<-100
  cts_lbl<-c()
  cts_lbl7<-c()
  cts_kdist<-c()
  for (c in 1:nrow(cts_embd)){
    #Compute distance
    d<-sqrt(rowSums(sweep(ref_embd,2, cts_embd[c,],FUN = "-")^2))
    #KNN to assign clusgter
    kd<-sort(d)[1:k]
    cts_kdist[c]<-mean(kd)
    cts_lbl7[c]<-names(which.max(table(lbls_clust7[names(kd)])))
    cts_lbl[c]<-names(which.max(table(lbls_clustSeurat[names(kd)])))
    
  }
  seurat_cts <- AddMetaData(object = seurat_cts, metadata = cts_lbl, col.name = 'knn_clust')
  seurat_cts <- AddMetaData(object = seurat_cts, metadata = cts_lbl7, col.name = 'knn_clust7')
  seurat_cts <- AddMetaData(object = seurat_cts, metadata = cts_kdist, col.name = 'knn_dist')
  rm(kd,d,cts_kdist,cts_lbl)
  
  #Save png
  png(file=paste0(log_folder,'/mapping7cell/',name,'_mapping_7celltype.png'), width=1500, height=500)
  p1 <- DimPlot(seurat_ref, reduction = "umap", group.by = "clusters7cells", label = TRUE, label.size = 3,
                repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
  p2 <- DimPlot(seurat_cts, reduction = "ref.umap", group.by = "clust7.predicted.id", label = TRUE,
                label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
  p3 <- DimPlot(seurat_cts, reduction = "ref.umap", group.by = "knn_clust7", label = TRUE,
                label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
  p <- p1 + p2 + p3
  print(p)
  dev.off()
  
  png(file=paste0(log_folder,'/mapping/',name,'_mapping.png'), width=1500, height=500)
  p1 <- DimPlot(seurat_ref, reduction = "umap", group.by = "seurat_clusters", label = TRUE, label.size = 3,
                repel = TRUE) + NoLegend() + ggtitle("Reference annotations")
  p2 <- DimPlot(seurat_cts, reduction = "ref.umap", group.by = "clustSeur.predicted.id", label = TRUE,
                label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
  p3 <- DimPlot(seurat_cts, reduction = "ref.umap", group.by = "knn_clust", label = TRUE,
                label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
  p <- p1 + p2 + p3
  print(p)
  dev.off()
  
  #2.1.8 Save cell info
  #num_layers <- num_layers + 1
  write.csv(cbind(seurat_cts@reductions$ref.umap@cell.embeddings,seurat_cts@meta.data[,c("individualID","libraryBatch","sequencingBatch","clust7.predicted.id","clust7.prediction.score.Exc", "clust7.prediction.score.Inh","clust7.prediction.score.VN","clust7.prediction.score.OPCs","clust7.prediction.score.Astrc","clust7.prediction.score.Oligo","clust7.prediction.score.Micr","clust7.prediction.score.max","clustSeur.predicted.id","clustSeur.prediction.score.max","knn_clust","knn_clust7","knn_dist")]), file= paste0(out_folder, "/cells_types/",name,".csv" ))
  
  rm(seurat_cts,p,p1,p2,p3,anchors,ref_embd,cts_embd,lbls_clust7,lbls_clustSeurat)
 
}


#2.2 Compare our cell types and web
cell_types<-c()
for (i in 1:length(files)){
  name <- gsub('.matrix.mtx.gz','', files[i]) #Basename of file
  str<-sprintf("[%d/%d] %s", i, length(files), name)
  print(str)
  if  (file.exists(paste0(out_folder, "/cells_types/",name,".csv"))==FALSE ) {
    next
  }
  
  celltps <- read.csv( paste0(out_folder, "/cells_types/",name,".csv") , header=T,stringsAsFactors = F)
  cell_types<-rbind(cell_types,celltps)
  rm(celltps)
}
rownames(cell_types) <- cell_types$X


cell_anot_file <- "../../data/Gene_Expression/snRNAseq/cell-annotation.csv"
cell_anot <- read.csv(cell_anot_file, header=T,stringsAsFactors = F)
rownames(cell_anot) <- cell_anot$barcode

#uctps <- unique(cell_types$predicted.id)
uctps <- unique(cell_types$knn_clust7)
lbls_cell <- integer(nrow(cell_types))
names(lbls_cell) <- rownames(cell_types)
for (i in 1:length(uctps)){
  ct<- uctps[i]
  print(ct)
  #a<-intersect(cell_types[cell_types$clust7.predicted.id==ct,"X"],rownames(cell_anot))
  a<-intersect(cell_types[cell_types$knn_clust7==ct,"X"],rownames(cell_anot))
  b<-table(cell_anot[a,]$grouping.by)
  b<-100*b/sum(b)
  print(table(cell_anot[a,]$grouping.by))
  
  #cell_types[rownames(cell_anot[cell_anot[a,]$grouping.by!="Microglia",]),]
  
  print(sprintf("%.2f",b))
  lbls_cell[a]=1 # Correct
  
  lbls_cell[a[cell_anot[a,]$grouping.by!=names(which(b==max(b)))]]=2 # incorrect
        

  
  rm(a,b)
}

# unique(cell_anot$cell.type)
rm(cell_anot_file,cell_anot,ct,i,str,uctps,name)


#2.3 Plot integrated data UMAP
library(data.table)
df <- data.frame(x = cell_types$refUMAP_1, y = cell_types$refUMAP_2,
                 Species = cell_types$knn_clust7, classification = lbls_cell )

centroids <- aggregate(cbind(x,y)~Species,df,mean)
ggplot(df, aes(x, y, colour = Species)) +
  geom_point(size=0.25,alpha=0.1) +
  geom_text(data = centroids, aes(label = Species),color = "black") +
  xlab('UMAP1') + ylab('UMAP2') + 
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


#Plot incorrectly labeled
ggplot() + 
  geom_point(aes(x = x, y = y), data = df[df$classification==1,], size=0.25,alpha=0.1, color= rgb(0.75,0.75,0.75) ) + # correctly classified
  geom_point(aes(x = x, y = y), data = df[df$classification==0,], size=0.25,alpha=0.1, color= rgb(55,126,184,maxColorValue = 255) ) + # missing cells
  geom_point(aes(x = x, y = y), data = df[df$classification==2,], size=0.25,alpha=0.1, color= rgb(228,26,28,maxColorValue = 255) ) + # incorrectly classified
  xlab('UMAP1') + ylab('UMAP2') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


#sum(df$classification==1)/sum(df$classification!=0)
df <- data.frame(x = cell_types$refUMAP_1, y = cell_types$refUMAP_2,
                 Species = cell_types$clust7.predicted.id)

centroids <- aggregate(cbind(x,y)~Species,df,mean)
ggplot(df, aes(x, y, colour = Species)) +
  geom_point(size=0.25,alpha=0.1) +
  geom_text(data = centroids, aes(label = Species),color = "black") +
  xlab('UMAP1') + ylab('UMAP2') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



#2.4 Discrepancies between knn and default
df <- data.frame(x = cell_types$refUMAP_1, y = cell_types$refUMAP_2,
                 Species = cell_types$knn_clust7, classification = cell_types$clust7.predicted.id!=cell_types$knn_clust7 )
ggplot() + 
  geom_point(aes(x = x, y = y), data = df[df$classification==FALSE,], size=0.25,alpha=0.1, color= rgb(0.75,0.75,0.75) ) + # correctly classified
  geom_point(aes(x = x, y = y), data = df[df$classification==TRUE,], size=0.25,alpha=0.1, color= rgb(228,26,28,maxColorValue = 255) ) + # incorrectly classified
  xlab('UMAP1') + ylab('UMAP2') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

#sum(df$classification==FALSE)/nrow(df)

for (i in 1:length(uctps)){
  ct<- uctps[i]
  print(ct)
  df<-data.frame(cell_types[cell_types$knn_clust7==ct,c("refUMAP_1","refUMAP_2","knn_clust7","clust7.predicted.id")])
  
  ggplot() + 
    geom_point(aes(x = refUMAP_1, y = refUMAP_2), data = df[df$knn_clust7==df$clust7.predicted.id,], size=0.25,alpha=0.1, color= rgb(0.75,0.75,0.75) ) + # correctly classified
    geom_point(aes(x = refUMAP_1, y = refUMAP_2), data = df[df$knn_clust7!=df$clust7.predicted.id,], size=0.25,alpha=0.1, color= rgb(228,26,28,maxColorValue = 255) ) + # incorrectly classified
    xlab('UMAP1') + ylab('UMAP2') + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  b<-table(df$clust7.predicted.id)
  
  b<-100*b/sum(b)
 
   rownames(df[df$knn_clust7!=df$clust7.predicted.id,])

  
  rm(b,df)
}

