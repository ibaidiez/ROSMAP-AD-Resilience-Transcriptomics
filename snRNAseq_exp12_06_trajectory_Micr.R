#0. clear workspace and Set working directory
rm(list = ls())
setwd('/Users/ibaidiez/Desktop/projects/ROSMAP/analysis/NatComm')

library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(ggplot2)


seurat_file <- 'snRAseq_exp12/seurat_integrated/12.Rds'
rosmap_metadata <- 'rosmap_data.csv';

################################
#1. Load data
################################

#1.1 Load preprocessed snRNAseq data
seurat_merged <- LoadSeuratRds(seurat_file)

#1.2 Load metadata
metadata <- read.csv(rosmap_metadata,sep = "\t", header=T,stringsAsFactors = F)
rownames(metadata) <- metadata$individualID


#1.3 Select metadata of participant with snRNASEQ
common_ids <- intersect(rownames(metadata), unique(seurat_merged@meta.data$individualID))
metadata <- metadata[common_ids, ]
rm(common_ids)

#1.4 Prepare variables and groups
metadata$sex <- as.factor(metadata$msex)
metadata$age_death[metadata$age_death == "90+"] <- "90" # Replace "90+" with "90"
metadata$age_death <- as.numeric(metadata$age_death)    # Convert the entire column to numbers
metadata$hip_scl_3reg_yn <- factor(metadata$hip_scl_3reg_yn)
metadata$lb_7reg <- factor(metadata$lb_7reg)
metadata$tdp_st4 <- factor(metadata$tdp_st4)
metadata$caa_4gp <- factor(metadata$caa_4gp)
metadata$gpath <- as.numeric(metadata$gpath)
metadata$tangles <- as.numeric(metadata$tangles)
metadata$amyloid <- as.numeric(metadata$amyloid)
metadata$cogn_global <- as.numeric(metadata$cogn_global)
metadata$cts_mmse30_lv <- as.numeric(metadata$cts_mmse30_lv)
metadata$cogn_wo <- as.numeric(metadata$cogn_wo)
metadata$gpath_z <- as.numeric(scale(metadata$gpath))
metadata$amyloid_z <- as.numeric(scale(metadata$amyloid))
metadata$tangles_z <- as.numeric(scale(metadata$tangles))
metadata$age_death_z <- as.numeric(scale(metadata$age_death))
metadata$pmi_z <- as.numeric(scale(metadata$pmi))

metadata$group <- "Other"

#1.4.1 AD
ad_condition <- !is.na(metadata$braaksc) & metadata$braaksc >= 5 &   !is.na(metadata$cogdx) & metadata$cogdx %in% c(4, 5)
metadata$group[ad_condition] <- "AD"

#1.4.2 Control
control_condition <- !is.na(metadata$braaksc) & metadata$braaksc <= 3 &   !is.na(metadata$cogdx) & metadata$cogdx == 1 &   !is.na(metadata$ad_reagan) & metadata$ad_reagan == 0
metadata$group[control_condition] <- "CONTROL"

#1.4.3 IMP_O
frail_condition <- !is.na(metadata$braaksc) & metadata$braaksc <= 3 &   !is.na(metadata$cogdx) & metadata$cogdx %in% c(4, 5) &   !is.na(metadata$ad_reagan) & metadata$ad_reagan == 0
metadata$group[frail_condition] <- "IMP_O"

#1.4.4 RES
resil_condition <- !is.na(metadata$braaksc) & metadata$braaksc >= 5 &    !is.na(metadata$cogdx) & metadata$cogdx == 1  &   !is.na(metadata$ad_reagan) & metadata$ad_reagan == 1
metadata$group[resil_condition] <- "RES"

metadata$group <- as.factor(metadata$group)
metadata$group <- relevel(metadata$group, ref = "CONTROL")
table(metadata$group)

rm(ad_condition,control_condition,frail_condition)


################################
#2. Trajectory analysis
################################

#2.1 Clusters
seurat_merged <- FindClusters(seurat_merged, resolution = 0.2)
DimPlot(seurat_merged, reduction = "umap", group.by = c("seurat_clusters") ,label = TRUE, pt.size = 0.05)

#2.2 Convert data to monocle 3
cds <- as.cell_data_set(seurat_merged)
rowData(cds)$gene_short_name <- row.names(rowData(cds))


#2.3 Transfer projected feature loadings
seurat_merged <- ProjectDim(object = seurat_merged, reduction = "pca")
cds@reduce_dim_aux[['PCA']][['model']][['svd_v']] <- seurat_merged@reductions[["pca"]]@feature.loadings.projected
cds@reduce_dim_aux[['PCA']][['model']][['svd_sdev']] <- seurat_merged@reductions$pca@stdev
reducedDim(cds, type = "PCA") <-   seurat_merged@reductions$pca@cell.embeddings
rm(seurat_merged)


#2.4 Find cluster and partitions
cds <-cluster_cells(cds, reduction_method = "UMAP", resolution=2e-5) # Create partitions (UMAP based)
#plot_cells(cds)

#2.5 Compute trajectory
cds <- learn_graph(cds, use_partition = FALSE) # Learn graph (UMAP based)

#2.6 Plot trajectory
plot_cells(cds,  color_cells_by = "cluster",
           label_groups_by_cluster=FALSE, label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           label_cell_groups = FALSE,)


#2.7 Choose root node
cds <- order_cells(cds)


#2.8 Analyse branches
cds_sub <- choose_graph_segments(cds)
cds_sub <- cds[,colnames(cds_sub)]



#2.9 Dynamic expression over path
x<-pseudotime(cds_sub)
cds_sub <- estimate_size_factors(cds_sub)
cds_exprs <- counts(cds_sub)/size_factors(cds_sub)


#cds_sub_micr5 <- cds_sub
#cds_sub_micr_ant <- cds_sub

gdata <- data.frame(x = x, 
                    y1 = cds_exprs['CD83',], 
                    y2 = cds_exprs['C1QC',],
                    y3 = cds_exprs['DUSP1',]) 


ggplot(gdata, aes(x)) +   
  geom_smooth(aes(y = y1), color = "black") + 
  geom_smooth(aes(y = y2), color = "red", fill='red') +
  geom_smooth(aes(y = y3), color = "blue", fill='blue') +
  labs(x = "pseudotime", y="expression", title="Dynamic Expression") +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))



#2.10 Check moran statistics of the genes in the trajectory
cds_sub_graph_test_results <- graph_test(cds_sub,neighbor_graph = "principal_graph", cores = 8)
trajectory_genes <- cds_sub_graph_test_results[ order(cds_sub_graph_test_results$morans_I, decreasing = TRUE), ]

sig_trajectory_genes <- subset(trajectory_genes, q_value < 0.05)
trajectory_genes[c("CD83", "DUSP1", "C1QC"), ]



