#0. clear workspace and Set working directory
rm(list = ls())
setwd('/Users/ibaidiez/Desktop/projects/ROSMAP/analysis/NatComm')

library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)

seurat_file <- 'snRAseq_exp12/seurat_integrated/12.Rds'
rosmap_metadata <- 'rosmap_data.csv';
groups_csv_file <- 'groups.csv'

num_dimensions_umap <- 50

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

#1.5 Load groups of available postmortem tissue
groups <- read.csv(groups_csv_file, header=T,stringsAsFactors = F)



################################
#2. Compute cluster and proportions
################################
#2.1  Find microglia subclusters for small resolution = 0.2
seurat_merged <- FindClusters(seurat_merged, resolution = 0.2 )
#seurat_merged <- FindClusters(seurat_merged, resolution = 1 )
#seurat_merged <- FindClusters(seurat_merged, resolution = 2.75 )

#seurat_merged <- RunUMAP(seurat_merged, dims = 1:num_dimensions_umap, reduction = "integrated.cca")
DimPlot(seurat_merged, reduction = "umap", group.by = c("seurat_clusters") ,label = TRUE, pt.size = 0.05)

#2.2 Export UMAP info to a table
#umap_coordinates <- Embeddings(seurat_merged, reduction = "umap")
#cluster_assignments <- seurat_merged$seurat_clusters
#export_table <- data.frame(
#  Cell_ID = rownames(umap_coordinates),
#  UMAP_1 = umap_coordinates[, "umap_1"],
#  UMAP_2 = umap_coordinates[, "umap_2"],
#  Cluster = cluster_assignments
#)
#write.csv(export_table, file = "umap_cluster_coordinates.csv", row.names = FALSE)
#rm(umap_coordinates,cluster_assignments,export_table)


#2.3 Plot markers
features_qc <- c( 'TOP2A', #proliferative
                  "P2RY12","CX3CR1",#Homeostatic
                  'IFI6','IFIT3',  #interferon
                  "SPP1","TMEM163", # Reactive
                  "HSPH1","DNAJB1", # stress response
                  'FTH1', 'FTL', #enhanced redox
                  'APOE', 'GPNMB',  #lipid associated
                  "MRC1", #Macrophages
                  "FCN1", #Monocytes
                  "CD247" #T cells
)
features_qc_sep <- c(1.5,3.5,5.5,7.5, 9.5,
                     11.5,13.5,14.5, 15.5)

#Reorder and plot
my_levels <- c(7, 6,0,  4,10, 1,5,2,3,9,11,8) 
Idents(seurat_merged) <- factor(Idents(seurat_merged), levels= my_levels)

DotPlot(seurat_merged, features = features_qc) +
  geom_vline(xintercept = features_qc_sep, linetype="dashed") +
  RotatedAxis()


#2.4 Get each state marker
for (i in 0:(nlevels(Idents(seurat_merged))-1)) {
  message(paste("Processing cluster:", i))
  markers_cluster <- FindMarkers(seurat_merged, ident.1 = i)
  write.csv(markers_cluster, file =  paste0("micr_markers_cluster_", i, ".csv"))
}


################################
#3. Compare with markers from other studies
################################

#3.1 Markers Paper Green et al: https://www.nature.com/articles/s41586-024-07871-6
features_qc <- c( 'TOP2A', #proliferative
                  'CX3CR1','P2RY12',#Surveiling
                  "TMEM163","IL4R", # Reactive
                  "HSPH1","DNAJB1","NLRP1", # stress response
                  'FTH1', 'FTL', #enhanced redox
                  'APOE', 'ADAM10', 'GPNMB','CPM', 'SPP1', 'TREM2' ,  #lipid associated
                  'IFI6','IFIT3','IFIT2',  #interferon
                  'CD83','CCL3','CCL4', #inflammatory
                  "SERPINE1" #serpinE
)

features_qc_sep <- c(1.5,3.5,5.5,8.5, 10.5,16.5,19.5,22.5,25.5)

DotPlot(seurat_merged, features = features_qc) +
  geom_vline(xintercept = features_qc_sep, linetype="dashed") +
  RotatedAxis()


#3.2 Markers Paper Gabitto et al: https://www.nature.com/articles/s41593-024-01774-5
features_qc <- c( 'MRC1', 'RGL1', 'AGFG1', #Micro-PVM_1
                  'FCN1', 'CD300E', 'LILRA5',#Monocyte
                  'FAM177B', 'DLEU7', 'ITM2B', #Micro-PVM_2
                  'SPC25', 'CIT', 'NSD2',#Micro-PVM_2_1
                  'SKAP1', 'SLC38A1', 'PARP8', # Lymphocyte
                  'CSMD1', 'KCNIP4', 'SH3RF3', #Micro-PVM_2_3
                  'PTPRG', 'TRHDE', 'CADM1',  #Micro-PVM_3
                  'RAPGEF1', 'PDE8A', 'MAP4K3' #Micro-PVM_4 
)

features_qc_sep <- c(3.5,6.5,9.5,12.5, 15.5,18.5,21.5)

DotPlot(seurat_merged, features = features_qc) +
  geom_vline(xintercept = features_qc_sep, linetype="dashed") +
  RotatedAxis()


#3.3 Markers Paper Cain et al: 
features_qc <- c(
  'FRMD4A', 'P2RY12', 'SYNDIG1', 'CX3CR1', # homeostatic
  'ARHGAP11B', 'HMMR', 'TOP2A',  # proliferative
  'APOE', 'FTL', 'RPL19', 'RPS6', 'FTH1', 'RPS19', 'C1QC', 'C1QB', 'C1QA', #enchanced redox
  'IFIT3', 'IFI44L', 'IFI6', 'PTPRG', #Interferon response
  'DDIT4', 'SLC2A3', 'TMEM163', 'SPP1' #Anti inflamatory
)
features_qc_sep <- c(4.5,7.5,16.5,20.5)

DotPlot(seurat_merged, features = features_qc) +
  geom_vline(xintercept = features_qc_sep, linetype="dashed") +
  RotatedAxis()



#3.4 Markers Paper Gerritis et al
features_qc <- c( 
  'APOE', 'GPNMB', 'CPM', 'STARD13', 'EYA2', #AB associated c10
  'ADAMTS17', 'IPCEF1', 'SOCS6', 'DSCAM', #AB associated c7
  'SPATS2L', 'DTNA', 'TPRG1', 'PTPRG', 'MYO1E', 'GLDN', 'ITGAX', 'LPL',  #AB associated c9
  'DDX17', 'FOXP2', 'GRID2', 'NAV2', 'DPP10', 'ADGRB3', #Tangles associated
  'TMEM163', 'SH3RF3', 'P2RY12', 'CX3CR1', 'SOX5', #homeostatic
  'NFKB1', 'CCL3', 'CCL4', 'SLC2A3', 'SRGN', 'IL1B', 'GPR183', 'CD83', #inflammation
  'FANCI', 'BRIP1', 'TOP2A', 'MKI67', #proliferative
  'JUNB', 'FOS',   'HIST1H2AC', 'HSP90AA1', 'HSPA1A', 'HSPA1B' #stress response
  
)
features_qc_sep <- c(5.5,9.5,17.5,23.5, 28.5,36.5,40.5)
DotPlot(seurat_merged, features = features_qc) +
  geom_vline(xintercept = features_qc_sep, linetype="dashed") +
  RotatedAxis()



#3.5 Markers Paper Martins-Ferreira et al # https://www.nature.com/articles/s41467-025-56124-1
features_qc <- c( 
  'P2RY12', 'KBTBD12', 'CX3CR1', #Homeostatic 1
  'GRID2', 'PRDM11',# Homeostatic 2
  'SERPINE1', 'SELPLG', #Homeostatic 3
  'SLC2A3', 'CD83', #DIMs
  'TMEM163', 'HAMP', #Inflam.DAM
  'SYT1', 'PCDH9',#Ribo.DAM1
  'PLEKHA7', 'MECOM', 'OOEP', #Ribo.DAM2
  'GPNMB', 'PTPRG', 'CD9', #Lipo.DAM
  'F13A1', 'CD163', 'IQGAP2' #MAC
)
features_qc_sep <- c(3.5,5.5,7.5,9.5,11.5, 13.5,16.5,19.5)
DotPlot(seurat_merged, features = features_qc) +
  geom_vline(xintercept = features_qc_sep, linetype="dashed") +
  RotatedAxis()



################################################################
#4. Normalized Substate proportions
################################################################

#4.1 compute cluster proportions
df<-t(table(seurat_merged@meta.data$seurat_clusters, seurat_merged@meta.data$individualID))
df_normalized <-  df / rowSums(df)

#4.2  Rank-based inverse normal transformation within each column
for (p in seq_len(ncol(df_normalized))) {
  y_ranked <- rank(df_normalized[, p], ties.method = "average")
  df_normalized[, p] <- qnorm((y_ranked - 0.5) / length(y_ranked))
}

rm(df)

################################################################
#5.  Subclusters association with pathology
################################################################

#Prepare design matrix and contrast for statistics
df <- metadata[rownames(df_normalized),]
df$caa_4gp <- as.integer(df$caa_4gp > 1)
covariates <- df[,c("age_death", "msex", "pmi", "hip_scl_3reg_yn", "lb_7reg", "tdp_st4", "ci_num2_gct", "caa_4gp" )]
XX <- cbind(matrix(1, nrow = nrow(covariates), ncol = 2), covariates)
contrast <- c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0)

#5.1 Amyloid
pathology <- df$amyloid
#Normalize pathology data
y_ranked <- rank(pathology, ties.method = "average")
pathology <- qnorm((y_ranked - 0.5) / length(y_ranked))

#Use complete data
indx <- which(complete.cases(cbind(pathology, XX)))
X <- as.matrix(XX[indx, , drop = FALSE])
X[, 2] <- pathology[indx] - mean(pathology[indx])
X[is.na(X[, 2]), 2] <- 0
df_residual <- nrow(X) - ncol(X)

#Compute statistics
y <- as.matrix(df_normalized[indx, , drop = FALSE])
betas <- MASS::ginv(X) %*% y
residuals <- y - X %*% betas
MSE <- colSums(residuals^2) / df_residual
contrast_variance <- as.numeric( t(contrast) %*% solve(crossprod(X)) %*% contrast )
SE <- sqrt(MSE * contrast_variance)
T <- as.numeric(t(contrast) %*% betas) / SE
pval <- 2 * pt(abs(T), df = df_residual, lower.tail = FALSE)
Z <- sign(T) * qnorm(pval / 2, lower.tail = FALSE)
for (c in seq_along(pval)) {
  cat(sprintf("Amyloid - Cluster %d T=%.2f Z=%.2f p=%.4g\n", c - 1, T[c], Z[c], pval[c]))
}


#5.2 Tangles
pathology <- df$tangles
#Normalize pathology data
y_ranked <- rank(pathology, ties.method = "average")
pathology <- qnorm((y_ranked - 0.5) / length(y_ranked))

#Use complete data
indx <- which(complete.cases(cbind(pathology, XX)))
X <- as.matrix(XX[indx, , drop = FALSE])
X[, 2] <- pathology[indx] - mean(pathology[indx])
X[is.na(X[, 2]), 2] <- 0
df_residual <- nrow(X) - ncol(X)

#Compute statistics
y <- as.matrix(df_normalized[indx, , drop = FALSE])
betas <- MASS::ginv(X) %*% y
residuals <- y - X %*% betas
MSE <- colSums(residuals^2) / df_residual
contrast_variance <- as.numeric( t(contrast) %*% solve(crossprod(X)) %*% contrast )
SE <- sqrt(MSE * contrast_variance)
T <- as.numeric(t(contrast) %*% betas) / SE
pval <- 2 * pt(abs(T), df = df_residual, lower.tail = FALSE)
Z <- sign(T) * qnorm(pval / 2, lower.tail = FALSE)
for (c in seq_along(pval)) {
  cat(sprintf("Tangles - Cluster %d T=%.2f Z=%.2f p=%.4g\n", c - 1, T[c], Z[c], pval[c]))
}



#5.3 Global pathology
pathology <- df$gpath
#Normalize pathology data
y_ranked <- rank(pathology, ties.method = "average")
pathology <- qnorm((y_ranked - 0.5) / length(y_ranked))

#Use complete data
indx <- which(complete.cases(cbind(pathology, XX)))
X <- as.matrix(XX[indx, , drop = FALSE])
X[, 2] <- pathology[indx] - mean(pathology[indx])
X[is.na(X[, 2]), 2] <- 0
df_residual <- nrow(X) - ncol(X)

#Compute statistics
y <- as.matrix(df_normalized[indx, , drop = FALSE])
betas <- MASS::ginv(X) %*% y
residuals <- y - X %*% betas
MSE <- colSums(residuals^2) / df_residual
contrast_variance <- as.numeric( t(contrast) %*% solve(crossprod(X)) %*% contrast )
SE <- sqrt(MSE * contrast_variance)
T <- as.numeric(t(contrast) %*% betas) / SE
pval <- 2 * pt(abs(T), df = df_residual, lower.tail = FALSE)
Z <- sign(T) * qnorm(pval / 2, lower.tail = FALSE)
for (c in seq_along(pval)) {
  cat(sprintf("Global pathology - Cluster %d T=%.2f Z=%.2f p=%.4g\n", c - 1, T[c], Z[c], pval[c]))
}



################################################################
#6.  Group differences
################################################################

#Prepare design matrix and contrast for statistics
df <- metadata[rownames(df_normalized),]
df$caa_4gp <- as.integer(df$caa_4gp > 1)
covariates <- df[,c("age_death", "msex", "pmi", "hip_scl_3reg_yn", "lb_7reg", "tdp_st4", "ci_num2_gct", "caa_4gp" )]
XX <- cbind(matrix(0, nrow = nrow(covariates), ncol = 3), covariates)
XX[df$group=="AD",1] <-1
XX[df$group=="IMP_O",2] <-1
XX[df$group=="CONTROL",3] <-1
indx <- which(complete.cases(XX))
X <- as.matrix(XX[indx, , drop = FALSE])
df_residual <- nrow(X) - ncol(X)

#6.1 AD vs control
contrast <- c(1,0, -1, 0, 0, 0, 0, 0, 0, 0, 0)
y <- as.matrix(df_normalized[indx, , drop = FALSE])
betas <- MASS::ginv(X) %*% y
residuals <- y - X %*% betas
MSE <- colSums(residuals^2) / df_residual
contrast_variance <- as.numeric( t(contrast) %*% solve(crossprod(X)) %*% contrast )
SE <- sqrt(MSE * contrast_variance)
T <- as.numeric(t(contrast) %*% betas) / SE
pval <- 2 * pt(abs(T), df = df_residual, lower.tail = FALSE)
Z <- sign(T) * qnorm(pval / 2, lower.tail = FALSE)
for (c in seq_along(pval)) {
  cat(sprintf("AD vs Control - Cluster %d T=%.2f Z=%.2f p=%.4g\n", c - 1, T[c], Z[c], pval[c]))
}


#6.2 IMP-O vs control
contrast <- c( 0, 1, -1, 0, 0, 0, 0, 0, 0, 0, 0)
y <- as.matrix(df_normalized[indx, , drop = FALSE])
betas <- MASS::ginv(X) %*% y
residuals <- y - X %*% betas
MSE <- colSums(residuals^2) / df_residual
contrast_variance <- as.numeric( t(contrast) %*% solve(crossprod(X)) %*% contrast )
SE <- sqrt(MSE * contrast_variance)
T <- as.numeric(t(contrast) %*% betas) / SE
pval <- 2 * pt(abs(T), df = df_residual, lower.tail = FALSE)
Z <- sign(T) * qnorm(pval / 2, lower.tail = FALSE)
for (c in seq_along(pval)) {
  cat(sprintf("IMP-O vs Control - Cluster %d T=%.2f Z=%.2f p=%.4g\n", c - 1, T[c], Z[c], pval[c]))
}




################################################################
#7.  Moderation resilience
################################################################

#7.1 Generate design matrix
df<-t(table(seurat_merged@meta.data$seurat_clusters, seurat_merged@meta.data$individualID))
df_normalized <-  df / rowSums(df)

df <- metadata[rownames(df_normalized),]
pathology <- df$gpath
#y_ranked <- rank(pathology, ties.method = "average")
#pathology <- qnorm((y_ranked - 0.5) / length(y_ranked))

y <- df$cts_mmse30_lv
XX <- cbind(matrix(1, nrow = nrow(covariates), ncol = 4), covariates)
indx <- which(complete.cases(cbind(pathology,y, XX)))
X <- as.matrix(XX[indx, , drop = FALSE]) #Use complete data
X[, 2] <- pathology[indx] - mean(pathology[indx])
X[is.na(X[, 2]), 2] <- 0
df_residual <- nrow(X) - ncol(X)
y <- y[indx]
contrast <- c( 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0 )


#7.2 For each cluster compute statistics
for (c in seq_len(ncol(df_normalized)) - 1L) {
  #Generate X and df
  X[, 3] <- df_normalized[indx,c+1] -  mean(df_normalized[indx, c + 1], na.rm = TRUE)
  X[is.na(X[, 3]), 3] <- 0
  X[, 4] <- X[, 2]  * X[, 3] 
  
  betas <- MASS::ginv(X) %*% y
  residuals <- y - X %*% betas
  MSE <- colSums(residuals^2) / df_residual
  contrast_variance <- as.numeric( t(contrast) %*% solve(crossprod(X)) %*% contrast )
  SE <- sqrt(MSE * contrast_variance)
  T <- as.numeric(t(contrast) %*% betas) / SE
  #pval <- 2 * pt(abs(T), df = df_residual, lower.tail = FALSE)
  #Z <- sign(T) * qnorm(pval / 2, lower.tail = FALSE)
  pval <- pt(T, df = df_residual, lower.tail = FALSE)
  Z <- qnorm(pval, lower.tail = FALSE)
  
  cat(sprintf("Moderation - Cluster %d T=%.2f Z=%.2f p=%.4g\n", c, T, Z, pval))
  
}




################################
#8. CD68
################################

#8.1 Load CD68 module correlations
cd68_cluster <- read.csv('/Users/ibaidiez/Desktop/papers/kumar/Nature_com/Review/bulk/CD68_Module_Correlations.csv', header=TRUE,stringsAsFactors = F)

#8.2 Create an aggregate CD68 clsuter expression
seurat_merged <- AddModuleScore(object = seurat_merged, features = list(cd68_cluster$Gene), name = "CD68_cluster_score" )


#8.3 Project in a UMAP the CD68 projected expression
# Extract UMAP coordinates
umap_df <- as.data.frame(Embeddings(seurat_merged, reduction = "umap"))
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$cell <- rownames(umap_df)
umap_df$CD68_cluster_score1 <- seurat_merged$CD68_cluster_score1
q50 <- quantile(umap_df$CD68_cluster_score1, 0.50, na.rm = TRUE)
q95 <- quantile(umap_df$CD68_cluster_score1, 0.95, na.rm = TRUE)

umap_df$CD68_plot <- pmin(pmax(umap_df$CD68_cluster_score1, q50), q95 )
umap_df$CD68_size <- scales::rescale(umap_df$CD68_plot, to = c(0.1, 1.8) )

# Plot
ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = CD68_plot)) +
  geom_point(size = 0.08, alpha = 0.8) +
  scale_color_gradient(low = "grey90", high = "blue") +
  theme_classic() +
  labs(
    title = "CD68-associated microglia score",
    subtitle = "Seurat UMAP; color clipped at q50–q95",
    color = "CD68 score")



#8.4 CD68 score per subject

#8.4.1. Pull the module score and the subject ID into a data frame
subject_scores_df <- FetchData(
  object = seurat_merged, 
  vars = c("CD68_cluster_score1", "individualID") 
)

#8.4.2. Group by subject and calculate the mean (or median) score
subject_summary <- subject_scores_df %>%
  group_by(individualID) %>%
  summarize(
    Mean_CD68_Score = mean(CD68_cluster_score1),
    Median_CD68_Score = median(CD68_cluster_score1),
    Cell_Count = n() # Handy to track how many cells each subject contributed
  )

rownames(subject_summary) <- subject_summary$individualID

#8.5 Statistics of CD68 vs cluster proportions
XX <- cbind(matrix(1, nrow = nrow(covariates), ncol = 2), covariates)
contrast <- c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0)

pathology <- subject_summary[rownames(covariates),]$Mean_CD68_Score

#Use complete data
indx <- which(complete.cases(cbind(pathology, XX)))
X <- as.matrix(XX[indx, , drop = FALSE])
X[, 2] <- pathology[indx] - mean(pathology[indx])
X[is.na(X[, 2]), 2] <- 0
df_residual <- nrow(X) - ncol(X)

#Compute statistics
y <- as.matrix(df_normalized[indx, , drop = FALSE])
betas <- MASS::ginv(X) %*% y
residuals <- y - X %*% betas
MSE <- colSums(residuals^2) / df_residual
contrast_variance <- as.numeric( t(contrast) %*% solve(crossprod(X)) %*% contrast )
SE <- sqrt(MSE * contrast_variance)
T <- as.numeric(t(contrast) %*% betas) / SE
pval <- 2 * pt(abs(T), df = df_residual, lower.tail = FALSE)
Z <- sign(T) * qnorm(pval / 2, lower.tail = FALSE)
for (c in seq_along(pval)) {
  cat(sprintf("CD68 - Cluster %d T=%.2f Z=%.2f p=%.4g\n", c - 1, T[c], Z[c], pval[c]))
}



################################
#9. High resolution populations
################################

#9.1 Higher resolution cluster
seurat_merged <- FindClusters(seurat_merged, resolution = 2.75 )
DimPlot(seurat_merged, reduction = "umap", group.by = c("seurat_clusters") ,label = TRUE, pt.size = 0.05)


#9.2 AD vs Resilient
df<-t(table(seurat_merged@meta.data$seurat_clusters, seurat_merged@meta.data$individualID))
df_normalized <- 100 * df / rowSums(df)
for ( i in 1:ncol(df) )
{
  AD_data <- df_normalized[ rownames(df) %in% groups[groups$Group=="AD",c("id")], i]
  RESIL_data <- df_normalized[ rownames(df) %in% groups[groups$Group=="RESIL",c("id")], i]
  
  wtest_result <- wilcox.test(AD_data, RESIL_data, alternative = "two.sided", conf.level = 0.95)
  print(paste0(i-1,". p=", wtest_result$p.value))
}


#9.3 Plot significant cluster
sig_clust <- 26
df <- data.frame(x = seurat_merged@reductions$umap@cell.embeddings[,c("umap_1")], y = seurat_merged@reductions$umap@cell.embeddings[,c("umap_2")], individualID = seurat_merged$individualID, cluster = seurat_merged$seurat_clusters )

ggplot() +
  geom_point(aes(x = x, y = y), data = df, size=0.05,alpha=0.5, color= rgb(0.75,0.75,0.75) ) + # correctly classified
  geom_point(aes(x = x, y = y), data = df[df$cluster==sig_clust,], size=0.05,alpha=0.5, color= rgb(228,26,28,maxColorValue = 255) ) + # Cluster 5
  xlab('UMAP1') + ylab('UMAP2') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



#9.4 Boxplot
AD_data <- df_normalized[ rownames(df) %in% groups[groups$Group=="AD",c("id")], i]
RESIL_data <- df_normalized[ rownames(df) %in% groups[groups$Group=="RESIL",c("id")], i]
FRAIL_data <- df_normalized[ rownames(df) %in% groups[groups$Group=="FRAIL",c("id")], i]
CONTRL_data <- df_normalized[ rownames(df) %in% groups[groups$Group=="CONTRL",c("id")], i]

indx_outliers<-c("R1521786","R2787688","R1627792","R4937504", "R6913259") #"R7516349"

y <-df_normalized[ rownames(df) %in% groups$id,i]
rownames(groups) <- groups$id
df2 <- data.frame(x = groups[names(y),c("Group")], y = y )


df2<-df2[ !(rownames(df2) %in% indx_outliers),  ] # remove outliers
AD_data <- df2[ df2$x=="AD", c("y")]
RESIL_data <-  df2[ df2$x=="RESIL", c("y")]
FRAIL_data <- df2[ df2$x=="FRAIL", c("y")]
CONTRL_data <- df2[ df2$x=="CONTRL", c("y")]


ggplot(df2, aes(x = x, y = y)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, height = 0) +
  labs(title = "Boxplot with All Data Points", x = "Group", y = "Value")




#9.5 Find markers
markers_cluster <- FindMarkers(seurat_merged, ident.1 = 26)
write.csv(markers_cluster, file="resil1_markers.csv")


#9.6 Plot distribution of cells in each subgroup
group <- "AD" 
group_color <- rgb(228,26,28,maxColorValue = 255) # red
ggplot() +
  geom_point(aes(x = x, y = y), data = df, size=0.05,alpha=0.5, color= rgb(0.75,0.75,0.75) ) + # correctly classified
  geom_point(aes(x = x, y = y), data = df[df$individualID %in% groups[groups$Group==group,c("id")], ], size=0.5,alpha=1, color= group_color ) + # Cluster 5
  xlab('UMAP1') + ylab('UMAP2') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


group <- "RESIL" 
group_color <-rgb(77,175,74,maxColorValue = 255) # green
ggplot() +
  geom_point(aes(x = x, y = y), data = df, size=0.05,alpha=0.5, color= rgb(0.75,0.75,0.75) ) + # correctly classified
  geom_point(aes(x = x, y = y), data = df[df$individualID %in% groups[groups$Group==group,c("id")], ], size=0.5,alpha=1, color= group_color ) + # Cluster 5
  xlab('UMAP1') + ylab('UMAP2') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



group <- "FRAIL"
group_color <-rgb(152,78,163,maxColorValue = 255) # purple
ggplot() +
  geom_point(aes(x = x, y = y), data = df, size=0.05,alpha=0.5, color= rgb(0.75,0.75,0.75) ) + # correctly classified
  geom_point(aes(x = x, y = y), data = df[df$individualID %in% groups[groups$Group==group,c("id")], ], size=0.5,alpha=1, color= group_color ) + # Cluster 5
  xlab('UMAP1') + ylab('UMAP2') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


group <- "CONTRL"
group_color <-rgb(33,102,172,maxColorValue = 255) # blue
ggplot() +
  geom_point(aes(x = x, y = y), data = df, size=0.05,alpha=0.5, color= rgb(0.75,0.75,0.75) ) + # correctly classified
  geom_point(aes(x = x, y = y), data = df[df$individualID %in% groups[groups$Group==group,c("id")], ], size=0.5,alpha=1, color= group_color ) + # Cluster 5
  xlab('UMAP1') + ylab('UMAP2') + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))








