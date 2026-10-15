#0. clear workspace and Set working directory
rm(list = ls())
setwd('/Users/ibaidiez/Desktop/projects/ROSMAP/analysis/NatComm')

library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(miloR)
library(SingleCellExperiment)

seurat_file <- 'snRAseq_exp12/seurat_integrated/12.Rds'
rosmap_metadata <- 'rosmap_data.csv';
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
#2. Milo
################################
#It takes your UMAP/PCA graph, breaks it down into thousands of tiny overlapping "neighborhoods" of single cells, 
# and checks which specific neighborhoods are changing between AD and CONTROL.

DimPlot(seurat_merged, reduction = "umap", group.by = c("seurat_clusters") ,label = TRUE, pt.size = 0.05)


#2.1 Prepare seural data for milo
metadata$group <- relevel(metadata$group, ref = "CONTROL")
metadata$group <- factor( metadata$group, levels = c("CONTROL", "AD", "IMP_O", "RES", "Other"))

covars <- c( "group", "age_death", "msex", "pmi","hip_scl_3reg_yn", "lb_7reg", "tdp_st4", "caa_4gp")


keep_samples <- rownames(metadata)[
  metadata$group %in% c("CONTROL", "AD", "IMP_O","RES") & 
    complete.cases(metadata[, covars])
]


table(metadata[keep_samples, "group"])
length(keep_samples)


seurat_milo <- subset( seurat_merged, subset = individualID %in% keep_samples)
rm(seurat_merged)


# 2.2 Convert your Seurat object to a SingleCellExperiment object and create milo object
sce <- as.SingleCellExperiment(seurat_milo)
milo_obj <- Milo(sce)

# 2.3 Build the neighborhood graph  (Uses the integrated graph you already computed)
milo_obj <- buildGraph(milo_obj, k = 30, d = 30, reduced.dim = "INTEGRATED.CCA")
milo_obj <- makeNhoods(milo_obj, prop = 0.1, k = 30, d = 30, refined = TRUE, reduced_dims = "INTEGRATED.CCA")
milo_obj <- calcNhoodDistance(milo_obj,d = 30,reduced.dim = "INTEGRATED.CCA")

# 2.4 Count cells per neighborhood per individual
milo_obj <- countCells(milo_obj, meta.data = as.data.frame(colData(milo_obj)), samples = "individualID")


#2.5 Define  experimental design using  full formula
sample_ids <- colnames(nhoodCounts(milo_obj))
design_df <- metadata[match(sample_ids, rownames(metadata)), ]
rownames(design_df) <- sample_ids

design_df$group <- droplevels(design_df$group)
stopifnot(all(rownames(design_df) == colnames(nhoodCounts(milo_obj))))
stopifnot(!any(is.na(design_df[, covars])))
table(design_df$group)

#2.6 Compute statistics
library(statmod)

design_formula <- ~ 0 + group + age_death + msex + pmi + hip_scl_3reg_yn + lb_7reg + tdp_st4 + caa_4gp

milo_AD_vs_CONTROL <- testNhoods(
  milo_obj,
  design = design_formula,
  design.df = design_df,
  model.contrasts = "groupAD - groupCONTROL",
  robust = TRUE
)


milo_IMPO_vs_CONTROL <- testNhoods(
  milo_obj,
  design = design_formula,
  design.df = design_df,
  model.contrasts = "groupIMP_O - groupCONTROL",
  robust = TRUE
)

milo_RES_vs_AD <- testNhoods(
  milo_obj,
  design = design_formula,
  design.df = design_df,
  model.contrasts = "groupRES - groupAD",
  robust = TRUE
)


#2.7 Plot results
milo_obj <- buildNhoodGraph(milo_obj)

#plotNhoodGraphDA(milo_obj, milo_AD_vs_CONTROL,alpha = 0.1, layout = "UMAP")
#plotNhoodGraphDA(milo_obj, milo_IMPO_vs_CONTROL,alpha = 0.1, layout = "UMAP")

#annotate significant neighborhoods by Seurat cluster
milo_AD_vs_CONTROL <- annotateNhoods( milo_obj, milo_AD_vs_CONTROL, coldata_col = "seurat_clusters")
milo_AD_vs_CONTROL[order(milo_AD_vs_CONTROL$SpatialFDR), ][1:20, ]

#how many significant nodes in each cluster
table(
  milo_AD_vs_CONTROL$seurat_clusters[
    milo_AD_vs_CONTROL$SpatialFDR < 0.1
  ]
)


#2.7.1  Extract cell-level UMAP from the filtered Seurat object used for Milo
umap_df <- as.data.frame(Embeddings(seurat_milo, reduction = "umap"))
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$cell <- rownames(umap_df)
nhood_index <- nhoodIndex(milo_obj)
nhood_index <- unlist(nhood_index, use.names = FALSE)


#2.7.2 Get UMAP position of each neighborhood based on its index cell
nhood_umap <- umap_df[nhood_index, c("UMAP_1", "UMAP_2")]
nhood_umap$Nhood <- seq_len(nrow(nhood_umap))

#2.7.3 Plot AD va CONTROL
milo_res_df <- as.data.frame(milo_AD_vs_CONTROL)
nhood_umap_clean <- nhood_umap[, c("UMAP_1", "UMAP_2")] # Remove Nhood from coordinate table to avoid duplicate

milo_plot_df <- cbind(nhood_umap_clean, milo_res_df)
milo_plot_df$sig <- milo_plot_df$SpatialFDR < 0.1
milo_plot_df$neglog10SpatialFDR <- -log10(pmax(milo_plot_df$SpatialFDR, 1e-300))

ggplot() +
  geom_point(data = umap_df, aes(UMAP_1, UMAP_2),
    color = "grey90", size = 0.08, alpha = 0.35
  ) +
  geom_point(
    data = milo_plot_df[milo_plot_df$sig, ],
    aes(UMAP_1,UMAP_2,
      color = -logFC,size = neglog10SpatialFDR ),
    alpha = 0.9
  ) +
  scale_color_gradient2(midpoint = 0) + scale_size_continuous(range = c(1, 5)) +
  theme_classic() +
  labs(
    title = "Milo AD-associated microglia neighborhoods",
    subtitle = "Nodes only, overlaid on Seurat UMAP; SpatialFDR < 0.1",
    color = "logFC\nAD vs CONTROL",
    size = "-log10\nSpatialFDR"
  )

sig_ADvsControl_milo_data <- milo_plot_df[milo_plot_df$sig, ]


#2.7.4 Plot IMP-O vs CONTROL
milo_res_df <- as.data.frame(milo_IMPO_vs_CONTROL)
nhood_umap_clean <- nhood_umap[, c("UMAP_1", "UMAP_2")] # Remove Nhood from coordinate table to avoid duplicate

milo_plot_df <- cbind(nhood_umap_clean, milo_res_df)
milo_plot_df$sig <- milo_plot_df$SpatialFDR < 0.1
milo_plot_df$neglog10SpatialFDR <- -log10(pmax(milo_plot_df$SpatialFDR, 1e-300))

ggplot() +
  geom_point(
    data = umap_df,
    aes(UMAP_1, UMAP_2),
    color = "grey90",
    size = 0.08,
    alpha = 0.35
  ) +
  geom_point(
    data = milo_plot_df[milo_plot_df$sig, ],
    aes(
      UMAP_1,
      UMAP_2,
      color = -logFC,
      size = neglog10SpatialFDR
    ),
    alpha = 0.9
  ) +
  scale_color_gradient2(midpoint = 0) +
  scale_size_continuous(range = c(1, 5)) +
  theme_classic() +
  labs(
    title = "Milo IMP_O-associated microglia neighborhoods",
    subtitle = "Nodes only, overlaid on Seurat UMAP; SpatialFDR < 0.1",
    color = "logFC\nIMP-O vs CONTROL",
    size = "-log10\nSpatialFDR"
  )

sig_IMPOvsControl_milo_data <- milo_plot_df[milo_plot_df$sig, ]




#2.7.5 Plot RES vs AD
milo_res_df <- as.data.frame(milo_RES_vs_AD)
nhood_umap_clean <- nhood_umap[, c("UMAP_1", "UMAP_2")] # Remove Nhood from coordinate table to avoid duplicate

milo_plot_df <- cbind(nhood_umap_clean, milo_res_df)
#milo_plot_df$sig <- milo_plot_df$PValue < 0.05
#milo_plot_df$neglog10SpatialFDR <- -log10(pmax(milo_plot_df$PValue, 1e-300))
milo_plot_df$sig <- milo_plot_df$SpatialFDR < 0.2
milo_plot_df$neglog10SpatialFDR <- -log10(pmax(milo_plot_df$SpatialFDR, 1e-300))

ggplot() +
  geom_point(
    data = umap_df,
    aes(UMAP_1, UMAP_2),
    color = "grey90",
    size = 0.08,
    alpha = 0.35
  ) +
  geom_point(
    data = milo_plot_df[milo_plot_df$sig, ],
    aes(
      UMAP_1,
      UMAP_2,
      color = -logFC,
      size = neglog10SpatialFDR
    ),
    alpha = 0.9
  ) +
  scale_color_gradient2(midpoint = 0) +
  scale_size_continuous(range = c(1, 5)) +
  theme_classic() +
  labs(
    title = "Milo IMP_O-associated microglia neighborhoods",
    subtitle = "Nodes only, overlaid on Seurat UMAP; SpatialFDR < 0.1",
    color = "logFC\nRES vs AD",
    size = "-log10\nSpatialFDR"
  )

sig_RESvsAD_milo_data <- milo_plot_df[milo_plot_df$sig, ]



#2.8  Plot Resilience  with residual approach milo

#2.8.1 Compute resilience measure
fit_resilience <- lm(cogn_global ~ gpath_z, data = metadata, na.action = na.exclude)
metadata$resilience_gcog_gpath <- residuals(fit_resilience)
fit_resilience <- lm(cogn_global ~ tangles_z, data = metadata, na.action = na.exclude)
metadata$resilience_gcog_tangles <- residuals(fit_resilience)

fit_resilience <- lm(cts_mmse30_lv ~ gpath_z, data = metadata, na.action = na.exclude)
metadata$resilience_mmse_gpath <- residuals(fit_resilience)
fit_resilience <- lm(cts_mmse30_lv ~ tangles_z, data = metadata, na.action = na.exclude)
metadata$resilience_mmse_tangles <- residuals(fit_resilience)

fit_resilience <- lm(cogn_wo ~ gpath_z, data = metadata, na.action = na.exclude)
metadata$resilience_wo_gpath <- residuals(fit_resilience)
fit_resilience <- lm(cogn_wo ~ tangles_z, data = metadata, na.action = na.exclude)
metadata$resilience_wo_tangles <- residuals(fit_resilience)


covars <- c(
  "cogn_global", "cts_mmse30_lv", "cogn_wo", "gpath_z", "amyloid_z", "tangles_z",
  "resilience_gcog_gpath", "resilience_gcog_tangles", "resilience_mmse_gpath", "resilience_mmse_tangles", "resilience_wo_gpath", "resilience_wo_tangles",
  "age_death_z", "msex", "pmi_z", "hip_scl_3reg_yn", "lb_7reg", "tdp_st4", "caa_4gp"
)



keep_samples <- rownames(metadata)[complete.cases(metadata[, covars])]
metadata_subset <- metadata[keep_samples, ]

seurat_milo <- subset( seurat_milo, subset = individualID %in% keep_samples)


#2.8.2 Convert your Seurat object to a SingleCellExperiment object and create milo object
sce <- as.SingleCellExperiment(seurat_milo)
milo_obj <- Milo(sce)

#2.8.3 Build the neighborhood graph  (Uses the integrated graph you already computed)
milo_obj <- buildGraph(milo_obj, k = 30, d = 30, reduced.dim = "INTEGRATED.CCA")
milo_obj <- makeNhoods(milo_obj, prop = 0.1, k = 30, d = 30, refined = TRUE, reduced_dims = "INTEGRATED.CCA")
milo_obj <- calcNhoodDistance(milo_obj,d = 30,reduced.dim = "INTEGRATED.CCA")


#2.8.4 extract Milo neighborhood abundance per individual
milo_obj <- countCells( milo_obj, meta.data = as.data.frame(colData(milo_obj)), samples = "individualID" )
nhood_counts <- nhoodCounts(milo_obj) # extract the neighborhood count matrix:

common_ids <- intersect(colnames(nhood_counts), rownames(metadata_subset))
nhood_counts <- nhood_counts[, common_ids]
metadata_mod <- metadata_subset[common_ids, ]



#2.8.5 Compute statistics 
#library(statmod)

design_formula <- ~ 0 + resilience_gcog_gpath + age_death_z + msex + pmi_z +  hip_scl_3reg_yn + lb_7reg + tdp_st4 + caa_4gp
milo_resilience_gcog_gpath <- testNhoods( milo_obj, design = design_formula, design.df = metadata_mod, model.contrasts = "resilience_gcog_gpath", robust = TRUE )

design_formula <- ~ 0 + resilience_gcog_tangles + age_death_z + msex + pmi_z +  hip_scl_3reg_yn + lb_7reg + tdp_st4 + caa_4gp
milo_resilience_gcog_tangles <- testNhoods( milo_obj, design = design_formula, design.df = metadata_mod, model.contrasts = "resilience_gcog_tangles", robust = TRUE )

design_formula <- ~ 0 + resilience_mmse_gpath + age_death_z + msex + pmi_z +  hip_scl_3reg_yn + lb_7reg + tdp_st4 + caa_4gp
milo_resilience_mmse_gpath <- testNhoods( milo_obj, design = design_formula, design.df = metadata_mod, model.contrasts = "resilience_mmse_gpath", robust = TRUE )

design_formula <- ~ 0 + resilience_mmse_tangles + age_death_z + msex + pmi_z +  hip_scl_3reg_yn + lb_7reg + tdp_st4 + caa_4gp
milo_resilience_mmse_tangles <- testNhoods( milo_obj, design = design_formula, design.df = metadata_mod, model.contrasts = "resilience_mmse_tangles", robust = TRUE )

design_formula <- ~ 0 + resilience_wo_gpath + age_death_z + msex + pmi_z +  hip_scl_3reg_yn + lb_7reg + tdp_st4 + caa_4gp
milo_resilience_wo_gpath <- testNhoods( milo_obj, design = design_formula, design.df = metadata_mod, model.contrasts = "resilience_wo_gpath", robust = TRUE )

design_formula <- ~ 0 + resilience_wo_tangles + age_death_z + msex + pmi_z +  hip_scl_3reg_yn + lb_7reg + tdp_st4 + caa_4gp
milo_resilience_wo_tangles <- testNhoods( milo_obj, design = design_formula, design.df = metadata_mod, model.contrasts = "resilience_wo_tangles", robust = TRUE )



milo_obj <- buildNhoodGraph(milo_obj)


#2.8.6 plot
umap_df <- as.data.frame(Embeddings(seurat_milo, reduction = "umap"))
colnames(umap_df) <- c("UMAP_1", "UMAP_2")
umap_df$cell <- rownames(umap_df)
nhood_index <- nhoodIndex(milo_obj)
nhood_index <- unlist(nhood_index, use.names = FALSE)
nhood_umap <- umap_df[nhood_index, c("UMAP_1", "UMAP_2")]
nhood_umap$Nhood <- seq_len(nrow(nhood_umap))


# Plot resilience: milo_resilience_gcog_gpath, milo_resilience_gcog_tangles, milo_resilience_mmse_gpath, milo_resilience_mmse_tangles, milo_resilience_wo_gpath, milo_resilience_wo_tangles 

milo_res_df <- as.data.frame(milo_resilience_mmse_tangles)
nhood_umap_clean <- nhood_umap[, c("UMAP_1", "UMAP_2")] # Remove Nhood from coordinate table to avoid duplicate

milo_plot_df <- cbind(nhood_umap_clean, milo_res_df)
milo_plot_df$sig <- milo_plot_df$SpatialFDR < 0.1
milo_plot_df$neglog10SpatialFDR <- -log10(pmax(milo_plot_df$SpatialFDR, 1e-300))
#milo_plot_df$sig <- milo_plot_df$PValue < 0.05
#milo_plot_df$neglog10SpatialFDR <- -log10(pmax(milo_plot_df$PValue, 1e-300))

ggplot() +
  geom_point(
    data = umap_df,
    aes(UMAP_1, UMAP_2),
    color = "grey90",
    size = 0.08,
    alpha = 0.35
  ) +
  geom_point(
    data = milo_plot_df[milo_plot_df$sig, ],
    aes(
      UMAP_1,
      UMAP_2,
      color = -logFC,
      size = neglog10SpatialFDR
    ),
    alpha = 0.9
  ) +
  scale_color_gradient2(midpoint = 0) +
  scale_size_continuous(range = c(1, 5)) +
  theme_classic() +
  labs(
    title = "Milo AD-associated microglia neighborhoods",
    subtitle = "Nodes only, overlaid on Seurat UMAP; SpatialFDR < 0.1",
    color = "logFC\nGlobal pathology",
    size = "-log10\nSpatialFDR"
  )
sig_resh_milo_data <- milo_plot_df[milo_plot_df$sig, ]

