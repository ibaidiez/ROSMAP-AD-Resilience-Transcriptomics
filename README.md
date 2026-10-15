# Tau Pathology in Synapses, Microglial Gene Expression Trajectories, and Resilience to Alzheimer Disease

## Overview

This repository contains R scripts supporting the bulk RNA-sequencing (bulk RNA-seq) and single-nucleus RNA-sequencing (snRNA-seq) analyses performed for the study **“Tau Pathology in Synapses, Microglial Gene Expression Trajectories, and Resilience to Alzheimer Disease”** (Kumar et al.).

The analyses use data from the **Religious Orders Study and Rush Memory and Aging Project (ROSMAP)** to investigate transcriptomic signatures associated with Alzheimer disease neuropathology, microglial states, and cognitive resilience. The single-nucleus workflow integrates two ROSMAP snRNA-seq datasets and focuses on immune-cell and microglial heterogeneity. The bulk workflow includes differential gene expression and the identification of a CD68-associated co-expression module.

**Important:** This repository contains analysis code only. **Individual-level ROSMAP data, expression matrices, clinical metadata, and derived participant-level files are not provided.** Access to the underlying ROSMAP data requires the applicable authorization and **data-use agreement (DUA) with Rush University / the Rush Alzheimer's Disease Center**. See [Data availability and access](#data-availability-and-access).

## Repository files

### Single-nucleus RNA-seq workflow

The scripts are numbered in their intended processing order.

| Script | Description |
| --- | --- |
| `snRNAseq_exp12_01_mapping.R` | Processes the Experiment 1 reference, annotates major brain cell types, and projects/maps cells from Experiment 2 onto the reference using Seurat. |
| `snRNAseq_exp12_02_cell_subtyping_qc_Micr.R` | Selects the microglial/immune-cell compartment and performs quality-control checks, including outlier and doublet assessment. |
| `snRNAseq_exp12_03_cell_integration_Micr.R` | Integrates microglial/immune-cell data across experimental libraries using Seurat, performs dimensionality reduction and clustering, and exports integrated objects and cell counts. |
| `snRNAseq_exp12_04_substates_Micr.R` | Identifies and characterizes immune/microglial states at different clustering resolutions; assesses their associations with pathology, diagnosis, cognition, and resilience; and projects the CD68 co-expression module. |
| `snRNAseq_exp12_05_milo_Micr.R` | Performs Milo neighborhood-level differential abundance analyses for diagnostic groups, pathology, and residual-based cognitive resilience. |
| `snRNAseq_exp12_06_trajectory_Micr.R` | Uses Monocle 3 to infer microglial transcriptional trajectories and examine pseudotime-associated gene expression. |

### Bulk RNA-seq workflow

| Script | Description |
| --- | --- |
| `RNAseq_bulk_differential_expression.R` | Performs bulk RNA-seq differential gene expression analysis, including the comparison between resilient and demented AD cases, using DESeq2. |
| `RNAseq_bulk_coexpression_network_cd68.R` | Processes bulk RNA-seq data, estimates a gene co-expression network, identifies the CD68-containing module, and performs functional enrichment analyses. |

## Software requirements

The scripts were written in **R**. Major packages used include:

- **Single-nucleus analyses:** `Seurat`, `SeuratWrappers`, `Matrix`, `data.table`, `tidyverse`, `DoubletFinder`, `dbscan`, `ggdist`, `patchwork`, `miloR`, `SingleCellExperiment`, `monocle3`, `statmod`, `MASS`, and `scales`.
- **Bulk analyses:** `DESeq2`, `limma`, `dynamicTreeCut`, `clusterProfiler`, `org.Hs.eg.db`, `dplyr`, `tibble`, and `ggplot2`.

Some packages are distributed through Bioconductor or their respective development repositories rather than CRAN. Install compatible package versions for your R environment before execution. **A lockfile and a validated, one-command installation environment are not supplied.** Recording the installed package versions with `sessionInfo()` is recommended when reproducing results.



## Citation

If you use these scripts, please cite the associated study:

**Kumar et al. _Tau Pathology in Synapses, Microglial Gene Expression Trajectories, and Resilience to Alzheimer Disease._** [Publication details/DOI to be added when available.]

## Contact

For questions about the analysis code, please contact the corresponding authors listed in the associated manuscript.
