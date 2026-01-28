### based on TRIGR_mergedCells_OE_new_groups.R

### merged ATAC, multiome RNA new_groups

#.libPaths("/home/kwinkley/02_analysis_software_libraries/seurat4/")
# Updated for Azure cloud_projects

#project_dir <- "/home/wacheung/analysis/TRIGR/merged_scATAC_DA"
project_dir <- "/analysis/cloud_projects/research/TRIGR"
setwd(project_dir)

# source("SRC/TRIGR_mergedCells_OE_new_groups.R")

#THREADS=25

library(tibble)
library(dplyr)
library(Seurat)
library(Signac)
library(sctransform)
library(future)
library(ggplot2)
library(readr)
library(tidyr)
library(stringr)
library(GenomicRanges)
library(EnsDb.Hsapiens.v86)
library(GenomeInfoDb)
library(SingleCellExperiment)
library(ExperimentHub)
#library(celldex)
#library(SeuratDisk)
library(glmGamPoi)
library(future)
library(SingleR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

#minCells=6
#minSamples=0
minCells=10
#minSamples=10
minSamples=5
maxFDR=0.1


subanalysis_name <- 
  paste0("2024_11_20_TRIGR_mergedCells_OE_new_groups_cells", minCells,
         "_samples", minSamples)
#  paste0("2024_11_07_TRIGR_mergedCells_OE_new_groups_cells", minCells,
#         "_samples", minSamples)


options(
  future.globals.maxSize=100*1e9
)

# Output folder for this document
options(knitr.figure_dir = paste0(
  project_dir,
  "/RESULT/",
  subanalysis_name
  )
)

#source("/home/kwinkley/03_scripts/seurat_rmarkdown_scripts/knitr_common.r")
#source("/home/wacheung/analysis/single_cell/Laurel_scPig/SRC/knitr_common.r")
source("SRC/knitr_common.r")

so_ATAC <- faster_readRDS(
  paste(
    "DATA",
#    "2021_04_16_TRIGR_scATAC_scMulti_merging",
    "integrated_ATAC_so.rds",
    sep = "/"
  )
)

so_ATAC$CellType <- ifelse(
  so_ATAC$assay_type == "scMulti", 
  so_ATAC$CellType, 
  so_ATAC$predicted.id
)

new_names <- ifelse(
  so_ATAC$assay_type == "scMulti", 
  rownames(so_ATAC@meta.data), 
  so_ATAC$full_droplet_name
)

so_ATAC <- RenameCells(
  object = so_ATAC, 
  new.names = new_names
)

so <- faster_readRDS(
  paste(
    "DATA",
#    "2021_03_15_TRIGR_scMulti_so_prep",
    "so_integrated_no_batch_var.rds",
    sep = "/"
  )
)

so@meta.data <- faster_readRDS(
  paste(
    "DATA",
    "metadata.rds",
    sep = "/"
  )
)

so@commands <- faster_readRDS(
  paste(
    "DATA",
    "commands.rds",
    sep = "/"
  )
)

so@reductions$umap.wnn <- faster_readRDS(
  paste(
    "DATA",
    "wnn_umap_reduction.rds",
    sep = "/"
  )
)

so@reductions$umap.rna <- faster_readRDS(
  paste(
    "DATA",
    "rna_umap_reduction.rds",
    sep = "/"
  )
)

so@reductions$umap.atac <- faster_readRDS(
  paste(
    "DATA",
    "atac_umap_reduction.rds",
    sep = "/"
  )
)

so@reductions$CITE_seq_umap <- faster_readRDS(
  paste(
    "DATA",
    "CITE_seq_umap_reduction.rds",
    sep = "/"
  )
)

so@assays$ATAC <- faster_readRDS(
  paste(
    "DATA",
    "blood_peaks_ATAC_assay_w_links_cluster_run.rds",
    sep = "/"
  )
)

so@assays$macs2Peaks <- faster_readRDS(
  paste(
    "DATA",
    "macs2_peaks_ATAC_assay_w_links_cluster_run.rds",
    sep = "/"
  )
)


so$CellType <- faster_readRDS(
  paste(
    "DATA",
    "cell_labels.rds",
    sep = "/"
  )
)

so_calendar_age = read_csv("DATA/calendar_age_TRIGR.csv")
so_calendar_age = so_calendar_age %>% dplyr::mutate(dplyr::across(DNA, as.character))

so_ATAC$nCount_macs2Peaks <- colSums(so_ATAC@assays$macs2Peaks@counts)
so_ATAC@assays$macs2Peaks@fragments <- so_ATAC@assays$ATAC@fragments

ignore_cell_types <- c("Unknown #2", "Unknown #1", "Not Assigned")

link_correlations <- read.csv(
  paste(
    "DATA",
    "2021_05_04_TRIGR_scRNA_scATAC_link_correlation",
    "link_external_correlations.csv",
    sep = "/"
  )
)


gencode.genes <- read_tsv("/reference/GRCh38/gencode38/gencode.v38.annotation.genes.bed", col_names=c("gene_chr", "gene_start", "gene_end", "gene", "strand"))

so_ATAC@meta.data$CellType_merged=
  with(so_ATAC@meta.data,
       ifelse(CellType %in% c("CD4 T", "Proliferating T", "GD T"),
              "CD4 T Cells",
              ifelse(CellType %in% c("CD8 Memory Eff T", "CD8 T"),
                     "CD8 T Cells",
                     ifelse(CellType %in% c("B", "Differentiated B"),
                            "B Cells",
                            ifelse(CellType %in% c("NK CD56 Bright", "NK CD56 Dim", "NKT"),
                                   "NK Cells",
                                   ifelse(CellType %in% c("CD14 Mono", "CD16 Mono", "cDC", "pDC"),
                                          "Monoctye DC", "Not Assigned"
                                          ))))))

so@meta.data$CellType_merged=
  with(so@meta.data,
       ifelse(CellType %in% c("CD4 T", "Proliferating T", "GD T"),
              "CD4 T Cells",
              ifelse(CellType %in% c("CD8 Memory Eff T", "CD8 T"),
                     "CD8 T Cells",
                     ifelse(CellType %in% c("B", "Differentiated B"),
                            "B Cells",
                            ifelse(CellType %in% c("NK CD56 Bright", "NK CD56 Dim", "NKT"),
                                   "NK Cells",
                                   ifelse(CellType %in% c("CD14 Mono", "CD16 Mono", "cDC", "pDC"),
                                          "Monoctye DC", "Not Assigned"
                                          ))))))


assignments <- readxl::read_excel(
  paste(
    "DATA",
    "HLA_assignments.xlsx",
    sep = "/"
  )
)

to_assign <- so@meta.data %>%
  tibble::rownames_to_column(var = "full_barcode_name") %>%
  dplyr::left_join(y = assignments, by = "Localid") %>% 
  dplyr::select(full_barcode_name, `HLA genetic effect test group`)

so$TP_HLA_grouping <- data.frame(
  to_assign$`HLA genetic effect test group`,
  row.names = to_assign$full_barcode_name
)

to_assign <- so_ATAC@meta.data %>%
  tibble::rownames_to_column(var = "full_barcode_name") %>%
  dplyr::left_join(y = assignments, by = "Localid") %>% 
  dplyr::select(full_barcode_name, `HLA genetic effect test group`)

so_ATAC$TP_HLA_grouping <- data.frame(
  to_assign$`HLA genetic effect test group`,
  row.names = to_assign$full_barcode_name
)

new_groups=read_tsv("DATA/TRIGR_2024_new_groups.filtered.tsv")

to_assign = so@meta.data %>%
  tibble::rownames_to_column(var = "full_barcode_name") %>%
  dplyr::select(full_barcode_name, Localid) %>%
  dplyr::left_join(new_groups %>% dplyr::select(Localid, new_grouping), by="Localid")

so$new_grouping <- data.frame(
  to_assign$`new_grouping`,
  row.names = to_assign$full_barcode_name
)

to_assign = so_ATAC@meta.data %>%
  tibble::rownames_to_column(var = "full_barcode_name") %>%
  dplyr::select(full_barcode_name, Localid) %>%
  dplyr::left_join(new_groups %>% dplyr::select(Localid, new_grouping), by="Localid")

so_ATAC$new_grouping <- data.frame(
  to_assign$`new_grouping`,
  row.names = to_assign$full_barcode_name
)

##################
# Plots for 2024 paper
#show main clusters
#show superclusters
#overlay case control
#timepoint 1/2/3 as separate graphs
#"new case/control" IAA IAA control,  GADA/GADA control
so@meta.data$tp=so@meta.data$timepoint #weird array size error?
so@meta.data$cond=so@meta.data$condition #weird array size error?

# 2024-01-14 - remove the IAA_lc and IAA_lc_Control

so@meta.data$new_grouping[so@meta.data$new_grouping == "IAA_lc"] = NA
so@meta.data$new_grouping[so@meta.data$new_grouping == "IAA_lc_Control"] = NA

so_ATAC@meta.data$new_grouping[so_ATAC@meta.data$new_grouping == "IAA_lc"] = NA
so_ATAC@meta.data$new_grouping[so_ATAC@meta.data$new_grouping == "IAA_lc_Control"] = NA

pdf("RESULT/paper_2024_12.mergedCells_OE.pdf", width=16, height=9)
MODE="ATAC"
DimPlot(so_ATAC, group.by = "condition") + ggtitle(paste(MODE, "Case/Control status"))
DimPlot(so_ATAC, group.by = "CellType") + ggtitle(paste(MODE, "Cell Type Clustering"))
DimPlot(so_ATAC, group.by = "CellType_merged") + ggtitle(paste(MODE, "Merged Cell Types"))
DimPlot(so_ATAC, group.by = "CellType_merged", split.by = "timepoint") + ggtitle(paste(MODE, "Merged Cell Types, Timepoint"))
DimPlot(so_ATAC, group.by = "CellType_merged", split.by = "new_grouping") + ggtitle(paste(MODE, "Merged Cell Types, New Grouping"))

MODE="snRNA"
DimPlot(so, group.by = "cond") + ggtitle(paste(MODE, "Case/Control status"))
DimPlot(so, group.by = "CellType") + ggtitle(paste(MODE, "Cell Type Clustering"))
DimPlot(so, group.by = "CellType_merged") + ggtitle(paste(MODE, "Merged Cell Types"))
DimPlot(so, group.by = "CellType_merged", split.by = "tp") + ggtitle(paste(MODE, "Merged Cell Types, Timepoint"))
DimPlot(so, group.by = "CellType_merged", split.by = "new_grouping") + ggtitle(paste(MODE, "Merged Cell Types, New Grouping"))

dev.off()

#p1 <- DimPlot(seurat_object, reduction = "umap", group.by = "cell_type")
#p2 <- DimPlot(seurat_object, reduction = "umap", group.by = "sample")
#CombinePlots(plots = list(p1, p2))
