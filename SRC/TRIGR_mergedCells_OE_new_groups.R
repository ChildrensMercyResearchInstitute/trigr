### merged ATAC, multiome RNA new_groups

# Updated for Azure cloud_projects

# top level output directory
project_dir <- "..."

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
library(celldex)
library(SeuratDisk)
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

# merged ATAC
seurat_scatac_merged <- ".../TRIGR_scATAC_scMulti_merging"
# multiome
seurat_scmulti_dir <- ".../TRIGR_scMulti_so_prep"

# multiome cell types
seurat_scmulti_cells <- ".../TRIGR_scMulti_cell_type_curation"

# link correlation result dir
trigr_cor <- ".../TRIGR_scRNA_scATAC_link_correlation"

subanalysis_name <- 
  paste0("2024_11_20_TRIGR_mergedCells_OE_new_groups_cells", minCells,
         "_samples", minSamples)


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

source("SRC/knitr_common.r")

so_ATAC <- faster_readRDS(
  paste(
    seurat_scatac_merged,
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
    seurat_scmulti_dir,
    "so_integrated_no_batch_var.rds",
    sep = "/"
  )
)

so@meta.data <- faster_readRDS(
  paste(
    seurat_scmulti_dir,
    "metadata.rds",
    sep = "/"
  )
)

so@commands <- faster_readRDS(
  paste(
    seurat_scmulti_dir,
    "commands.rds",
    sep = "/"
  )
)

so@reductions$umap.wnn <- faster_readRDS(
  paste(
    seurat_scmulti_dir,
    "wnn_umap_reduction.rds",
    sep = "/"
  )
)

so@reductions$umap.rna <- faster_readRDS(
  paste(
    seurat_scmulti_dir,
    "rna_umap_reduction.rds",
    sep = "/"
  )
)

so@reductions$umap.atac <- faster_readRDS(
  paste(
    seurat_scmulti_dir,
    "atac_umap_reduction.rds",
    sep = "/"
  )
)

so@reductions$CITE_seq_umap <- faster_readRDS(
  paste(
    seurat_scmulti_dir,
    "CITE_seq_umap_reduction.rds",
    sep = "/"
  )
)

so@assays$ATAC <- faster_readRDS(
  paste(
    seurat_scmulti_dir,
    "blood_peaks_ATAC_assay_w_links_cluster_run.rds",
    sep = "/"
  )
)

so@assays$macs2Peaks <- faster_readRDS(
  paste(
    seurat_scmulti_dir,
    "macs2_peaks_ATAC_assay_w_links_cluster_run.rds",
    sep = "/"
  )
)


so$CellType <- faster_readRDS(
  paste(
    seurat_scmulti_cells,
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
    trigr_cor,
    "link_external_correlations.csv",
    sep = "/"
  )
)


gencode.genes <- read_tsv("DATA/gencode.v38.annotation.genes.bed", col_names=c("gene_chr", "gene_start", "gene_end", "gene", "strand"))

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

### helper function
# keep as global to view tmp results?
all_results_tmp = data.frame()
subset_ident_1 <- NA
subset_ident_2 <- NA

diff_OE <- function(test_assay,
                    test_id1, test_id2, test_var, filter_condition=FALSE) {
    # reset
    if(test_assay == "ATAC") {
      assay_name = "macs2Peaks"
      so_name = "so_ATAC"
    } else if (test_assay == "RNA") {
      assay_name="SCT"
      so_name = "so"
    } else {
      print(paste("UNKNOWN ASSAY", test_assay))
      return(NA)
    }
  
    all_results_tmp = data.frame()

    subset_ident_1 <- test_id1
    subset_ident_2 <- test_id2
    
    print(c("Testing", test_id1, test_id2, test_var, filter_condition))
    
    # Use CellType_merged
    cell_types_to_test <- get(so_name)@meta.data %>%
      dplyr::filter(!CellType_merged %in% ignore_cell_types) %>%
      dplyr::pull(CellType_merged) %>%
      unique()
    
    # Fisher test based on:
    #     num_reads,
    ##        total_reads_in_group,
    ##        total_reads_in_peak - num_reads,
    ##        total_reads - total_reads_in_group
    
    all_results_tmp <- data.frame()
    
    for (cells in cell_types_to_test) {
      for (tp in (1:3)) {
        print(c("Testing", cells, tp))
        #print(dim(all_results_tmp))
        cells1 <- get(so_name)@meta.data %>%
          tibble::rownames_to_column(var = "full_barcode") %>%
          dplyr::filter(CellType_merged == cells) %>%
          dplyr::filter(timepoint == tp) %>%
          dplyr::filter(.data[[test_var]] == subset_ident_1) %>%
          dplyr::filter(isFALSE(filter_condition) | condition == filter_condition) %>%
          dplyr::group_by(genotype_ID) %>%
          dplyr::mutate(cells_per_indiv = dplyr::n()) %>%
          dplyr::filter(cells_per_indiv >= minCells) %>%
          dplyr::pull(full_barcode)
        indiv_group_1 <- get(so_name)@meta.data %>%
          tibble::rownames_to_column(var = "full_barcode") %>%
          dplyr::filter(full_barcode %in% cells1) %>%
          dplyr::pull(genotype_ID) %>%
          unique() %>%
          length()
        cells2 <- get(so_name)@meta.data %>%
          tibble::rownames_to_column(var = "full_barcode") %>%
          dplyr::filter(CellType_merged == cells) %>%
          dplyr::filter(timepoint == tp) %>%
          dplyr::filter(.data[[test_var]] == subset_ident_2) %>%
          dplyr::filter(isFALSE(filter_condition) | condition == filter_condition) %>%
          dplyr::group_by(genotype_ID) %>%
          dplyr::mutate(cells_per_indiv = dplyr::n()) %>%
          dplyr::filter(cells_per_indiv >= minCells) %>%
          dplyr::pull(full_barcode)
        indiv_group_2 <- get(so_name)@meta.data %>%
          tibble::rownames_to_column(var = "full_barcode") %>%
          dplyr::filter(full_barcode %in% cells2) %>%
          dplyr::pull(genotype_ID) %>%
          unique() %>%
          length()
    #    if((length(cells1) > 0) & (length(cells2) > 0)) {
        #print(c(length(cells1), length(cells2)))
        #print(c(indiv_group_1, indiv_group_2))
        if(indiv_group_1 >= minSamples & indiv_group_2 >= minSamples) {
          de_results = data.frame(num_reads1=rowSums(get(so_name)[[assay_name]]@counts[,cells1])) %>% tibble::rownames_to_column(var="peak")
          de_results$num_reads2=rowSums(get(so_name)[[assay_name]]@counts[,cells2])
          
          total_reads1=sum(de_results$num_reads1)
          total_reads2=sum(de_results$num_reads2)
          de_results=de_results %>%
            dplyr::mutate(rest_reads1=total_reads1-num_reads1,
                   rest_reads2=total_reads2-num_reads2,
                   timepoint=tp, CellType=cells,
                   indiv_group1=indiv_group_1,
                   indiv_group2=indiv_group_2) %>%
            dplyr::group_by(peak) %>%
            dplyr::mutate(FC=(num_reads2/(num_reads2+rest_reads2))/(num_reads1/(num_reads1+rest_reads1)), p_val=fisher.test(matrix(c(num_reads1, rest_reads1, num_reads2, rest_reads2), ncol=2))$p.value)
    
          all_results_tmp <- rbind(
              all_results_tmp, de_results
            )
          
        }
      }
    }
    
    # Try all timepoints merged
    tp=NA
    for (cells in cell_types_to_test) {
       cells1 <- get(so_name)@meta.data %>%
          tibble::rownames_to_column(var = "full_barcode") %>%
          dplyr::filter(CellType_merged == cells) %>%
          dplyr::filter(.data[[test_var]] == subset_ident_1) %>%
          dplyr::filter(isFALSE(filter_condition) | condition == filter_condition) %>%
          dplyr::group_by(genotype_ID) %>%
          dplyr::mutate(cells_per_indiv = dplyr::n()) %>%
          dplyr::filter(cells_per_indiv >= minCells) %>%
          dplyr::pull(full_barcode)
        indiv_group_1 <- get(so_name)@meta.data %>%
          tibble::rownames_to_column(var = "full_barcode") %>%
          dplyr::filter(full_barcode %in% cells1) %>%
          dplyr::pull(genotype_ID) %>%
          unique() %>%
          length()
        cells2 <- get(so_name)@meta.data %>%
          tibble::rownames_to_column(var = "full_barcode") %>%
          dplyr::filter(CellType_merged == cells) %>%
          dplyr::filter(.data[[test_var]] == subset_ident_2) %>%
          dplyr::filter(isFALSE(filter_condition) | condition == filter_condition) %>%
          dplyr::group_by(genotype_ID) %>%
          dplyr::mutate(cells_per_indiv = dplyr::n()) %>%
          dplyr::filter(cells_per_indiv >= minCells) %>%
          dplyr::pull(full_barcode)
        indiv_group_2 <- get(so_name)@meta.data %>%
          tibble::rownames_to_column(var = "full_barcode") %>%
          dplyr::filter(full_barcode %in% cells2) %>%
          dplyr::pull(genotype_ID) %>%
          unique() %>%
          length()
    #    if((length(cells1) > 0) & (length(cells2) > 0)) {
        if(indiv_group_1 >= minSamples & indiv_group_2 >= minSamples) {
          de_results = data.frame(num_reads1=rowSums(get(so_name)[[assay_name]]@counts[,cells1])) %>% tibble::rownames_to_column(var="peak")
          de_results$num_reads2=rowSums(get(so_name)[[assay_name]]@counts[,cells2])
          
          total_reads1=sum(de_results$num_reads1)
          total_reads2=sum(de_results$num_reads2)
          de_results=de_results %>%
            dplyr::mutate(rest_reads1=total_reads1-num_reads1,
                   rest_reads2=total_reads2-num_reads2,
                   timepoint=tp, CellType=cells,
                   indiv_group1=indiv_group_1,
                   indiv_group2=indiv_group_2) %>%
            dplyr::group_by(peak) %>%
            dplyr::mutate(FC=(num_reads2/(num_reads2+rest_reads2))/(num_reads1/(num_reads1+rest_reads1)), p_val=fisher.test(matrix(c(num_reads1, rest_reads1, num_reads2, rest_reads2), ncol=2))$p.value)
    
          #print(dim(de_results))
          all_results_tmp <- rbind(
              all_results_tmp, de_results
            )
          
        }
    }
    
    print(dim(all_results_tmp))
    #print("FDR")
    all_results_tmp <- all_results_tmp %>%
      dplyr::group_by(CellType, timepoint) %>%
      dplyr::mutate(fdr_p_val = p.adjust(p_val, method = "fdr"))

    #print("Writing output")
    if(isFALSE(filter_condition)) {
      OUTFILE=paste0(figure_path(), "all_", 
                     test_assay, "_OE_diff_access_", 
                     test_var, "_", test_id1, "_", test_id2,
                     ".csv")
      OUTFILE_FULL=paste0(figure_path(), "all_", 
                     test_assay, "_OE_diff_access_", 
                     test_var, "_", test_id1, "_", test_id2,
                     ".full.csv.gz")

    } else {
      OUTFILE=paste0(figure_path(), "all_",
                     test_assay, "_OE_diff_access_", 
                      test_var, "_", test_id1, "_", test_id2,
                      ".only_", filter_condition, ".csv")
      OUTFILE_FULL=paste0(figure_path(), "all_",
                     test_assay, "_OE_diff_access_", 
                      test_var, "_", test_id1, "_", test_id2,
                      ".only_", filter_condition, ".full.csv.gz")
    }
    print(paste("OUTPUT FILE:", OUTFILE))
    print(paste("FULL OUTPUT FILE:", OUTFILE_FULL))
    if(test_assay == "ATAC") {
      write.csv(
        x = all_results_tmp %>% 
          dplyr::filter(fdr_p_val <= maxFDR) %>%
          dplyr::arrange(fdr_p_val) %>%
          dplyr::mutate(fdr_p_val=signif(fdr_p_val),
                        p_val=signif(p_val),
                        FC=signif(FC)) %>%
          dplyr::left_join(link_correlations),
        file = OUTFILE,
        row.names = FALSE, quote=FALSE
      )

      write.csv(
        x = all_results_tmp %>% 
          dplyr::arrange(fdr_p_val) %>%
          dplyr::mutate(fdr_p_val=signif(fdr_p_val),
                        p_val=signif(p_val),
                        FC=signif(FC)) %>%
          dplyr::left_join(link_correlations),
        file = OUTFILE_FULL,
        row.names = FALSE, quote=FALSE
      )

    } else if (test_assay == "RNA") {
      all_results_tmp <- all_results_tmp %>%
        dplyr::rename(gene=peak)
      write.csv(
        x = all_results_tmp %>% 
          dplyr::filter(fdr_p_val <= maxFDR) %>%
          dplyr::arrange(fdr_p_val) %>%
          dplyr::mutate(fdr_p_val=signif(fdr_p_val),
                        p_val=signif(p_val),
                        FC=signif(FC)) %>%
          dplyr::left_join(gencode.genes),
        file = OUTFILE,
        row.names = FALSE, quote=FALSE
      )

      write.csv(
        x = all_results_tmp %>% 
          dplyr::arrange(fdr_p_val) %>%
          dplyr::mutate(fdr_p_val=signif(fdr_p_val),
                        p_val=signif(p_val),
                        FC=signif(FC)) %>%
          dplyr::left_join(gencode.genes),
        file = OUTFILE_FULL,
        row.names = FALSE, quote=FALSE
      )

    }
    return(all_results_tmp)
}


####
#Analyses 1:
#      3 GADA lc
#      8 GADA sc
#Vs.
#	21 IAA sc
 
#Analyses 2:
#      3 GADA lc
#      8 GADA sc
#Vs.
#	11 paired controls for the above samples
	 
#Analyses 3:
#	21 IAA sc
#	Vs.
#	21 paired controls for the above IAA sc samples

#	In the spreadsheet each case is followed by its’ age matched control

all_results_GADA_GADA_control_ATAC_OE <- diff_OE("ATAC", 
                                       "GADA",
                                       "GADA_Control",
                                       "new_grouping")


all_results_GADA_IAA_sc_ATAC_OE <- diff_OE("ATAC", 
                                       "GADA",
                                       "IAA_sc",
                                       "new_grouping")

all_results_GADA_IAA_sc_RNA_OE <- diff_OE("RNA", 
                                       "GADA",
                                       "IAA_sc",
                                       "new_grouping")


all_results_GADA_GADA_control_RNA_OE <- diff_OE("RNA", 
                                       "GADA",
                                       "GADA_Control",
                                       "new_grouping")


all_results_IAA_sc_IAA_sc_Control_ATAC_OE <- diff_OE("ATAC", 
                                       "IAA_sc",
                                       "IAA_sc_Control",
                                       "new_grouping")


all_results_IAA_sc_IAA_sc_Control_RNA_OE <- diff_OE("RNA", 
                                       "IAA_sc",
                                       "IAA_sc_Control",
                                       "new_grouping")

