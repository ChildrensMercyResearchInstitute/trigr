## older scRNA new_groupings

# check calendar age
# add no calendar age
#setwd("research/TRIGR/")
#source("SRC/TRIGR_scRNA.diff_oe_new_groups.R", echo=T)

# conda activate signac
#library(SeuratDisk) # needs hdf5
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
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SingleR)
library(future)
#library(glmGamPoi)
# off CRAN (tests fail?) - recommended older version?
# devtools::install_url("https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.7.tar.gz")
library(Matrix.utils)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(parallel)

PREFIX="seurat_analysis_2020-09-28"

# scRNA seurat prep dir
RNA_DIR=".../TRIGR_scRNA_so_prep"

OUT_DIR="RESULT/2024_11_11_scRNA_OE_diff_express_new_grouping"

padj_cutoff <- 0.1

faster_readRDS <- function(file) {
  stopifnot(file.exists(file))
  pigz_in <- pipe(paste0("pigz -dc ", file), "rb")
  on.exit(close(pigz_in))
  readRDS(pigz_in)
}

### RNA PCA
so <- faster_readRDS(file.path(RNA_DIR, "so.rds"))
sc <- faster_readRDS(file.path(RNA_DIR, "singler.rds"))
umap_df <- faster_readRDS(file.path(RNA_DIR, "umap_df.rds"))

# "cluster" is the factor version of seurat_clusters?
umap_df %<>%
  dplyr::rename(UMAP1=UMAP_1,
         UMAP2=UMAP_2) %>%
  dplyr::mutate(cluster = factor(seurat_clusters, levels=sort(unique(as.integer(seurat_clusters)))))

so_calendar_age = read_csv("DATA/calendar_age_TRIGR.csv")
so_calendar_age = so_calendar_age %>% dplyr::mutate(dplyr::across(DNA, as.character))
so_calendar_age$glue_timepoint = paste0("timepoint_", so_calendar_age$Timepoint)

gencode.genes <- read_tsv("DATA/gencode.v32.annotation.genes.bed",
                             col_names=c("gene_chr", "gene_start", "gene_end",
                                         "gene", "strand"))


#new_groups=read_tsv("DATA/TRIGR_2024_new_groups.filtered.tsv")
#meta=faster_readRDS("DATA/metadata.rds")
#write_tsv(meta |> dplyr::select(Reg_Code, Localid, genotype_ID) |> unique() |> left_join(new_groups), "DATA/TRIGR_2024_new_groups.filtered.geno.tsv")

new_groups=read_tsv("DATA/TRIGR_2024_new_groups.filtered.geno.tsv",
  col_types = list( genotype_ID = col_character()) )


to_assign = so@meta.data %>%
  tibble::rownames_to_column(var = "full_barcode_name") %>%
  dplyr::select(full_barcode_name, genotype_ID=genotype) %>%
  dplyr::left_join(new_groups %>% dplyr::select(genotype_ID, new_grouping), by="genotype_ID")

so$new_grouping <- data.frame(
  to_assign$`new_grouping`,
  row.names = to_assign$full_barcode_name
)

######

OUT_PREFIX=paste0(PREFIX, ".new_groups")

dir.create(OUT_DIR)

assay_name = "SCT"
DefaultAssay(so)=assay_name

ignore_cell_types <- c("Unknown #2", "Unknown #1", "Not Assigned")
minCells=10
#minSamples=10
minSamples=5
maxFDR=0.1

stopifnot(all.equal(sc$HPCA$cluster, sc$BP$cluster))

totals <- umap_df %>% 
  group_by(cluster) %>%
  summarize(total_cells = n(),
            .groups='drop')

cell_types <- data.frame(cluster = sc$HPCA$cluster,
                        HPCA    = sc$HPCA$labels,
                        BP      = sc$BP$labels) 


celltypes_tomi=read_tsv("DATA/TRIGR_scRNA_cell_types.tomi.tsv")
cell_types = cell_types |> left_join(celltypes_tomi)

# need to transform so that the rownames don't get jumbled
# do this AFTER joining with celltypes_tomi since the join can reorder the rows
rownames(cell_types)=as.character(cell_types$cluster)

# Tomi reconcicled BluePrint and HPCA annotations to CellType_merged
so@meta.data$CellType_merged=cell_types[so@meta.data$seurat_clusters,]$CellType_merged
so@meta.data$CellType_merged_cluster=cell_types[so@meta.data$seurat_clusters,]$cluster

# sanity check that the merge was successful
print(head(so@meta.data[, c("CellType_merged", "seurat_clusters", "CellType_merged_cluster")]))
stopifnot(all.equal(as.integer(so@meta.data$seurat_clusters), as.integer(so@meta.data$CellType_merged_cluster)))

cell_types_to_test <- so@meta.data %>%
    dplyr::filter(!CellType_merged %in% ignore_cell_types) %>%
    dplyr::pull(CellType_merged) %>%
    unique() %>% sort()

groups_to_test <- so@meta.data %>%
    dplyr::filter(!CellType_merged %in% ignore_cell_types) %>%
    dplyr::pull(group) %>%
    unique() %>% sort()

timepoints_to_test <- so@meta.data |>
    dplyr::filter(!CellType_merged %in% ignore_cell_types) %>%
    dplyr::pull(timepoint) %>%
    unique() %>% sort()

#variable naming diffs TRIGR merge vs scRNA
#genotype = genotype_ID
#group = condition

total_reads=data.frame(cell_reads=colSums(so[[assay_name]]@counts)) %>%
    tibble::rownames_to_column(var = "full_barcode") %>%
    dplyr::left_join(so@meta.data %>%
                     tibble::rownames_to_column(var = "full_barcode") %>%
                     dplyr::select(full_barcode, genotype, timepoint)) %>%
    dplyr::group_by(genotype, timepoint) %>%
    dplyr::summarise(total_reads=sum(cell_reads))

total_celltypes=data.frame(cell_reads=colSums(so[[assay_name]]@counts)) %>%
        tibble::rownames_to_column(var = "full_barcode") %>%
        dplyr::left_join(so@meta.data %>%
                         tibble::rownames_to_column(var = "full_barcode") %>%
                         dplyr::select(full_barcode, CellType_merged, timepoint)) %>%
        dplyr::group_by(CellType_merged, timepoint) %>%
        dplyr::summarise(total_reads=sum(cell_reads), total_cells=n())
total_cluster=data.frame(cell_reads=colSums(so[[assay_name]]@counts)) %>%
        tibble::rownames_to_column(var = "full_barcode") %>%
        dplyr::left_join(so@meta.data %>%
                         tibble::rownames_to_column(var = "full_barcode") %>%
                         dplyr::select(full_barcode, seurat_clusters, CellType_merged, timepoint)) %>%
        dplyr::group_by(seurat_clusters, CellType_merged, timepoint) %>%
        dplyr::summarise(total_reads=sum(cell_reads), total_cells=n())

cell_types |> 
  left_join(totals, by="cluster") %>%
  arrange(as.integer(cluster)) %>%
  dplyr::select(cluster, total_cells, everything()) 

### These stats are unchanged
#write_tsv(cell_types |> 
#  left_join(totals, by="cluster") %>%
#  arrange(as.integer(cluster)) %>%
#  dplyr::select(cluster, total_cells, everything()) , 
#  file.path(OUT_DIR, paste0(OUT_PREFIX, ".scRNA_OE_diff_express_case_control.total_celltypes.tsv")))

#write_tsv(total_celltypes, 
#  file.path(OUT_DIR, paste0(OUT_PREFIX, ".scRNA_OE_diff_express_case_control.total_celltypes_readstp.tsv")))
#write_tsv(total_cluster, 
#  file.path(OUT_DIR, paste0(OUT_PREFIX, ".scRNA_OE_diff_express_case_control.total_clusters_readstp.tsv")))


#### 2024 Nov comparisons
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

#subset_ident_1 <- "Case"
#subset_ident_2 <- "Control"
#assay_name = "SCT"

diff_OE <- function(assay_name, subset_ident_1, subset_ident_2) {
  all_results_tmp <- data.frame()

  for (cells in cell_types_to_test) {
    for (tp in timepoints_to_test) {
      print(c("TESTING", cells, tp))
      cells1 <- so@meta.data %>%
        tibble::rownames_to_column(var = "full_barcode") %>%
        dplyr::filter(CellType_merged == cells) %>%
        dplyr::filter(timepoint == tp) %>%
        dplyr::filter(new_grouping == subset_ident_1) %>%
        dplyr::group_by(genotype) %>%
        dplyr::mutate(cells_per_indiv = dplyr::n()) %>%
        dplyr::filter(cells_per_indiv >= minCells) %>%
        dplyr::pull(full_barcode)
      indiv_group_1 <- so@meta.data %>%
        tibble::rownames_to_column(var = "full_barcode") %>%
        dplyr::filter(full_barcode %in% cells1) %>%
        dplyr::pull(genotype) %>%
        unique() %>%
        length()
      cells2 <- so@meta.data %>%
        tibble::rownames_to_column(var = "full_barcode") %>%
        dplyr::filter(CellType_merged == cells) %>%
        dplyr::filter(timepoint == tp) %>%
        dplyr::filter(new_grouping == subset_ident_2) %>%
        dplyr::group_by(genotype) %>%
        dplyr::mutate(cells_per_indiv = dplyr::n()) %>%
        dplyr::filter(cells_per_indiv >= minCells) %>%
        dplyr::pull(full_barcode)
      indiv_group_2 <- so@meta.data %>%
        tibble::rownames_to_column(var = "full_barcode") %>%
        dplyr::filter(full_barcode %in% cells2) %>%
        dplyr::pull(genotype) %>%
        unique() %>%
        length()
  #    if((length(cells1) > 0) & (length(cells2) > 0)) {
      #print(c("g1", indiv_group_1, "g2", indiv_group_2)) ## DEBUG
      if(indiv_group_1 >= minSamples & indiv_group_2 >= minSamples) {
        de_results = data.frame(num_reads1=rowSums(so[[assay_name]]@counts[,cells1])) %>% tibble::rownames_to_column(var="gene")
        de_results$num_reads2=rowSums(so[[assay_name]]@counts[,cells2])
        
        total_reads1=sum(de_results$num_reads1)
        total_reads2=sum(de_results$num_reads2)
        de_results=de_results %>%
          dplyr::mutate(rest_reads1=total_reads1-num_reads1,
                rest_reads2=total_reads2-num_reads2,
                timepoint=tp, CellType=cells,
                indiv_group1=indiv_group_1,
                indiv_group2=indiv_group_2) %>%
          dplyr::group_by(gene) %>%
          dplyr::mutate(FC=(num_reads2/(num_reads2+rest_reads2))/(num_reads1/(num_reads1+rest_reads1)), p_val=fisher.test(matrix(c(num_reads1, rest_reads1, num_reads2, rest_reads2), ncol=2))$p.value)

        all_results_tmp <- rbind(
            all_results_tmp, de_results
          )
        print(dim(de_results))
        print(dim(all_results_tmp)) ## DEBUG
        
      } else {
        print(c(indiv_group_1, indiv_group_2))
      }
    }
  }
  print(dim(all_results_tmp)) ## DEBUG

  # Try all timepoints merged
  print("ALL TP merged")
  tp=NA
  for (cells in cell_types_to_test) {
    cells1 <- so@meta.data %>%
        tibble::rownames_to_column(var = "full_barcode") %>%
        dplyr::filter(CellType_merged == cells) %>%
        dplyr::filter(new_grouping == subset_ident_1) %>%
        dplyr::group_by(genotype) %>%
        dplyr::mutate(cells_per_indiv = dplyr::n()) %>%
        dplyr::filter(cells_per_indiv >= minCells) %>%
        dplyr::pull(full_barcode)
      indiv_group_1 <- so@meta.data %>%
        tibble::rownames_to_column(var = "full_barcode") %>%
        dplyr::filter(full_barcode %in% cells1) %>%
        dplyr::pull(genotype) %>%
        unique() %>%
        length()
      cells2 <- so@meta.data %>%
        tibble::rownames_to_column(var = "full_barcode") %>%
        dplyr::filter(CellType_merged == cells) %>%
        dplyr::filter(new_grouping == subset_ident_2) %>%
        dplyr::group_by(genotype) %>%
        dplyr::mutate(cells_per_indiv = dplyr::n()) %>%
        dplyr::filter(cells_per_indiv >= minCells) %>%
        dplyr::pull(full_barcode)
      indiv_group_2 <- so@meta.data %>%
        tibble::rownames_to_column(var = "full_barcode") %>%
        dplyr::filter(full_barcode %in% cells2) %>%
        dplyr::pull(genotype) %>%
        unique() %>%
        length()
  #    if((length(cells1) > 0) & (length(cells2) > 0)) {
      if(indiv_group_1 >= minSamples & indiv_group_2 >= minSamples) {
        de_results = data.frame(num_reads1=rowSums(so[[assay_name]]@counts[,cells1])) %>% tibble::rownames_to_column(var="gene")
        de_results$num_reads2=rowSums(so[[assay_name]]@counts[,cells2])
        
        total_reads1=sum(de_results$num_reads1)
        total_reads2=sum(de_results$num_reads2)
        de_results=de_results %>%
          dplyr::mutate(rest_reads1=total_reads1-num_reads1,
                rest_reads2=total_reads2-num_reads2,
                timepoint=tp, CellType=cells,
                indiv_group1=indiv_group_1,
                indiv_group2=indiv_group_2) %>%
          dplyr::group_by(gene) %>%
          dplyr::mutate(FC=(num_reads2/(num_reads2+rest_reads2))/(num_reads1/(num_reads1+rest_reads1)), p_val=fisher.test(matrix(c(num_reads1, rest_reads1, num_reads2, rest_reads2), ncol=2))$p.value)

        all_results_tmp <- rbind(
            all_results_tmp, de_results
          )
      } else {
        print(c(indiv_group_1, indiv_group_2))
      }
  }
  print(dim(all_results_tmp)) ## DEBUG

  print("FDR")
  colnames(all_results_tmp)
  all_results_tmp <- all_results_tmp %>%
    dplyr::group_by(CellType, timepoint) %>%
    dplyr::mutate(fdr_p_val = p.adjust(p_val, method = "fdr"))  %>%
    dplyr::arrange(fdr_p_val)

  all_results_tmp_filtered <- 
    all_results_tmp %>% 
    dplyr::filter(fdr_p_val <= maxFDR) %>%
    dplyr::arrange(fdr_p_val)
  #%>%
  #  dplyr::filter(indiv_group1 >= minSamples & indiv_group2 >= minSamples) %>%

  write_csv(
    x = all_results_tmp_filtered %>%
      dplyr::left_join(gencode.genes) %>%
      dplyr::mutate(fdr_p_val=signif(fdr_p_val),
                    p_val=signif(p_val),
                    FC=signif(FC)),
    file = file.path(OUT_DIR, paste0(OUT_PREFIX, ".scRNA_OE_diff_express_",subset_ident_1,"_", subset_ident_2, ".csv"))
    )

  write_csv(
    x = all_results_tmp %>%
      dplyr::left_join(gencode.genes) %>%
      dplyr::mutate(fdr_p_val=signif(fdr_p_val),
                    p_val=signif(p_val),
                    FC=signif(FC)),
    file = file.path(OUT_DIR, paste0(OUT_PREFIX, ".scRNA_OE_diff_express_",subset_ident_1,"_", subset_ident_2, ".full.csv.gz"))
    )

  return(all_results_tmp)
}

result_GADA_IAA_sc = diff_OE("SCT",  "GADA", "IAA_sc")
result_GADA_GADA_control = diff_OE("SCT",  "GADA", "GADA_Control")
result_IAA_sc_IAA_sc_Control = diff_OE("SCT",  "IAA_sc", "IAA_sc_Control")
