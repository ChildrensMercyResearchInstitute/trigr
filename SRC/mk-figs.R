#setwd("/analysis/cloud_projects/research/TRIGR")

# Additional paper figures
# Tomi notes:
# The independent scRNA/ATAC you found would go into a supplementary figure and we use the multiome umap/wnn for figure 1.
# Generate UMAP for scRNA and scATAC similar to existing UMAP for scMulti (ATAC & RNA)
# generate multiome coverage plot

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
library(ggforce)

ATAC_DIR="DATA/"
RNA_DIR="DATA/"
RNA_DIR2="DATA/seurat_analysis_ribodepletion_2020-10-30"
RNA_DIR3="DATA/seurat_analysis_2020-09-28"
OUT_DIR="RESULT/"

ignore_cell_types <- c("Unknown #2", "Unknown #1", "Not Assigned")

faster_readRDS <- function(file) {
        stopifnot(file.exists(file))
        pigz_in <- pipe(paste0("pigz -dc ", file), "rb")
        on.exit(close(pigz_in))
        readRDS(pigz_in)
}

faster_saveRDS <- function(object, file, cores=8) {
  pigz_out <- pipe(paste0("pigz -p", cores, " > ", file), "wb")
  on.exit(close(pigz_out))
  saveRDS(object, pigz_out)
}


# ATAC so
INFILE=file.path(file.path(ATAC_DIR, "integrated_ATAC_so.rds"))
MODE="ATAC"
so <- faster_readRDS(INFILE)
so$CellType <- faster_readRDS(file.path(ATAC_DIR, "so_ATAC_projected_cell_labels.rds"))


plt1 <- DimPlot(
  so, group.by = "CellType") +
  NoLegend() + 
  ggtitle(MODE)

ggsave(
  plot = plt1,
  filename = file.path("RESULT", "paper", paste0(MODE, "_umap.jpeg")),
  units = "in",
  width = 10,
  height = 8
)

# "old" scRNA
MODE="RNA2020-09-28"
INFILE=file.path(RNA_DIR3, "so.rds")
so <- faster_readRDS(INFILE)
plt1 <- DimPlot(
  so, group.by = "seurat_clusters") +
  NoLegend() + 
  ggtitle(MODE)

ggsave(
  plot = plt1,
  filename = file.path("RESULT", "paper", paste0(MODE, "_umap.jpeg")),
  units = "in",
  width = 10,
  height = 8
)

## multiome coverage plot
## BROKEN this doesn't work as it needs direct access to the multiome fragment files
## Use dnanexus heatmap generation code
## (tried updatePath but that doesn't work?)
## would need like download/relink them
MODE="scMulti"
INFILE=file.path(RNA_DIR, "so_integrated_no_batch_var.rds")
so <- faster_readRDS(INFILE)
so@meta.data <- faster_readRDS(file.path(RNA_DIR, "metadata.rds"))
so$CellType <- faster_readRDS(file.path(RNA_DIR, "cell_labels.rds"))
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
cell_types_to_test <- so@meta.data %>%
    dplyr::filter(!CellType_merged %in% ignore_cell_types) %>%
    dplyr::pull(CellType_merged) %>%
    unique() %>% sort()

so@assays$ATAC <- faster_readRDS(file.path(RNA_DIR, "blood_peaks_ATAC_assay_w_links_cluster_run.rds"))

# fix the Fragments
new_fragments=list()
for(i in Fragments(so@assays$ATAC)) {
    frag_path=i@path %>%
        str_replace("/gpfs/data/analysis/projects/xenon/cellranger-arc/",
            "/analysis/projects/xenon/cellranger-arc/") 
    # for the bulk ATAC aggregation
    frag_path=frag_path %>%
      str_replace("/gpfs/data/analysis/projects/xenon/", "/analysis/projects/xenon/")

    cell_ids=i@cells
    print(frag_path)
    frags <- CreateFragmentObject(
        path = frag_path,
        cells = cell_ids,
        validate.fragments = FALSE)
    new_fragments=append(new_fragments, frags)
}

Fragments(so@assays$ATAC) <- NULL # clear the old fragments
Fragments(so@assays$ATAC) <- new_fragments # add the new fragments

DefaultAssay(so)="ATAC"
faster_saveRDS(so, file.path(RNA_DIR, "so_integrated_no_batch_var_fixazpath.rds"))
so=faster_readRDS(file.path(RNA_DIR, "so_integrated_no_batch_var_fixazpath.rds"))

plot_region="chr17-35671491-36312878"
gene_names=c("CCL3L1", "CCL3", "CCL5", "CCL4L2")

# ATAC needs the fragments
cov_plot=CoveragePlot(
      object = so,
      region = plot_region,
      features = gene_names,
      group.by="CellType_merged",
      annotation = TRUE,
      peaks = FALSE,
      tile = TRUE,
      tile.size= 5000,
      tile.cells=20,
      window = 5000,
      links = TRUE
  ) + theme_gray(base_size = 22) # font size
pdf("RESULT/paper/coveragePlot.scMulti.CCL3L1.express.pdf", width=8.5, height=6)
print(cov_plot) 
dev.off()

########## FeaturePlot the RNA Mono

INFILE="RESULT/recluster_mono/TRIGR_recluster_mono.RNA.res0.8.so.rds"
so=faster_readRDS(INFILE)
pdf("RESULT/paper/FeaturePlot.recluster_mono.RNA.res0.8.genes.pdf")
for(f in c("TNFAIP3", "NFKB1", "IL1B", "TNFAIP6", "IL1A")) {
  print(FeaturePlot(so, f, alpha=0.75, raster=TRUE))
  print(FeaturePlot(so, f, alpha=0.75, split.by="condition", raster=TRUE))
}
dev.off()

##### gene DE

so$seurat_clusters.condition.tp=paste(so$seurat_clusters, so$condition, so$timepoint, sep=".")
Idents(so)="seurat_clusters.condition.tp"

top_markers=NULL
for(i in unique(so@meta.data$seurat_clusters)) {
  print(i)
  for(tp in unique(so$timepoint)) {
    print(tp)
    top_markers=top_markers |>
      bind_rows(FindMarkers(so, ident.1=paste(i, "case", tp, sep="."), 
                    ident.2=paste(i, "control", tp, sep=".")) |>
                head(n=100) |>
                mutate(group=i, timepoint=tp) |>
                rownames_to_column("feature"))
  }
}
write_tsv(top_markers |> arrange(group,timepoint,p_val), "RESULT/paper/TRIGR_recluster_mono.RNA.res0.8.case_control.topMarkers.tsv")
########### gene DE o/e analysis

maxFDR=0.1
#minCells=10
#minSamples=10
#OUTFILE=paste0("RESULT/paper/TRIGR_recluster_mono.RNA.res0.8.case_control.oe_results.fdr.", maxFDR, ".tsv")

minCells=0
minSamples=0
OUTFILE=paste0("RESULT/paper/TRIGR_recluster_mono.RNA.res0.8.case_control.oe_results.fdr.", maxFDR, ".nofilter.tsv")

assay_name="SCT"

all_results_case_control_RNA_OE <- data.frame()
for(i in unique(so@meta.data$seurat_clusters)) {
  print(i)
  for(tp in unique(so$timepoint)) {
    print(tp)
    cells_case <- so@meta.data %>%
      tibble::rownames_to_column(var = "full_barcode") %>%
      dplyr::filter(seurat_clusters == i) %>%
      dplyr::filter(timepoint == tp) %>%
      dplyr::filter(condition == "case") %>%
      dplyr::group_by(genotype_ID) %>%
      dplyr::mutate(cells_per_indiv = dplyr::n()) %>%
      dplyr::filter(cells_per_indiv >= minCells) %>%
      dplyr::pull(full_barcode)
    n_case <- so@meta.data %>%
      tibble::rownames_to_column(var = "full_barcode") %>%
      dplyr::filter(full_barcode %in% cells_case) %>%
      dplyr::pull(genotype_ID) %>%
      unique() %>%
      length()

    cells_control <- so@meta.data %>%
      tibble::rownames_to_column(var = "full_barcode") %>%
      dplyr::filter(seurat_clusters == i) %>%
      dplyr::filter(timepoint == tp) %>%
      dplyr::filter(condition == "control") %>%
      dplyr::group_by(genotype_ID) %>%
      dplyr::mutate(cells_per_indiv = dplyr::n()) %>%
      dplyr::filter(cells_per_indiv >= minCells) %>%
      dplyr::pull(full_barcode)
    n_control <- so@meta.data %>%
      tibble::rownames_to_column(var = "full_barcode") %>%
      dplyr::filter(full_barcode %in% cells_control) %>%
      dplyr::pull(genotype_ID) %>%
      unique() %>%
      length()

    # basic 5 read case/control for the fisher test
    if(n_case >= minSamples && n_control >= minSamples && length(cells_case) >= 5 && length(cells_control) >= 5) {

      de_results = data.frame(num_reads1=rowSums(so[[assay_name]]@counts[,cells_case])) %>% tibble::rownames_to_column(var="gene")
      de_results$num_reads2=rowSums(so[[assay_name]]@counts[,cells_control])
      
      total_reads1=sum(de_results$num_reads1)
      total_reads2=sum(de_results$num_reads2)
      de_results=de_results %>%
        dplyr::mutate(rest_reads1=total_reads1-num_reads1,
               rest_reads2=total_reads2-num_reads2,
               timepoint=tp, cluster_id=i,
               n_case=n_case,
               n_control=n_control) %>%
        dplyr::group_by(gene) %>%
        dplyr::mutate(FC=(num_reads2/(num_reads2+rest_reads2))/(num_reads1/(num_reads1+rest_reads1)), p_val=fisher.test(matrix(c(num_reads1, rest_reads1, num_reads2, rest_reads2), ncol=2))$p.value)

      all_results_case_control_RNA_OE <<- rbind(
          all_results_case_control_RNA_OE, de_results
        )
      
    }
  }
}

# all timepoints merged
tp=NA
for(i in unique(so@meta.data$seurat_clusters)) {
  print(i)
  cells_case <- so@meta.data %>%
    tibble::rownames_to_column(var = "full_barcode") %>%
    dplyr::filter(seurat_clusters == i) %>%
    dplyr::filter(condition == "case") %>%
    dplyr::group_by(genotype_ID) %>%
    dplyr::mutate(cells_per_indiv = dplyr::n()) %>%
    dplyr::filter(cells_per_indiv >= minCells) %>%
    dplyr::pull(full_barcode)
  n_case <- so@meta.data %>%
    tibble::rownames_to_column(var = "full_barcode") %>%
    dplyr::filter(full_barcode %in% cells_case) %>%
    dplyr::pull(genotype_ID) %>%
    unique() %>%
    length()

  cells_control <- so@meta.data %>%
    tibble::rownames_to_column(var = "full_barcode") %>%
    dplyr::filter(seurat_clusters == i) %>%
    dplyr::filter(condition == "control") %>%
    dplyr::group_by(genotype_ID) %>%
    dplyr::mutate(cells_per_indiv = dplyr::n()) %>%
    dplyr::filter(cells_per_indiv >= minCells) %>%
    dplyr::pull(full_barcode)
  n_control <- so@meta.data %>%
    tibble::rownames_to_column(var = "full_barcode") %>%
    dplyr::filter(full_barcode %in% cells_control) %>%
    dplyr::pull(genotype_ID) %>%
    unique() %>%
    length()

  if(n_case >= minSamples && n_control >= minSamples && length(cells_case) >= 5 && length(cells_control) >= 5) {

    de_results = data.frame(num_reads1=rowSums(so[[assay_name]]@counts[,cells_case])) %>% tibble::rownames_to_column(var="gene")
    de_results$num_reads2=rowSums(so[[assay_name]]@counts[,cells_control])
    
    total_reads1=sum(de_results$num_reads1)
    total_reads2=sum(de_results$num_reads2)
    de_results=de_results %>%
      dplyr::mutate(rest_reads1=total_reads1-num_reads1,
              rest_reads2=total_reads2-num_reads2,
              timepoint=tp, cluster_id=i,
              n_case=n_case,
              n_control=n_control) %>%
      dplyr::group_by(gene) %>%
      dplyr::mutate(FC=(num_reads2/(num_reads2+rest_reads2))/(num_reads1/(num_reads1+rest_reads1)), p_val=fisher.test(matrix(c(num_reads1, rest_reads1, num_reads2, rest_reads2), ncol=2))$p.value)

    all_results_case_control_RNA_OE <<- rbind(
        all_results_case_control_RNA_OE, de_results
      )
    
  }
}

all_results_case_control_RNA_OE <- all_results_case_control_RNA_OE %>%
  dplyr::group_by(cluster_id, timepoint) %>%
  dplyr::mutate(fdr_p_val = p.adjust(p_val, method = "fdr"))

write_tsv(all_results_case_control_RNA_OE |> 
            dplyr::filter(fdr_p_val <= maxFDR) |> 
            mutate(across(contains("p_val"), signif, 4)) |>
            mutate(across(contains("FC"), signif, 4)) |>
            arrange(cluster_id,timepoint,fdr_p_val), 
          OUTFILE)

print(OUTFILE)
##############

# pull out cells from cluster 11, 4 and 0 
# with high IL1B expression (>5 and >7.5 in your feature plot) 
# list those cells counts in case and control per timepoint. 

# FeaturePlot uses the normalised data from the @data slot
# https://github.com/satijalab/seurat/issues/1485

# default assay is "integratedRNA" NOT RNA
assay_name="integratedRNA"
#assay_name="RNA"
gene="IL1B"
gene_cellcounts=NULL
for(i in unique(so@meta.data$seurat_clusters)) {
  for(tp in unique(so$timepoint)) {
      print(c(i,tp))
      cells_case <- so@meta.data %>%
        tibble::rownames_to_column(var = "full_barcode") %>%
        dplyr::filter(seurat_clusters == i) %>%
        dplyr::filter(timepoint == tp) %>%
        dplyr::filter(condition == "case") %>%
        dplyr::group_by(genotype_ID) %>%
        dplyr::pull(full_barcode)

      n_case=length(cells_case)
      n_case0=sum(so[[assay_name]]@data[gene,cells_case] > 0)
      n_case5=sum(so[[assay_name]]@data[gene,cells_case] > 5)
      n_case7.5=sum(so[[assay_name]]@data[gene,cells_case] > 7.5)

      cells_control <- so@meta.data %>%
        tibble::rownames_to_column(var = "full_barcode") %>%
        dplyr::filter(seurat_clusters == i) %>%
        dplyr::filter(timepoint == tp) %>%
        dplyr::filter(condition == "control") %>%
        dplyr::group_by(genotype_ID) %>%
        dplyr::pull(full_barcode)

      n_control=length(cells_control)
      n_control0=sum(so[[assay_name]]@data[gene,cells_control] > 0)
      n_control5=sum(so[[assay_name]]@data[gene,cells_control] > 5)
      n_control7.5=sum(so[[assay_name]]@data[gene,cells_control] > 7.5)

    gene_cellcounts<- gene_cellcounts |> rbind(data.frame(
      cluster_id=i, timepoint=tp, 
      n_case=n_case, n_case0=n_case0, n_case5=n_case5, n_case7.5=n_case7.5,
      n_control=n_control, n_control0=n_control0, n_control5=n_control5, n_control7.5=n_control7.5
    ))
  }
}

write_tsv(gene_cellcounts |> arrange(cluster_id,timepoint), 
    paste0("RESULT/paper/TRIGR_recluster_mono.", assay_name, ".res0.8.gene_cellcounts.", gene, ".tsv"))



#####
INFILE="RESULT/recluster_mono/TRIGR_recluster_mono.macs2Peaks.res0.8.so.rds"
so_ATAC=faster_readRDS(INFILE)
atac_peaks=rownames(so_ATAC)
closest_genes=ClosestFeature(so_ATAC, atac_peaks)
pdf(paste0("RESULT/paper/FeaturePlot.recluster_mono.macs2Peaks.res0.8.peaks.pdf"))
for(f in c("TNFAIP3", "NFKB1", "IL1B", "TNFAIP6", "IL1A")) {
  for (p in closest_genes |> dplyr::filter(gene_name == f & distance==0) |> pull(query_region) |> unique()) {
    print(FeaturePlot(so_ATAC, p, alpha=0.75, raster=TRUE)+ggtitle(paste(p, "(", f,")")))
    print(FeaturePlot(so_ATAC, p, alpha=0.75, split.by="condition", raster=TRUE)+ggtitle(paste(p, "(", f,")")))
  }
}
dev.off()

##############
# Find Expressed genes
# columns are CellType_Timepoint
MODE="scMulti"
INFILE=file.path(RNA_DIR, "so_integrated_no_batch_var.rds")
so <- faster_readRDS(INFILE)
so@meta.data <- faster_readRDS(file.path(RNA_DIR, "metadata.rds"))
so$CellType <- faster_readRDS(file.path(RNA_DIR, "cell_labels.rds"))
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
cell_types_to_test <- so@meta.data %>%
    dplyr::filter(!CellType_merged %in% ignore_cell_types) %>%
    dplyr::pull(CellType_merged) %>%
    unique() %>% sort()

so.expressed=AggregateExpression(so, group.by=c("CellType_merged", "timepoint"))

so.expressed.hist=NULL
for(i in c(1,5,10,20)) {
  so.expressed.hist= so.expressed.hist |> bind_rows(c(RNAminReads=i, colSums(so.expressed$RNA >= i))
  )
}
write_tsv(so.expressed.hist, "RESULT/paper/TRIGR_scMulti_RNA_expressed_hist.tsv")

# "old" scRNA
MODE="RNA2020-09-28"
INFILE=file.path(RNA_DIR3, "so.rds")
so_rna <- faster_readRDS(INFILE)

# get Tomi merged celltypes
cell_types <- data.frame(cluster = sc$HPCA$cluster,
                     HPCA    = sc$HPCA$labels,
                     BP      = sc$BP$labels) 

celltypes_tomi=read_tsv("DATA/TRIGR_scRNA_cell_types.tomi.tsv")

cell_types = cell_types |> left_join(celltypes_tomi)
so_rna@meta.data$CellType_merged=cell_types[so_rna@meta.data$seurat_clusters,]$CellType_merged

so_rna.expressed=AggregateExpression(so_rna, group.by=c("CellType_merged", "timepoint"))

so_rna.expressed.hist=NULL
for(i in c(1,5,10,20)) {
  so_rna.expressed.hist= so_rna.expressed.hist |> bind_rows(c(RNAminReads=i, colSums(so_rna.expressed$RNA >= i))
  )
}
write_tsv(so_rna.expressed.hist, "RESULT/paper/TRIGR_RNA2020-09-28_RNA_expressed_hist.tsv")
