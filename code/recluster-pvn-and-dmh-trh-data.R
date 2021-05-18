## It breaks on FindIntegrationAnchors (line 61) in both sequential and parallel mode
## I've tried to decrease number of anchors for filtering and scaling anchors
## TODO complete debug or/and search for description of error code below
# Error in irlba(A = mat3, nv = num.cc) :
# BLAS/LAPACK routine 'DLASCL' gave error code -4


# Load tidyverse infrastructure packages
suppressPackageStartupMessages({
    library(here)
    library(tidyverse)
    library(magrittr)
    library(future)
})

if (!requireNamespace("Nebulosa", quietly = TRUE)) BiocManager::install("Nebulosa")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install('limma')
# Load packages for scRNA-seq analysis and visualisation
suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratWrappers)
    library(SeuratDisk)
    library(sctransform)
    library(glmGamPoi)
    library(UpSetR)
    library(patchwork)
    library(Nebulosa)
})

# set seed
set.seed(seed = reseed)

##### Load in data #####
samples_df <- read_tsv(here::here(data_dir, "samples.tsv"))
rar2020.srt.pvn.dmh <- readr::read_rds(file.path(tables_dir, "2021-05-13_001-086-446-605-669_adult_mm10-hypothalamus-pvn_alk-thinness-paper-subset-srt4-obj-umap.rds"))
colours <- readr::read_lines(here(data_dir, "colours_wtree.tsv"))
clrlev <- 1:45
names(colours) <- clrlev
pvn.cols <- colours[c(15, 24, 26, 31, 43)]
names(pvn.cols) %<>%
  plyr::mapvalues(x = .,
                  from = c(15, 24, 26, 31, 43),
                  to = c("pneTRH",
                         "pneCRH",
                         "mneVAS",
                         "pneSS",
                         "mneOXY"))


rar2020.srt.pd.list <-
    SplitObject(rar2020.srt.pvn.dmh, split.by = "samples")
for (i in 1:length(rar2020.srt.pd.list)) {
  rar2020.srt.pd.list[[i]] %<>%
        SCTransform(seed.use = reseed,
                    return.only.var.genes = FALSE,
                    verbose = FALSE)
}

features <- SelectIntegrationFeatures(object.list = rar2020.srt.pd.list, nfeatures = 3000)
rar2020.srt.pd.list <- PrepSCTIntegration(object.list = rar2020.srt.pd.list, anchor.features = features)
plan('sequential')
rar2020.srt.pd.anchors <-
  FindIntegrationAnchors(object.list          = rar2020.srt.pd.list,
                         k.filter             = 50,
                         k.score              = 10,
                         normalization.method = "SCT",
                         anchor.features      = features)
rar2020.pvn.dmh.combined.sct <- IntegrateData(anchorset = rar2020.srt.pd.anchors, normalization.method = "SCT")

## Code below is from derive-data.R script and has not been tested here yet
DefaultAssay(rar2020.srt.pvn) <- "RNA"
rar2020.srt.pvn <-
    SCTransform(rar2020.srt.pvn,
                variable.features.n = 3000,
                vars.to.regress = c("log_umi_per_gene"),
                return.only.var.genes = FALSE,
                seed.use = reseed,
                verbose = FALSE)
readr::write_rds(
    rar2020.srt.pvn,
    file = here::here(tables_dir,
                      '2021-05-13_001-086-446-605-669_adult_mm10-hypothalamus-pvn_bottom-lh-excl-srt4-obj-sct-umap.rds'))
SaveH5Seurat(object = rar2020.srt.pvn,
             filename = here::here(data_dir, "2021-05-13_001-086-446-605-669_adult_mm10-hypothalamus-pvn_bottom-lh-excl-srt4-obj-sct-umap.h5seurat"))

Idents(rar2020.srt.pvn) <- "ident"
all_markers_pvn_wtree_final <- FindAllMarkers(rar2020.srt.pvn,
                                              assay = "SCT",
                                              test.use = "wilcox",
                                              logfc.threshold = 0.05,
                                              min.pct = 0.05,
                                              random.seed = reseed,
                                              return.thresh = 0.01)

all_markers_pvn_wtree_final %>%
    group_by(cluster) %>%
    filter(p_val_adj < 0.01) %>%
    slice_max(n = 7, order_by = avg_log2FC) %>%
    print(., n = 35)
readr::write_csv(
    all_markers_pvn_wtree_final,
    file = here::here(tables_dir,
                      '2021-05-13_001-086-446-605-669_adult_mm10-hypothalamus-pvn_deg-mrk-gns-ctypes-by-sct-wilcox.csv'))

markers_wtree_final <-
    all_markers_pvn_wtree_final %>%
    dplyr::filter(gene %in% c(gene_int)) %>%
    group_by(cluster) %>%
    top_n(n = 12, wt = pct.1) %>%
    top_n(n = 7, wt = avg_log2FC) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(pct.1)) %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>%
    dplyr::arrange(cluster) %>%
    .$gene

srt_mcrs <-
    GetAssayData(rar2020.srt.pvn, 'data', 'SCT') %>%
    as.data.frame() %>%
    .[unique(c(markers_wtree_final)) %>%
          .[. %in% row.names(rar2020.srt.pvn@assays$SCT@data)], ] %>%
    as.data.frame() %>%
    t()
readr::write_csv(
    srt_mcrs %>% as_tibble(rownames = "cell_names"),
    file = here::here(tables_dir,
                      '2021-05-13_001-086-446-605-669_adult_mm10-hypothalamus-pvn_expr-mtx-cells-by-mrk-gns-int.csv'))

minim <- srt_mcrs %>% as_tibble() %>% dplyr::summarise_each(funs = function(x) quantile(x, .3))

plot_data <-
    (srt_mcrs > as.double(minim)) %>%
    as_tibble() %>%
    mutate_all(as.numeric) %>%
    dplyr::select(-Avp, -Oxt)
readr::write_csv(
    plot_data,
    file = here::here(tables_dir,
                      '2021-05-13_001-086-446-605-669_adult_mm10-hypothalamus-pvn_ncells-long-by-mrk-gns-int.csv'))

plot_data_filt <-
    (srt_mcrs > as.double(minim)) %>%
    as_tibble() %>% filter(Mc4r > 0) %>%
    mutate_all(as.numeric) %>%
    dplyr::select(-Avp, -Oxt, -Mc4r)
readr::write_csv(
    plot_data_filt,
    file = here::here(tables_dir,
                      '2021-05-13_001-086-446-605-669_adult_mm10-hypothalamus-pvn-mc4r-pos_ncells-long-by-mrk-gns-int.csv'))

plot_data2 <-
    (srt_mcrs > as.double(minim)) %>%
    as_tibble() %>%
    mutate_all(as.numeric) %>%
    mutate('Cell_type' = Idents(rar2020.srt.pvn),
           'Cell_name' = colnames(rar2020.srt.pvn)) %>%
    as_tibble()

plot_data3 <-
    plot_data2 %>%
    mutate(flat = 1) %>%
    spread(Cell_type, flat, fill = 0)
readr::write_csv(
    plot_data3,
    file = here::here(tables_dir,
                      '2021-05-13_001-086-446-605-669_adult_mm10-hypothalamus-pvn_ncells-long-by-mrk-gns-int-ctypes.csv'))

rar2020.srt.pub <-
    readr::read_rds(here::here(tables_dir,
                               '2021-05-13_001-086-446-605-669_mm10-hypothalamus_srt4-obj-umap.rds'))
rar2020.srt.pub <-
    readr::read_rds(here::here(tables_dir,
                               '2021-05-13_001-086-446-605-669_adult_mm10-hypothalamus_srt4-obj-umap.rds'))
rar2020.srt.pvn.bck <-
    readr::read_rds(here::here(tables_dir,
                               '2021-05-13_001-086-446-605-669_adult_mm10-hypothalamus-pvn_alk-thinness-paper-subset-srt4-obj-umap.rds'))
rar2020.srt.pvn <-
    readr::read_rds(here::here(tables_dir,
                               '2021-05-13_001-086-446-605-669_adult_mm10-hypothalamus-pvn_bottom-lh-excl-srt4-obj-sct-umap.rds'))
all_markers_pvn_wtree_final <-
    readr::read_csv(here::here(tables_dir,
                               '2021-05-13_001-086-446-605-669_adult_mm10-hypothalamus-pvn_deg-mrk-gns-ctypes-by-sct-wilcox.csv'))
srt_mcrs <-
    readr::read_csv(here::here(tables_dir,
                               '2021-05-13_001-086-446-605-669_adult_mm10-hypothalamus-pvn_expr-mtx-cells-by-mrk-gns-int.csv')) %>%
    column_to_rownames(var = "cell_names")
plot_data <-
    readr::read_csv(here::here(tables_dir,
                               '2021-05-13_001-086-446-605-669_adult_mm10-hypothalamus-pvn_ncells-long-by-mrk-gns-int.csv'))
plot_data_filt <-
    readr::read_csv(here::here(tables_dir,
                               '2021-05-13_001-086-446-605-669_adult_mm10-hypothalamus-pvn-mc4r-pos_ncells-long-by-mrk-gns-int.csv'))
plot_data3 <-
    readr::read_csv(here::here(tables_dir,
                               '2021-05-13_001-086-446-605-669_adult_mm10-hypothalamus-pvn_ncells-long-by-mrk-gns-int-ctypes.csv'))

#TODO attach table to convert between different IDs

