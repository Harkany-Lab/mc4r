


# Load tidyverse infrastructure packages
suppressPackageStartupMessages({
    library(here)
    library(tidyverse)
    library(magrittr)
    library(future)
})

if (!requireNamespace("glmGamPoi", quietly = TRUE)) BiocManager::install("glmGamPoi")
if (!requireNamespace("Nebulosa", quietly = TRUE)) BiocManager::install("Nebulosa")
if (!requireNamespace("limma", quietly = TRUE)) BiocManager::install('limma')
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) BiocManager::install("ComplexHeatmap")
if (!requireNamespace("schex", quietly = TRUE)) BiocManager::install("schex")
# Load packages for scRNA-seq analysis and visualisation
suppressPackageStartupMessages({
    library(Seurat)
    library(SeuratWrappers)
    library(SeuratDisk)
    library(conos)
    library(harmony)
    library(sctransform)
    library(UpSetR)
    library(patchwork)
    library(Nebulosa)
    library(schex)
})

# set seed
set.seed(seed = reseed)

n_bins <- 40

##### Load in data #####
samples_df <- read_tsv(here::here(data_dir, "samples.tsv"))
rar2020.srt.pub <- readr::read_rds(file.path(data_dir, "subset_cca_to_velo.rds"))
rar2020.srt.pub %<>% UpdateSeuratObject()
colnames(rar2020.srt.pub@reductions$umap@cell.embeddings) <- c('UMAP_1', 'UMAP_2')
# Consistent colours and clusters names
rar2020.srt.pub$wtree <- factor(rar2020.srt.pub$wtree,
                                levels  = 1:45,
                                ordered = TRUE)
colours <- readr::read_lines(here(data_dir, "colours_wtree.tsv"))
clrlev <- levels(rar2020.srt.pub$wtree)
names(colours) <- clrlev
rar2020.clusters.neuro.invert <-
    as.character(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 38, 42, 45))
# Dev. Stages
rar2020.ages.all <- c('E15', 'E17', 'P00', 'P02', 'P10', 'P23')
rar2020.ages.postnat <-                         c('P10', 'P23')

# Subset
# stage
rar2020.srt.pub$orig.ident <-
    rar2020.srt.pub %>%
    colnames() %>%
    str_split(pattern = ':', simplify = TRUE) %>% .[, 1] %>%
    plyr::mapvalues(x    = .,
                    from = samples_df$fullname,
                    to   = samples_df$sample)
rar2020.srt.pub$samples <- factor(rar2020.srt.pub$orig.ident)
rar2020.srt.pub$samples %>% summary()

readr::write_rds(
    rar2020.srt.pub,
    file = here::here(tables_dir,
                      '2021-05-18_001-086-446-605-669_mm10-hypothalamus_srt4-obj-umap.rds'))
SaveH5Seurat(object = rar2020.srt.pub,
             filename = here::here(data_dir, "2021-05-18_001-086-446-605-669_mm10-hypothalamus_srt4-obj-umap.h5seurat"))
rar2020.srt.pub$age %>% summary()
rar2020.srt.pub$postnatal <-
    rar2020.srt.pub$age %in% rar2020.ages.postnat
rar2020.srt.pub <- subset(rar2020.srt.pub,
                          subset = postnatal == TRUE)
Idents(rar2020.srt.pub) <- "wtree"
readr::write_rds(
    rar2020.srt.pub,
    file = here::here(tables_dir,
                      '2021-05-18_001-086-446-605-669_adult_mm10-hypothalamus_srt4-obj-umap.rds'))
SaveH5Seurat(object = rar2020.srt.pub,
             filename = here::here(data_dir, "2021-05-18_001-086-446-605-669_adult_mm10-hypothalamus_srt4-obj-umap.h5seurat"))

rar2020.srt.pvn <-
    subsetSrt(dat.srt = rar2020.srt.pub,
              srt.name = "Neuronal clusters of PVN and DMH P10/23 mice",
              idents = c(15, 24, 26, 31, 43),
              invert = FALSE)
rar2020.srt.pvn %<>% RenameIdents(object = ., `43` = "PVN-mneOXY")
rar2020.srt.pvn %<>% RenameIdents(object = ., `26` = "PVN-mneVAS")
rar2020.srt.pvn %<>% RenameIdents(object = ., `31` = "PVN-pneSst")
rar2020.srt.pvn %<>% RenameIdents(object = ., `24` = "PVN-pneCrhTrh")
rar2020.srt.pvn %<>% RenameIdents(object = ., `15` = "DMH-Trh")
rar2020.srt.pvn$ident <- Idents(rar2020.srt.pvn)
pvn.cols <- colours[c(15, 24, 26, 31, 43)]
names(pvn.cols) %<>%
    plyr::mapvalues(x = .,
                    from = c(15, 24, 26, 31, 43),
                    to = c("PVN-mneOXY",
                           "PVN-mneVAS",
                           "PVN-pneSst",
                           "PVN-pneCrhTrh",
                           "DMH-Trh"))

rar2020.srt.pvn.bck <- rar2020.srt.pvn
readr::write_rds(
    rar2020.srt.pvn.bck,
    file = here::here(tables_dir,
                      '2021-05-18_001-086-446-605-669_adult_mm10-hypothalamus-pvn_alk-thinness-paper-subset-srt4-obj-umap.rds'))
SaveH5Seurat(object = rar2020.srt.pvn.bck,
             filename = here::here(data_dir, "2021-05-18_001-086-446-605-669_adult_mm10-hypothalamus-pvn_alk-thinness-paper-subset-srt4-obj-umap.h5seurat"))

## FIXME integration doesn't work (everything below hasn't checked)
## Probably very first step is wrong and I should change default assay to RNA
################################################################################
################# HARMONY INTEGRATION ##########################################
################################################################################
hpth.harmony <- NormalizeData(rar2020.srt.pvn) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
hpth.harmony <- RunHarmony(hpth.harmony, group.by.vars = 'samples')
hpth.harmony <- RunUMAP(hpth.harmony, reduction = 'harmony', dims = 1:30)
hpth.harmony <- FindNeighbors(hpth.harmony, reduction = 'harmony', dims = 1:30) %>% FindClusters()
DimPlot(hpth.harmony, group.by = 'wtree', cols = colours) / DimPlot(hpth.harmony, group.by = 'samples')
hpth.harmony <- make_hexbin(hpth.harmony, nbins = n_bins,
                            dimension_reduction = "UMAP")
plot_hexbin_density(hpth.harmony)
################################################################################
readr::write_rds(
    hpth.harmony,
    file = here::here(tables_dir,
                      '2021-05-18_001-086-446-605-669_adult_mm10-hypothalamus-pvn_alk-thinness-paper-subset-srt4-obj-umap-harmony.rds'))
SaveH5Seurat(object = hpth.harmony,
             filename = here::here(data_dir, "2021-05-18_001-086-446-605-669_adult_mm10-hypothalamus-pvn_alk-thinness-paper-subset-srt4-obj-umap-harmony.h5seurat"))

################################################################################
################################################################################
################# CONOS INTEGRATION ############################################
################################################################################
hpth.conos.panel <- SplitObject(hpth, split.by = 'samples')
for (i in 1:length(hpth.conos.panel)) {
    hpth.conos.panel[[i]] <- NormalizeData(hpth.conos.panel[[i]]) %>% FindVariableFeatures() %>% ScaleData() %>%
        RunPCA(verbose = FALSE)
}
hpth.conos.con <- Conos$new(hpth.conos.panel)
hpth.conos.con$buildGraph(k = 15, k.self = 5, space = 'PCA', ncomps = 30, n.odgenes = 2000, matching.method = 'mNN',
                          metric = 'angular', score.component.variance = TRUE, verbose = TRUE)
hpth.conos.con$findCommunities()
hpth.conos.con$embedGraph()
hpth.conos <- as.Seurat(hpth.conos.con)
names(hpth.conos@reductions) %<>% tolower()
DimPlot(hpth.harmony, reduction = 'largevis', group.by = 'wtree', cols = colours) / DimPlot(hpth.harmony, reduction = 'largevis', group.by = 'samples')
hpth.conos <- make_hexbin(hpth.conos, nbins = n_bins,
                          dimension_reduction = "largeVis")
plot_hexbin_density(hpth.conos)
################################################################################
readr::write_rds(
    hpth.conos,
    file = here::here(tables_dir,
                      '2021-05-18_001-086-446-605-669_adult_mm10-hypothalamus-pvn_alk-thinness-paper-subset-srt4-obj-conos.rds'))
SaveH5Seurat(object = hpth.conos,
             filename = here::here(data_dir, "2021-05-18_001-086-446-605-669_adult_mm10-hypothalamus-pvn_alk-thinness-paper-subset-srt4-obj-conos.h5seurat"))

################################################################################
################################################################################

rar2020.srt.pvn %<>% subset(x = ., subset = UMAP_1 > 5)
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

