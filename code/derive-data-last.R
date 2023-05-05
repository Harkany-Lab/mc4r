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

# Set paths
src_dir    <- here::here('code')
data_dir   <- here::here('data')
output_dir <- here::here('output')
plots_dir  <- here::here(output_dir, 'figures')
tables_dir <- here::here(output_dir, 'tables')
source(here::here(src_dir, 'main.R'))
source(here::here(src_dir, 'genes.R'))

# set seed
reseed <- 42
set.seed(seed = reseed)

##### Load in data #####
samples_df <- read_tsv(here::here(data_dir, "samples.tsv"))
rar2020.srt.pub <- readr::read_rds(file.path(data_dir, "oldCCA_nae_srt.rds"))
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
rar2020.srt.pub$age <-
    plyr::mapvalues(x    = rar2020.srt.pub$orig.ident,
                    from = samples_df$sample,
                    to   = samples_df$age)
Idents(rar2020.srt.pub) <-
    factor(rar2020.srt.pub$age,
           levels  = rar2020.ages.all,
           ordered = TRUE)
rar2020.srt.pub$age <- Idents(rar2020.srt.pub)

rar2020.srt.pub$age %>% summary()
rar2020.srt.pub$postnatal <-
    rar2020.srt.pub$age %in% rar2020.ages.postnat
rar2020.srt.pub <- subset(rar2020.srt.pub,
                          subset = postnatal == TRUE)
Idents(rar2020.srt.pub) <- "wtree"

rar2020.srt.pvn <-
    subsetSrt(dat.srt = rar2020.srt.pub,
              srt.name = "Neuronal clusters of PVN P10/23 mice",
              idents = c(15, 24, 26, 31, 43),
              invert = FALSE)

pvn.cols <- colours[c(15, 24, 26, 31, 43)]
names(pvn.cols) %<>%
    plyr::mapvalues(x = .,
                    from = c(15, 24, 26, 31, 43),
                    to = c("pneTRH",
                           "pneCRH",
                           "mneVAS",
                           "pneSS",
                           "mneOXY"))

rar2020.srt.pvn %<>% subset(x = ., subset = UMAP_1 > 5)
DefaultAssay(rar2020.srt.pvn) <- "RNA"
## gene pairs correlations in pvn
# rna corrected only
p_corrs <- list(
    FeatureScatter(rar2020.srt.pvn, feature1 = "Alk", feature2 = "Mc4r", slot = "data"),
    FeatureScatter(rar2020.srt.pvn, feature1 = "Crh",  feature2 = "Alk", slot = "data"),
    FeatureScatter(rar2020.srt.pvn, feature1 = "Trh",  feature2 = "Alk", slot = "data"),
    FeatureScatter(rar2020.srt.pvn, feature1 = "Crh", feature2 = "Mc4r", slot = "data"),
    FeatureScatter(rar2020.srt.pvn, feature1 = "Trh", feature2 = "Mc4r", slot = "data"),
    FeatureScatter(rar2020.srt.pvn, feature1 = "Crh",  feature2 = "Trh", slot = "data")
)
n_corrs <- list(
    "pvn-def-rna-data-Alk-Mc4r",
    "pvn-def-rna-data-Crh-Alk" ,
    "pvn-def-rna-data-Trh-Alk" ,
    "pvn-def-rna-data-Crh-Mc4r",
    "pvn-def-rna-data-Trh-Mc4r",
    "pvn-def-rna-data-Crh-Trh"
)

purrr::walk2(n_corrs, p_corrs, sPlot, type = "corr-plt")

rar2020.srt.pvn <-
    SCTransform(rar2020.srt.pvn,
                variable.features.n = 3000,
                vars.to.regress = c("log_umi_per_gene"),
                return.only.var.genes = FALSE,
                seed.use = reseed,
                verbose = FALSE)

Idents(rar2020.srt.pvn) <- "wtree"
## gene pairs correlations in pvn
# sct corrected only
p_corrs <- list(
    FeatureScatter(rar2020.srt.pvn, feature1 = "Alk", feature2 = "Mc4r", slot = "data"),
    FeatureScatter(rar2020.srt.pvn, feature1 = "Crh",  feature2 = "Alk", slot = "data"),
    FeatureScatter(rar2020.srt.pvn, feature1 = "Trh",  feature2 = "Alk", slot = "data"),
    FeatureScatter(rar2020.srt.pvn, feature1 = "Crh", feature2 = "Mc4r", slot = "data"),
    FeatureScatter(rar2020.srt.pvn, feature1 = "Trh", feature2 = "Mc4r", slot = "data"),
    FeatureScatter(rar2020.srt.pvn, feature1 = "Crh",  feature2 = "Trh", slot = "data")
)
n_corrs <- list(
    "pvn-def-sct-data-Alk-Mc4r",
    "pvn-def-sct-data-Crh-Alk" ,
    "pvn-def-sct-data-Trh-Alk" ,
    "pvn-def-sct-data-Crh-Mc4r",
    "pvn-def-sct-data-Trh-Mc4r",
    "pvn-def-sct-data-Crh-Trh"
)

purrr::walk2(n_corrs, p_corrs, sPlot, type = "corr-plt")

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
                      '2021-09-24_001-086-446-605-669_adult_mm10-hypothalamus-pvn_deg-mrk-gns-ctypes-by-sct-wilcox.csv'))

markers_wtree_final <-
    all_markers_pvn_wtree_final %>%
    dplyr::filter(gene %in% c(gene_int)) %>%
    group_by(cluster) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(pct.1)) %>%
    dplyr::distinct(gene, .keep_all = TRUE) %>%
    dplyr::arrange(cluster) %>%
    .$gene

rar2020.srt.pub.pne <-
    rar2020.srt.pvn %>%
    subset(idents = c(24))
DefaultAssay(rar2020.srt.pub.pne) <- "RNA"
Idents(rar2020.srt.pub.pne) <- "wtree"

## gene pairs correlations in 24
# rna corrected only
p_corrs <- list(
    FeatureScatter(rar2020.srt.pub.pne, feature1 = "Alk", feature2 = "Mc4r", slot = "data"),
    FeatureScatter(rar2020.srt.pub.pne, feature1 = "Crh",  feature2 = "Alk", slot = "data"),
    FeatureScatter(rar2020.srt.pub.pne, feature1 = "Trh",  feature2 = "Alk", slot = "data"),
    FeatureScatter(rar2020.srt.pub.pne, feature1 = "Crh", feature2 = "Mc4r", slot = "data"),
    FeatureScatter(rar2020.srt.pub.pne, feature1 = "Trh", feature2 = "Mc4r", slot = "data"),
    FeatureScatter(rar2020.srt.pub.pne, feature1 = "Crh",  feature2 = "Trh", slot = "data")
)
n_corrs <- list(
    "24-def-rna-data-Alk-Mc4r",
    "24-def-rna-data-Crh-Alk" ,
    "24-def-rna-data-Trh-Alk" ,
    "24-def-rna-data-Crh-Mc4r",
    "24-def-rna-data-Trh-Mc4r",
    "24-def-rna-data-Crh-Trh"
)

purrr::walk2(n_corrs, p_corrs, sPlot, type = "corr-plt")

rar2020.srt.pub.pne <-
    SCTransform(rar2020.srt.pub.pne,
                variable.features.n = 3000,
                vars.to.regress = c("log_umi_per_gene"),
                return.only.var.genes = FALSE,
                seed.use = reseed,
                verbose = FALSE)

srt_mcrs_filtered_ages <-
    GetAssayData(rar2020.srt.pub.pne, 'data', 'SCT') %>%
    # as.data.frame() %>%
    # .[unique(c(markers_wtree_final)) %>%
    #       .[. %in% row.names(rar2020.srt.pub.pne@assays$SCT@data)], ] %>%
    as.data.frame() %>%
    t()
readr::write_csv(
    srt_mcrs_filtered_ages %>% as_tibble(rownames = "cell_names"),
    file = here::here(tables_dir,
                      '2021-09-24_001-086-446-605-669_adult_mm10-hypothalamus-pvn-pne24_expr-mtx-cells-by-mrk-gns-int.csv'))

filt_age_pvn_sum <-
    colSums(srt_mcrs_filtered_ages) %>%
    .[. > 5] %>%
    names()
srt_mcrs_filtered_ages %<>% .[, filt_age_pvn_sum]

minim_pvn_age_filt <-
    srt_mcrs_filtered_ages %>%
    as_tibble() %>%
    select(all_of(filt_age_pvn_sum)) %>%
    summarise(across(.fns = function(x) quantile(x, .1)))

minim_pvn_age_filt2 <-
    (srt_mcrs_filtered_ages > as.double(minim_pvn_age_filt)) %>%
    as_tibble() %>%
    mutate_all(as.numeric) %>%
    summarise(across(.fns = function(x) sum(x)/length(x))) %>%
    as_vector() %>%
    .[.['Mc4r'] <= .]

srt_mcrs_filtered_ages %<>% .[, names(minim_pvn_age_filt2)]
pvn_age_corr_genes  <- corrr::correlate(srt_mcrs_filtered_ages)

mcrs_age_corr_genes <-
    pvn_age_corr_genes %>%
    arrange(-Mc4r) %>%
    select(term, Mc4r, Trh)


markers_wtree_final_trh <- union(
    all_markers_pvn_wtree_final %>%
        dplyr::filter(cluster %in% c('15'),
                      p_val <= 1e-8) %>%
        mutate(pct.dif = pct.1 - pct.2) %>%
        arrange(-pct.dif) %>%
        filter(pct.dif > .1,
               pct.2 >= .5,
               avg_log2FC >= .5,
               gene %in% big_genes) %>%
        arrange(-avg_log2FC) %>%
        .$gene,
    all_markers_pvn_wtree_final %>%
        dplyr::filter(cluster %in% c('15'),
                      p_val_adj <= 1e-5) %>%
        mutate(pct.dif = pct.1 - pct.2) %>%
        arrange(-pct.dif) %>%
        filter(pct.dif > .1,
               pct.2 >= .5,
               avg_log2FC >= 1.5) %>%
        arrange(-avg_log2FC) %>%
        .$gene
)

markers_wtree_final_pvn <- union(
    all_markers_pvn_wtree_final %>%
        dplyr::filter(cluster %in% c('24'),
                      p_val <= 1e-4) %>%
        mutate(pct.dif = pct.1 - pct.2) %>%
        arrange(-pct.dif) %>%
        filter(pct.dif > .1,
               pct.1 >= .1,
               avg_log2FC >= .15,
               gene %in% big_genes) %>%
        arrange(-avg_log2FC) %>%
        .$gene,
    all_markers_pvn_wtree_final %>%
        dplyr::filter(cluster %in% c('24'),
                      p_val_adj <= 1e-3) %>%
        mutate(pct.dif = pct.1 - pct.2) %>%
        arrange(-pct.dif) %>%
        filter(pct.dif > .1,
               pct.1 >= .1,
               avg_log2FC >= 1) %>%
        arrange(-avg_log2FC) %>%
        .$gene
)

mcrs_age_corr_genes_intersect_3 <-
    unique(
        c(
            markers_wtree_final_trh,
            intersect(x = mcrs_age_corr_genes %>%
                          slice_max(order_by = Trh,
                                    n = 25),
                      y = mcrs_age_corr_genes %>%
                          slice_max(order_by = Mc4r,
                                    n = 25))$term,
            mcrs_age_corr_genes %>%
                slice_max(order_by = Trh, n = 3) %>%
                filter(Mc4r >= 0.1) %>%
                .$term,
            mcrs_age_corr_genes %>%
                slice_max(order_by = Mc4r, prop = 0.01) %>%
                filter(Trh >= 0.3) %>%
                .$term,
            markers_wtree_final_pvn,
            "Mc1r", "Mc2r", "Mc3r", "Mc4r",
            "Opcml", "Trh", "Nell2", "Brinp3",
            "Irs4", "Lrfn5", "Gpr101", "Mbnl2",
            "Nr3c2", "Pura", "Mbnl3",
            "Rora", "Fos", "Ghr", "Ebf3",
            "Cartpt", "Nrxn3", "Sstr2", "Slc1a4",
            "Alk", "Scgn", "Crh", "Pgf", "Ptn",
            "Mdk", "Alkal2", "Fam150b", "Tmem18",
            "Ptprz1", "Ptprf", "Ptprd", "Ptprt",
            "Ptpra", "Ptprb", "Oxt",
            "Slitrk1", "Slitrk2", "Slitrk4", "Slitrk5", "Slitrk6",
            pvn_genes,
            dmh_trh_g
        )
    )

srt_mcrs_filtered <-
    GetAssayData(rar2020.srt.pub.pne, 'data', 'SCT') %>%
    as.data.frame() %>%
    .[unique(c(mcrs_age_corr_genes_intersect_3)) %>%
          .[. %in% row.names(rar2020.srt.pub.pne@assays$SCT@data)], ] %>%
    as.data.frame() %>%
    t()

custom_genes <- c(corrr::correlate(srt_mcrs_filtered) %>% arrange(Mc4r) %>% select(term, Mc4r) %>% .$term)
custom_genes <-
  c("Alk", "Crh",
    custom_genes[!custom_genes %in% c("Trh", "Mc4r", "Alk", "Crh",
                                      "Ptprz1", "Ptprf", "Ptprd",
                                      "Ptprt", "Ptpra", "Avp", "Oxt",
                                      "Slitrk4","Slitrk5")],
    "Ptprz1", "Ptprf", "Ptprd", "Ptprt", "Ptpra",
    "Avp", "Oxt", "Slitrk4", "Slitrk5",
    "Trh", "Mc4r")

custom_order <-
    srt_mcrs_filtered %>%
    as.data.frame() %>%
    rownames_to_column(var = 'cname') %>%
    select(cname, Mc4r, Opcml, Trh, Nell2,
           Brinp3, Irs4, Lrfn5, Gpr101,
           Mbnl2, Nr3c2, Pura, Mbnl3,
           Ptn, Rora, Fos, Ghr, Ebf3,
           Cartpt, Nrxn3, Sstr2,
           Alk, Scgn,
           Crh, Pgf) %>%
    arrange(-Mc4r, -Opcml, -Trh, -Nell2,
            -Brinp3, -Irs4, -Lrfn5, -Gpr101,
            -Mbnl2, -Nr3c2, -Pura, -Mbnl3,
            Ptn, Rora, Fos, Ghr, Ebf3,
            Cartpt, Nrxn3, Sstr2,
            Alk, Scgn,
            Crh, Pgf) %>%
    .$cname

p <- DoHeatmap(rar2020.srt.pub.pne, cells = custom_order, features = custom_genes, size = 3, assay = 'SCT', slot = 'scale.data')
sPlot(name = "24-pne-adult",
      plt = p,
      type = "heatmap",
      asp = 7/6)

## gene pairs correlations in 24
# sct corrected only
p_corrs <- list(
    FeatureScatter(rar2020.srt.pub.pne, feature1 = "Alk", feature2 = "Mc4r", slot = "data"),
    FeatureScatter(rar2020.srt.pub.pne, feature1 = "Crh",  feature2 = "Alk", slot = "data"),
    FeatureScatter(rar2020.srt.pub.pne, feature1 = "Trh",  feature2 = "Alk", slot = "data"),
    FeatureScatter(rar2020.srt.pub.pne, feature1 = "Crh", feature2 = "Mc4r", slot = "data"),
    FeatureScatter(rar2020.srt.pub.pne, feature1 = "Trh", feature2 = "Mc4r", slot = "data"),
    FeatureScatter(rar2020.srt.pub.pne, feature1 = "Crh",  feature2 = "Trh", slot = "data")
)
n_corrs <- list(
    "24-def-sct-data-Alk-Mc4r",
    "24-def-sct-data-Crh-Alk" ,
    "24-def-sct-data-Trh-Alk" ,
    "24-def-sct-data-Crh-Mc4r",
    "24-def-sct-data-Trh-Mc4r",
    "24-def-sct-data-Crh-Trh"
)

purrr::walk2(n_corrs, p_corrs, sPlot, type = "corr-plt")

srt_mcrs_filtered <-
    GetAssayData(rar2020.srt.pub.pne, 'data', 'SCT') %>%
    as.data.frame() %>%
    .[unique(c(mcrs_age_corr_genes_intersect_3)) %>%
          .[. %in% row.names(rar2020.srt.pub.pne@assays$SCT@data)], ] %>%
    as.data.frame() %>%
    t()

minim <-
    srt_mcrs_filtered %>%
    as_tibble() %>%
    dplyr::summarise(across(.fns = function(x) quantile(x, .1)))

plot_data <-
    (srt_mcrs_filtered > as.double(minim)) %>%
    as_tibble() %>%
    mutate_all(as.numeric)

pdf(file=here::here(plots_dir, stringr::str_glue("upset", as.character("24-pne-adult_simple"), ".pdf", .sep = "_")), onefile=FALSE, height = 10, width = 1.618 * 10)
upset(as.data.frame(plot_data), order.by = 'freq',
      sets.x.label = 'Number of cells',
      text.scale = c(2, 1.6, 2, 1.3, 2, 3),
      nsets = 10,
      sets = c('Mc4r', 'Trh', 'Crh', 'Alk'),
      nintersects = 60,
      empty.intersections = "on")
dev.off()

pdf(file=here::here(plots_dir, stringr::str_glue("upset", as.character("24-pne-adult_test-compare"), ".pdf", .sep = "_")), onefile=FALSE, height = 10, width = 1.618 * 10)
upset(as.data.frame(plot_data), order.by = 'freq',
      sets.x.label = 'Number of cells',
      text.scale = c(2, 1.6, 2, 1.3, 2, 3),
      nsets = 10,
      sets = c('Trh', 'Crh', 'Mc4r',
               'Ptn', 'Pura', 'Rora',
               'Scgn', 'Ebf3', 'Nr3c2',
               'Irs4', 'Gpr101', 'Sstr2',
               'Ghr', 'Opcml', 'Alk', 'Slc1a4'),
      nintersects = 60)
dev.off()

pdf(file=here::here(plots_dir, stringr::str_glue("upset", as.character("24-pne-adult_wt-query"), ".pdf", .sep = "_")), onefile=FALSE, height = 10, width = 1.618 * 10)
upset(as.data.frame(plot_data), order.by = 'freq',
      sets.x.label = 'Number of cells',
      text.scale = c(2, 1.6, 2, 1.3, 2, 3),
      nsets = 10,
      sets = c('Mc4r', 'Trh', 'Crh', 'Alk'),
      queries = list(
          list(
              query = intersects,
              params = list('Mc4r', 'Trh'),
              active = T),
          list(
              query = intersects,
              params = list('Crh', 'Alk'))),
      nintersects = 60,
      empty.intersections = "on")
dev.off()


plot_data_extra <-
    plot_data %>%
    as_tibble() %>%
    dplyr::transmute("Trh+" = Trh, "Alk+" = Alk, "Mc4r+" = Mc4r, "Crh+" = Crh) %>%
    dplyr::bind_cols(as_tibble(srt_mcrs_filtered))

library(ggstatsplot)
scatter_corr <- function(mydata, x, y) {
    plot <- (ggstatsplot::ggscatterstats(
        data = mydata, # dataframe from which variables are taken
        x = x, # predictor/independent variable
        y = y, # dependent variable
        xfill = "#CC79A7", # fill for marginals on the x-axis
        yfill = "#009E73", # fill for marginals on the y-axis
    ))
}

matrix_corr <- function(mydata) {
    mydata <- dplyr::select(mydata, Trh, Alk, Mc4r, Crh)
    plot <- (ggstatsplot::ggcorrmat(
        data = mydata,
        ggcorrplot.args = list(outline.color = "black", hc.order = TRUE)
    ))
}

library(dabestr)
dabest_plot <- function(mydata, x, y) {
    two.group.unpaired <-
        mydata %>%
        mutate(groups = plyr::mapvalues(`Mc4r+`, from = c(0, 1), to = c("negative", "positive"))) %>%
        dabestr::dabest(x = groups, y = y,
                        # The idx below passes "Control" as the control group,
                        # and "Group1" as the test group. The mean difference
                        # will be computed as mean(Group1) - mean(Control1).
                        idx = c("negative", "positive"),
                        paired = FALSE)
    myplot <- (two.group.unpaired %>%
                   hedges_g() %>%
                   plot())
}


pdf(file=here::here(plots_dir, stringr::str_glue("upset", as.character("24-pne-adult_complex-wt-query"), ".pdf", .sep = "_")), onefile=FALSE, height = 14, width = 1.218 * 12)
upset(as.data.frame(plot_data_extra),
      sets.x.label = 'Number of cells',
      text.scale = c(2, 1.6, 2, 1.3, 2, 3),
      nsets = 4,
      sets = c('Trh+', 'Crh+', 'Mc4r+', 'Alk+'),
      nintersects = 20,
      main.bar.color = "black",
      mb.ratio = c(0.5, 0.5),
      queries = list(
          list(
              query = intersects,
              params = list('Mc4r+', 'Trh+'),
              active = T),
          list(
              query = intersects,
              params = list('Crh+', 'Alk+'))),
      attribute.plots = list(gridrows = 100,
                             plots = list(
                                 # list(plot = scatter_corr,
                                 list(plot = scatter_plot,
                                      x = "Trh",
                                      y = "Mc4r",
                                      queries = T),
                                 # list(plot = scatter_corr,
                                 list(plot = scatter_plot,
                                      x = "Crh",
                                      y = "Alk",
                                      queries = T)
                                 # list(plot = dabest_plot,
                                 #      x = "Mc4r+",
                                 #      y = "Alk",
                                 #      queries = T),
                                 # list(plot = matrix_corr,
                                 #      queries = F)
                             ),
                             ncols = 2))
dev.off()

#### gene pairs correlations ggstatsplot in 24 ####
#### sct corrected only ####
p_corrs <- list(
    ggstatsplot::ggscatterstats(as.data.frame(plot_data_extra), x = Alk,  y = Mc4r, xfill = "yellow", yfill = "magenta"),
    ggstatsplot::ggscatterstats(as.data.frame(plot_data_extra), x = Crh,  y = Alk,  xfill = "blue",   yfill = "yellow"),
    ggstatsplot::ggscatterstats(as.data.frame(plot_data_extra), x = Trh,  y = Alk,  xfill = "cyan",   yfill = "yellow"),
    ggstatsplot::ggscatterstats(as.data.frame(plot_data_extra), x = Crh,  y = Mc4r, xfill = "blue",   yfill = "magenta"),
    ggstatsplot::ggscatterstats(as.data.frame(plot_data_extra), x = Trh,  y = Mc4r, xfill = "cyan",   yfill = "magenta"),
    ggstatsplot::ggscatterstats(as.data.frame(plot_data_extra), x = Crh,  y = Trh,  xfill = "blue",   yfill = "cyan"),
    ggstatsplot::ggscatterstats(as.data.frame(plot_data_extra), y = Crh,  x = Alk,  yfill = "blue",   xfill = "yellow"),
    ggstatsplot::ggscatterstats(as.data.frame(plot_data_extra), y = Trh,  x = Alk,  yfill = "cyan",   xfill = "yellow"),
    ggstatsplot::ggscatterstats(as.data.frame(plot_data_extra), y = Sst,  x = Alk,  yfill = "pink",   xfill = "yellow"),
    ggstatsplot::ggscatterstats(as.data.frame(plot_data_extra), y = Oxt,  x = Alk,  yfill = "red",    xfill = "yellow")
)
n_corrs <- list(
    "24-sct-data-Alk-Mc4r",
    "24-sct-data-Crh-Alk" ,
    "24-sct-data-Trh-Alk" ,
    "24-sct-data-Crh-Mc4r",
    "24-sct-data-Trh-Mc4r",
    "24-sct-data-Crh-Trh",
    "24-sct-data-Alk-Crh" ,
    "24-sct-data-Alk-Trh" ,
    "24-sct-data-Alk-Sst" ,
    "24-sct-data-Alk-Oxt"
)

purrr::walk2(n_corrs, p_corrs, sPlot, type = "stat-corr-plt")
pSD <- matrix_corr(plot_data_extra)
sPlot(name = "24-sct-data",
      plt = pSD,
      type = "stat-corr-mtx",
      asp = 1.218)

## gene pairs correlations ggstatsplot in 24
#### sct corrected and scaled ####
srt_mcrs_filtered_scaled <-
    GetAssayData(rar2020.srt.pub.pne, 'scale.data', 'SCT') %>%
    as.data.frame() %>%
    .[unique(c(mcrs_age_corr_genes_intersect_3)) %>%
          .[. %in% row.names(rar2020.srt.pub.pne@assays$SCT@scale.data)], ] %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame()

p_corrs <- list(
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_scaled, x = Alk,  y = Mc4r, xfill = "yellow", yfill = "magenta"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_scaled, x = Crh,  y = Alk,  xfill = "blue",   yfill = "yellow"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_scaled, x = Trh,  y = Alk,  xfill = "cyan",   yfill = "yellow"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_scaled, x = Crh,  y = Mc4r, xfill = "blue",   yfill = "magenta"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_scaled, x = Trh,  y = Mc4r, xfill = "cyan",   yfill = "magenta"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_scaled, x = Crh,  y = Trh,  xfill = "blue",   yfill = "cyan"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_scaled, y = Crh,  x = Alk,  yfill = "blue",   xfill = "yellow"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_scaled, y = Trh,  x = Alk,  yfill = "cyan",   xfill = "yellow"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_scaled, y = Sst,  x = Alk,  yfill = "pink",   xfill = "yellow"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_scaled, y = Oxt,  x = Alk,  yfill = "red",    xfill = "yellow")
)
n_corrs <- list(
    "24-sct-scaledata-Alk-Mc4r",
    "24-sct-scaledata-Crh-Alk" ,
    "24-sct-scaledata-Trh-Alk" ,
    "24-sct-scaledata-Crh-Mc4r",
    "24-sct-scaledata-Trh-Mc4r",
    "24-sct-scaledata-Crh-Trh",
    "24-sct-scaledata-Alk-Crh" ,
    "24-sct-scaledata-Alk-Trh" ,
    "24-sct-scaledata-Alk-Sst" ,
    "24-sct-scaledata-Alk-Oxt"
)
purrr::walk2(n_corrs, p_corrs, sPlot, type = "stat-corr-plt")

pSS <- matrix_corr(srt_mcrs_filtered_scaled)
sPlot(name = "24-sct-scaledata",
      plt = pSS,
      type = "stat-corr-mtx",
      asp = 1.218)

## gene pairs correlations ggstatsplot in 24
#### rna corrected only ####
srt_mcrs_filtered_rna <-
    GetAssayData(rar2020.srt.pub.pne, 'data', 'RNA') %>%
    as.data.frame() %>%
    .[unique(c(mcrs_age_corr_genes_intersect_3)) %>%
          .[. %in% row.names(rar2020.srt.pub.pne@assays$RNA@data)], ] %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame()

p_corrs <- list(
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_rna, x = Alk,  y = Mc4r, xfill = "yellow", yfill = "magenta"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_rna, x = Crh,  y = Alk,  xfill = "blue",   yfill = "yellow"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_rna, x = Trh,  y = Alk,  xfill = "cyan",   yfill = "yellow"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_rna, x = Crh,  y = Mc4r, xfill = "blue",   yfill = "magenta"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_rna, x = Trh,  y = Mc4r, xfill = "cyan",   yfill = "magenta"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_rna, x = Crh,  y = Trh,  xfill = "blue",   yfill = "cyan"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_rna, y = Crh,  x = Alk,  yfill = "blue",   xfill = "yellow"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_rna, y = Trh,  x = Alk,  yfill = "cyan",   xfill = "yellow"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_rna, y = Sst,  x = Alk,  yfill = "pink",   xfill = "yellow"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_rna, y = Oxt,  x = Alk,  yfill = "red",    xfill = "yellow")
)
n_corrs <- list(
    "24-rna-data-Alk-Mc4r",
    "24-rna-data-Crh-Alk" ,
    "24-rna-data-Trh-Alk" ,
    "24-rna-data-Crh-Mc4r",
    "24-rna-data-Trh-Mc4r",
    "24-rna-data-Crh-Trh",
    "24-rna-data-Alk-Crh" ,
    "24-rna-data-Alk-Trh" ,
    "24-rna-data-Alk-Sst" ,
    "24-rna-data-Alk-Oxt"
)
purrr::walk2(n_corrs, p_corrs, sPlot, type = "stat-corr-plt")

pRD <- matrix_corr(srt_mcrs_filtered_rna)
sPlot(name = "24-rna-data",
      plt = pRD,
      type = "stat-corr-mtx",
      asp = 1.218)


## gene pairs correlations ggstatsplot in 24
#### rna corrected and scaled ####
srt_mcrs_filtered_scaled_rna <-
    GetAssayData(rar2020.srt.pub.pne, 'scale.data', 'RNA') %>%
    as.data.frame() %>%
    .[unique(c(mcrs_age_corr_genes_intersect_3)) %>%
          .[. %in% row.names(rar2020.srt.pub.pne@assays$RNA@scale.data)], ] %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame()

p_corrs <- list(
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_scaled_rna, x = Alk,  y = Mc4r, xfill = "yellow", yfill = "magenta"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_scaled_rna, x = Crh,  y = Alk,  xfill = "blue",   yfill = "yellow"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_scaled_rna, x = Trh,  y = Alk,  xfill = "cyan",   yfill = "yellow"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_scaled_rna, x = Crh,  y = Mc4r, xfill = "blue",   yfill = "magenta"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_scaled_rna, x = Trh,  y = Mc4r, xfill = "cyan",   yfill = "magenta"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_scaled_rna, x = Crh,  y = Trh,  xfill = "blue",   yfill = "cyan"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_scaled_rna, y = Crh,  x = Alk,  yfill = "blue",   xfill = "yellow"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_scaled_rna, y = Trh,  x = Alk,  yfill = "cyan",   xfill = "yellow"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_scaled_rna, y = Sst,  x = Alk,  yfill = "pink",   xfill = "yellow"),
    ggstatsplot::ggscatterstats(srt_mcrs_filtered_scaled_rna, y = Oxt,  x = Alk,  yfill = "red",    xfill = "yellow")
)
n_corrs <- list(
    "24-rna-scaledata-Alk-Mc4r",
    "24-rna-scaledata-Crh-Alk" ,
    "24-rna-scaledata-Trh-Alk" ,
    "24-rna-scaledata-Crh-Mc4r",
    "24-rna-scaledata-Trh-Mc4r",
    "24-rna-scaledata-Crh-Trh",
    "24-rna-scaledata-Alk-Crh" ,
    "24-rna-scaledata-Alk-Trh" ,
    "24-rna-scaledata-Alk-Sst" ,
    "24-rna-scaledata-Alk-Oxt"
)
purrr::walk2(n_corrs, p_corrs, sPlot, type = "stat-corr-plt")

pRS <- matrix_corr(srt_mcrs_filtered_scaled_rna)
sPlot(name = "24-rna-scaledata",
      plt = pRS,
      type = "stat-corr-mtx",
      asp = 1.218)

#### DABEst ####
pdf(file=here::here(plots_dir, stringr::str_glue("gardner-altman", as.character("24-pne-adult_Trh"), ".pdf", .sep = "_")), onefile=FALSE, height = 10, width = 1.618 * 10)
dabest_plot(as.data.frame(plot_data_extra), x="Mc4r", y = "Trh")
dev.off()

pdf(file=here::here(plots_dir, stringr::str_glue("gardner-altman", as.character("24-pne-adult_Crh"), ".pdf", .sep = "_")), onefile=FALSE, height = 10, width = 1.618 * 10)
dabest_plot(as.data.frame(plot_data_extra), x="Mc4r", y = "Crh")
dev.off()

pdf(file=here::here(plots_dir, stringr::str_glue("gardner-altman", as.character("24-pne-adult_Alk"), ".pdf", .sep = "_")), onefile=FALSE, height = 10, width = 1.618 * 10)
dabest_plot(as.data.frame(plot_data_extra), x="Mc4r", y = "Alk")
dev.off()


