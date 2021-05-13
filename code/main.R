tabMetaDatSrt <- function(dat.srt,
                          agg.dir = 2,
                          var1 = "age",
                          var2 = "wtree",
                          ...) {
    dat.srt@meta.data %>%
        dplyr::select(var1, var2) %>%
        table() %>%
        prop.table(margin = agg.dir) %>%
        .[, colSums(is.na(.)) != nrow(.)] %>% # this line deletes NaN columns
        `*`(100) %>%
        knitr::kable(format = "html",
                     digits = 2,
                     caption = sprintf("Proportion table (%s) of %s by %s for %s\n",
                                       ifelse(agg.dir == 1,
                                              "by rows",
                                              "by columns"),
                                       var1,
                                       var2,
                                       dat.srt@project.name))
}

tabMetaDatSrtPipe <- function(dat.srt,
                              agg.dir = 2,
                              var1 = "age",
                              var2 = "wtree",
                              ...) {
    dat.srt@meta.data %>%
        dplyr::select(var1, var2) %>%
        table() %>%
        prop.table(margin = agg.dir) %>%
        .[, colSums(is.na(.)) != nrow(.)] %>% # this line deletes NaN columns
        `*`(100) %>%
        knitr::kable(format = "pipe",
                     digits = 2,
                     caption = sprintf("Proportion table (%s) of %s by %s for %s\n",
                                       ifelse(agg.dir == 1,
                                              "by rows",
                                              "by columns"),
                                       var1,
                                       var2,
                                       dat.srt@project.name))
}

subsetSrt <- function(dat.srt, srt.name, checkmd = TRUE, ...) {
    init.num.cells <- dim(dat.srt)[2]
    dat.srt <- subset(dat.srt, ...)
    dat.srt@project.name <-
        sprintf("%s subset from %s", srt.name, dat.srt@project.name)

    if (checkmd) {
        tabMetaDatSrt(dat.srt, ...)
    }

    message(sprintf("number of cells for %s is %s;\n it is %s percent of initial %s\n",
                    dat.srt@project.name,
                    dim(dat.srt)[2],
                    round(dim(dat.srt)[2]/init.num.cells * 100, digits = 2),
                    init.num.cells))
    return(dat.srt)
}

## plot gene within cell population
pltUpSet <- function(gene = "Alk", pop = 'pneSS', pdata = plot_data3) {
    UpSetR::upset(as.data.frame(pdata),
          order.by = c("freq", "degree"),
          sets.x.label = 'Number of cells',
          sets = c(gene, pop),
          keep.order = T,
          empty.intersections = F,
          text.scale = c(2, 1.6, 2, 1.3, 2, 3),
          nsets = 6)
}

#' Available RAM in kB
CheckRAM <- function() as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern = TRUE))


#' Available cores
AvailableCores <- function(prop2use = .9) {
    max(1, floor(parallel::detectCores(logical = FALSE)*prop2use))
}


#' Get number of cores to fit RAM needs
Cores4RAM <- function(need) max(1, min(AvailableCores(), floor(CheckRAM() / need)))

readRDS <- function(name, file) {
    assign(envir = .GlobalEnv, x = name, value = readr::read_rds(here::here(tables_dir, file)))
}
readTable <- function(name, file) {
    assign(envir = .GlobalEnv, x = name, value = readr::read_csv(here::here(tables_dir, file)))
}

# try({
#     rar2020.srt.pub <-
#         readr::read_rds(here::here(tables_dir,
#                                    '2021-05-13_001-086-446-605-669_mm10-hypothalamus_srt4-obj-umap.rds'))
#
#     colours <- readr::read_lines(here(data_dir, "colours_wtree.tsv"))
#     clrlev <- levels(rar2020.srt.pub$wtree)
#     names(colours) <- clrlev
#     pvn.cols <- colours[c(15, 24, 26, 31, 43)]
#     names(pvn.cols) %<>%
#         plyr::mapvalues(x = .,
#                         from = c(15, 24, 26, 31, 43),
#                         to = c("pneTRH",
#                                "pneCRH",
#                                "mneVAS",
#                                "pneSS",
#                                "mneOXY"))
#
#     rar2020.srt.pub <-
#         readr::read_rds(here::here(tables_dir,
#                                    '2021-05-13_001-086-446-605-669_adult_mm10-hypothalamus_srt4-obj-umap.rds'))
#     rar2020.srt.pvn.bck <-
#         readr::read_rds(here::here(tables_dir,
#                                    '2021-05-13_001-086-446-605-669_adult_mm10-hypothalamus-pvn_alk-thinness-paper-subset-srt4-obj-umap.rds'))
#     rar2020.srt.pvn <-
#         readr::read_rds(here::here(tables_dir,
#                                    '2021-05-13_001-086-446-605-669_adult_mm10-hypothalamus-pvn_bottom-lh-excl-srt4-obj-sct-umap.rds'))
#     all_markers_pvn_wtree_final <-
#         readr::read_csv(here::here(tables_dir,
#                                    '2021-05-13_001-086-446-605-669_adult_mm10-hypothalamus-pvn_deg-mrk-gns-ctypes-by-sct-wilcox.csv'))
#     srt_mcrs <-
#         readr::read_csv(here::here(tables_dir,
#                                    '2021-05-13_001-086-446-605-669_adult_mm10-hypothalamus-pvn_expr-mtx-cells-by-mrk-gns-int.csv')) %>%
#         column_to_rownames(var = "cell_names")
#     plot_data <-
#         readr::read_csv(here::here(tables_dir,
#                                    '2021-05-13_001-086-446-605-669_adult_mm10-hypothalamus-pvn_ncells-long-by-mrk-gns-int.csv'))
#     plot_data_filt <-
#         readr::read_csv(here::here(tables_dir,
#                                    '2021-05-13_001-086-446-605-669_adult_mm10-hypothalamus-pvn-mc4r-pos_ncells-long-by-mrk-gns-int.csv'))
#     plot_data3 <-
#         readr::read_csv(here::here(tables_dir,
#                                    '2021-05-13_001-086-446-605-669_adult_mm10-hypothalamus-pvn_ncells-long-by-mrk-gns-int-ctypes.csv'))
#
# })

