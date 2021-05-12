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

#' Available RAM in kB
CheckRAM <- function() as.numeric(system("awk '/MemFree/ {print $2}' /proc/meminfo", intern = TRUE))


#' Available cores
AvailableCores <- function(prop2use = .9) {
    max(1, floor(parallel::detectCores(logical = FALSE)*prop2use))
}


#' Get number of cores to fit RAM needs
Cores4RAM <- function(need) max(1, min(AvailableCores(), floor(CheckRAM() / need)))

