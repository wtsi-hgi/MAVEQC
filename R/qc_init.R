#' load the default config file
#'
#' @name .onLoad
.onLoad <- function(libname, pkgname) {
    quietly <- getOption("quietly")
    options(quietly = TRUE)

    pkg_info <- "Welcome to MAVEQC"
    packageStartupMessage(pkg_info)
    options(quietly = quietly)

    config_path <- system.file("config.yaml", package = "MAVEQC")

    if (file.exists(config_path)) {
        maveqc_config <<- read.config(file = config_path)
    } else {
        warning("Configuration file \"config/config.yaml\" not found.")
    }
}

#' create the template of config file
#'
#' @export
#' @param config_dir  the directory of config file
#' @name create_config
create_config <- function(config_dir) {
    config_path <- paste0(config_dir, "/", "config.yaml")

    sink(report_path)

    cat("# user defined thresholds for QC", "\n", sep = "")
    cat("\n", sep = "")

    cat("# the sequence must have at least 5 counts in at least 25% of the samples", "\n", sep = "")
    cat("sqc_low_count: 5", "\n", sep = "")
    cat("sqc_low_sample_per: 0.25", "\n", sep = "")
    cat("\n", sep = "")

    cat("# the sample must have more than 1000000 reads after the low count filtering", "\n", sep = "")
    cat("sqc_filtered: 1000000", "\n", sep = "")
    cat("\n", sep = "")

    cat("# the sample must have more than 60% of reads aligned to the library including reference and PAM reads", "\n", sep = "")
    cat("sqc_mapping_per: 0.6", "\n", sep = "")
    cat("\n", sep = "")

    cat("# the sample must have more than 40% of library reads excluding reference and PAM reads", "\n", sep = "")
    cat("sqc_library_per: 0.4", "\n", sep = "")
    cat("\n", sep = "")

    cat("# the sample must have more than 100x average coverage", "\n", sep = "")
    cat("# the number of library reads divided by the number of sequences", "\n", sep = "")
    cat("sqc_library_cov: 100", "\n", sep = "")
    cat("\n", sep = "")

    cat("# the majority of the variants (>70%) distributed above the 0.005% cutoff for the reference samples", "\n", sep = "")
    cat("sqc_low_per: 0.00005", "\n", sep = "")
    cat("sqc_low_lib_per: 0.7", "\n", sep = "")

    sink()

}

#' load the user's config file
#'
#' @export
#' @param config_path  the path of config file
#' @name load_config
load_config <- function(config_path) {
    if (file.exists(config_path)) {
        user_config <- read.config(file = config_path)

        config_names <- names(maveqc_config)
        for (i in 1:length(config_names)) {
            if (!(is.null(user_config[[config_names[i]]]))) {
                maveqc_config[[config_names[i]]] <<- user_config[[config_names[i]]]
            }
        }
    } else {
        stop(paste0("====> Error: ", config_path, " is not found!"))
    }
}
