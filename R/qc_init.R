#' load the default config file
#'
#' @name .onLoad
.onLoad <- function(libname, pkgname) {
    quietly <- getOption("quietly")
    options(quietly = TRUE)

    pkg_info <- paste0("MAVEQC", "-v", packageVersion("MAVEQC"))
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
#' @name create_config
#' @param config_dir  the directory of config file
create_config <- function(config_dir) {
    config_path <- paste0(config_dir, "/", "config.yaml")

    sink(config_path)

    cat("# user defined thresholds for QC", "\n", sep = "")
    cat("\n", sep = "")

    cat("#-----------#", "\n", sep = "")
    cat("# Sample QC #", "\n", sep = "")
    cat("#-----------#", "\n", sep = "")
    cat("\n", sep = "")

    cat("# gini coefficient must be lower than 0.5", "\n", sep = "")
    cat("gini_coeff: 0.5", "\n", sep = "")
    cat("\n", sep = "")

    cat("# the sample must have more than 1000000 total reads", "\n", sep = "")
    cat("sqc_total: 1000000", "\n", sep = "")
    cat("\n", sep = "")

    cat("# the missing varaints in the library must be less than 10%", "\n", sep = "")
    cat("sqc_missing: 0.01", "\n", sep = "")
    cat("\n", sep = "")

    cat("# the sequence must have at least 5 counts in at least 25% of the samples", "\n", sep = "")
    cat("sqc_low_count: 5", "\n", sep = "")
    cat("sqc_low_sample_per: 0.25", "\n", sep = "")
    cat("\n", sep = "")

    cat("# the sample must have more than 1000000 reads after the low count filtering", "\n", sep = "")
    cat("sqc_accepted: 1000000", "\n", sep = "")
    cat("\n", sep = "")

    cat("# the sample must have more than 60% of reads aligned to the library including reference and PAM reads", "\n", sep = "")
    cat("sqc_mapping_per: 0.6", "\n", sep = "")
    cat("\n", sep = "")

    cat("# the sample must have less than 10% of reads aligned to reference sequence", "\n", sep = "")
    cat("sqc_ref_per: 0.1", "\n", sep = "")
    cat("\n", sep = "")

    cat("# the sample must have more than 40% of reads aligned to the library", "\n", sep = "")
    cat("sqc_library_per: 0.4", "\n", sep = "")
    cat("\n", sep = "")

    cat("# the sample must have more than 100x average coverage", "\n", sep = "")
    cat("# the number of library reads divided by the number of sequences", "\n", sep = "")
    cat("sqc_library_cov: 100", "\n", sep = "")
    cat("\n", sep = "")

    cat("# the majority of the variants (>70%) distributed above the 0.005% cutoff for the reference samples", "\n", sep = "")
    cat("sqc_low_per: 0.00005", "\n", sep = "")
    cat("sqc_low_lib_per: 0.7", "\n", sep = "")

    cat("#---------------#", "\n", sep = "")
    cat("# Experiment QC #", "\n", sep = "")
    cat("#---------------#", "\n", sep = "")
    cat("\n", sep = "")

    cat("# DESeq2 relevant cutoffs", "\n", sep = "")
    cat("expqc_padj: 0.05", "\n", sep = "")
    cat("expqc_lfc_depleted: 0", "\n", sep = "")
    cat("expqc_lfc_enriched: 0", "\n", sep = "")
    cat("expqc_top_variants: 500", "\n", sep = "")
    cat("\n", sep = "")

    cat("# Log2 Fold Change cutoffs", "\n", sep = "")
    cat("expqc_lfc_min: -6", "\n", sep = "")
    cat("expqc_lfc_max: 2", "\n", sep = "")

    sink()

}

#' load the user's config file
#'
#' @export
#' @name load_config
#' @param config_path  the path of config file
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
