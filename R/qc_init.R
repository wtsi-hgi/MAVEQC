#' load the configure file
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
