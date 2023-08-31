#' show basic info of a SGE object
#'
#' @export
#' @name show
#' @param object SGE object
setMethod(
    "show",
    signature = "SGE",
    definition = function(object) {
        cat("An object of class ", class(object), "\n", sep = "")
        cat("|--> sample name: ", object@sample, "\n", sep = "")
        cat("|--> library type: ", object@libtype, "\n", sep = "")
        cat("|--> library name: ", object@libname, "\n", sep = "")
        cat("    |--> 5' adaptor: ", object@adapt5, "\n", sep = "")
        cat("    |--> 3' adaptor: ", object@adapt3, "\n", sep = "")
        cat("    |--> ref seq: ", object@refseq, "\n", sep = "")
        cat("    |--> pam seq: ", object@pamseq, "\n", sep = "")
        cat("    |--> No. of library-dependent sequences: ", nrow(object@libcounts), "\n", sep = "")
        cat("    |--> No. of library-independent sequences: ", nrow(object@allcounts), "\n", sep = "")
        cat("|--> valiant meta: ", nrow(object@valiant_meta), " records and ", ncol(object@valiant_meta), " fields", "\n", sep = "")
        cat("    |--> ", sum(object@libcounts$name %in% object@valiant_meta$oligo_name), " library-dependent sequence IDs matched in valiant meta oligo names", "\n", sep = "")
    }
)

#' initialize function
setGeneric("show_stats", function(object, ...) {
  standardGeneric("show_stats")
})

#' show basic stats of a SGE object
#'
#' @export
#' @name show_stats
#' @param object SGE object
setMethod(
    "show_stats",
    signature = "SGE",
    definition = function(object) {
        colstrs <- colnames(object@libstats)
        dash_line <- paste0("|", strrep("-", 20 + nchar("type")), "-|-")
        dash_line <- paste0(dash_line, strrep("-", nchar("library dependent counts")), "-|-")
        dash_line <- paste0(dash_line, strrep("-", nchar("library independent counts")), "-|")
        cat("Basic stats of sample: ", object@sample, "\n", sep = "")
        cat(dash_line, "\n", sep = "")
        header_line <- paste0("|", strrep(" ", 20), "type | library dependent counts | library independent counts |")
        cat(header_line, "\n", sep = "")
        cat(dash_line, "\n", sep = "")
        for (i in 1:length(colstrs)) {
            spacestr <- strrep(" ", 20 - nchar(colstrs[i]) + nchar("type"))
            info_line <- paste0("|", spacestr, colstrs[i], " | ")
            spacestr <- strrep(" ", nchar("library dependent counts") - nchar(toString(object@libstats[, i])))
            info_line <- paste0(info_line, spacestr, object@libstats[, i], " | ")
            spacestr <- strrep(" ", nchar("library independent counts") - nchar(toString(object@allstats[, i])))
            info_line <- paste0(info_line, spacestr, object@allstats[, i], " | ")
            cat(info_line, "\n", sep = "")
        }
        cat(dash_line, "\n", sep = "")
    }
)

#' initialize function
setGeneric("show_qc_stats", function(object, ...) {
  standardGeneric("show_qc_stats")
})

#' show qc stats of a SGE object
#'
#' @export
#' @name show_qc_stats
#' @param object SGE object
setMethod(
    "show_qc_stats",
    signature = "SGE",
    definition = function(object) {
        colstrs <- colnames(object@libstats_qc)
        dash_line <- paste0("|", strrep("-", 20 + nchar("type")), "-|-")
        dash_line <- paste0(dash_line, strrep("-", nchar("library dependent counts")), "-|-")
        dash_line <- paste0(dash_line, strrep("-", nchar("library independent counts")), "-|")
        cat("QC stats of sample: ", object@sample, "\n", sep = "")
        cat(dash_line, "\n", sep = "")
        header_line <- paste0("|", strrep(" ", 20), "type | library dependent counts | library independent counts |")
        cat(header_line, "\n", sep = "")
        cat(dash_line, "\n", sep = "")
        for (i in 1:length(colstrs)) {
            spacestr <- strrep(" ", 20 - nchar(colstrs[i]) + nchar("type"))
            info_line <- paste0("|", spacestr, colstrs[i], " | ")
            spacestr <- strrep(" ", nchar("library dependent counts") - nchar(toString(object@libstats_qc[, i])))
            info_line <- paste0(info_line, spacestr, object@libstats_qc[, i], " | ")
            spacestr <- strrep(" ", nchar("library independent counts") - nchar(toString(object@allstats_qc[, i])))
            info_line <- paste0(info_line, spacestr, object@allstats_qc[, i], " | ")
            cat(info_line, "\n", sep = "")
        }
        cat(dash_line, "\n", sep = "")
    }
)

#' show basic info of a SGE object
#'
#' @export
#' @name show
#' @param object sampleQC object
setMethod(
    "show",
    signature = "sampleQC",
    definition = function(object) {
        cat("An object of class ", class(object), "\n", sep = "")
        cat("|--> samples: ", "\n", sep = "")
        for (s in object@samples) {
            cat("    |--> ", s@sample, "\n", sep = "")
        }

        cat("|--> reference samples: ", "\n", sep = "")
        if (length(object@samples_ref) == 0) {
            cat("    |--> no sample found", "\n", sep = "")
        } else {
            for (s in object@samples_ref) {
                cat("    |--> ", s@sample, "\n", sep = "")
            }
        }

        cat("|--> QC results: ", "\n", sep = "")
        for (i in 1:nrow(object@stats)) {
            cat("    |--> ", rownames(object@stats)[i], ": ", object@stats[i, ]$qcpass, "\n", sep = "")
        }
    }
)
