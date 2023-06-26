#' initialize function
setGeneric("qcout_bad_seqs", function(object, ...) {
  standardGeneric("qcout_bad_seqs")
})

#' create output file of bad seqs which fail filtering
#'
#' @export
#' @param object  sampleQC object
#' @param outdir  the output directory
setMethod(
    "qcout_bad_seqs",
    signature = "sampleQC",
    definition = function(object,
                          outdir) {
        if (length(outdir) == 0) {
            stop(paste0("====> Error: outdir is not provided, no output directory."))
        }

        write.table(object@bad_seqs$by_cluster, 
                    file = paste0(outdir, "/", "failed_sequences_by_cluster.txt"),
                    quote = FALSE,
                    sep = "\t",
                    row.names = FALSE,
                    col.names = FALSE)
        write.table(object@bad_seqs$by_count, 
                    file = paste0(outdir, "/", "failed_sequences_by_count.txt"),
                    quote = FALSE,
                    sep = "\t",
                    row.names = FALSE,
                    col.names = FALSE)
        write.table(object@bad_seqs$by_effect, 
                    file = paste0(outdir, "/", "failed_sequences_by_mapping.txt"),
                    quote = FALSE,
                    sep = "\t",
                    row.names = FALSE,
                    col.names = FALSE)
    }
)

#' initialize function
setGeneric("qcout_sampleqc_stats", function(object, ...) {
  standardGeneric("qcout_sampleqc_stats")
})

#' create output file of bad seqs which fail filtering
#'
#' @export
#' @param object  sampleQC object
#' @param outdir  the output directory
setMethod(
    "qcout_sampleqc_stats",
    signature = "sampleQC",
    definition = function(object,
                          outdir) {
        if (length(outdir) == 0) {
            stop(paste0("====> Error: outdir is not provided, no output directory."))
        }

        write.table(object@stats, 
                    file = paste0(outdir, "/", "sampleqc_stats.tsv"),
                    quote = FALSE,
                    sep = "\t",
                    row.names = TRUE,
                    col.names = TRUE)
    }
)