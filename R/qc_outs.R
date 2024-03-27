#' initialize function
setGeneric("qcout_samqc_all", function(object, ...) {
  standardGeneric("qcout_samqc_all")
})

#' create all the output files
#'
#' @export
#' @name qcout_samqc_all
#' @param object   sampleQC object
#' @param qc_type  qc type
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_all",
    signature = "sampleQC",
    definition = function(object,
                          qc_type = c("plasmid", "screen"),
                          out_dir = NULL) {
        if (is.null(out_dir)) {
            stop(paste0("====> Error: out_dir is not provided, no output directory."))
        }

        qc_type <- match.arg(qc_type)

        qcout_samqc_meta(object = object, out_dir = out_dir)
        qcout_samqc_cutoffs(object = object, out_dir = out_dir)
        qcout_samqc_readlens(object = object, out_dir = out_dir)
        qcout_samqc_total(object = object, out_dir = out_dir)
        qcout_samqc_missing(object = object, out_dir = out_dir)
        qcout_samqc_accepted(object = object, out_dir = out_dir)
        qcout_samqc_libcov(object = object, out_dir = out_dir)
        qcout_samqc_pos_cov(object = object, qc_type = qc_type, out_dir = out_dir)
        qcout_samqc_results(object = object, qc_type = qc_type, out_dir = out_dir)

        if (qc_type == "screen") {
            qcout_samqc_pos_anno(object = object, out_dir = out_dir)
        }

        qcout_samqc_badseqs(object = object, qc_type = qc_type, out_dir = out_dir)
    }
)

#' initialize function
setGeneric("qcout_samqc_meta", function(object, ...) {
  standardGeneric("qcout_samqc_meta")
})

#' create output file of bad seqs which fail filtering
#'
#' @export
#' @name qcout_samqc_meta
#' @param object   sampleQC object
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_meta",
    signature = "sampleQC",
    definition = function(object,
                          out_dir = NULL) {
        if (is.null(out_dir)) {
            stop(paste0("====> Error: out_dir is not provided, no output directory."))
        }

        write.table(object@samples_meta,
                    file = paste0(out_dir, "/", "sample_qc_meta.tsv"),
                    quote = FALSE,
                    sep = "\t",
                    row.names = FALSE,
                    col.names = TRUE)
    }
)

#' initialize function
setGeneric("qcout_samqc_cutoffs", function(object, ...) {
  standardGeneric("qcout_samqc_cutoffs")
})

#' create output file of bad seqs which fail filtering
#'
#' @export
#' @name qcout_samqc_cutoffs
#' @param object   sampleQC object
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_cutoffs",
    signature = "sampleQC",
    definition = function(object,
                          out_dir = NULL) {
        if (is.null(out_dir)) {
            stop(paste0("====> Error: out_dir is not provided, no output directory."))
        }

        write.table(object@cutoffs,
                    file = paste0(out_dir, "/", "sample_qc_cutoffs.tsv"),
                    quote = FALSE,
                    sep = "\t",
                    row.names = FALSE,
                    col.names = TRUE)
    }
)

#' initialize function
setGeneric("qcout_samqc_badseqs", function(object, ...) {
  standardGeneric("qcout_samqc_badseqs")
})

#' create output file of bad seqs which fail filtering for screen qc
#'
#' @export
#' @name qcout_samqc_badseqs
#' @param object   sampleQC object
#' @param qc_type  plasmid or screen
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_badseqs",
    signature = "sampleQC",
    definition = function(object,
                          qc_type = c("plasmid", "screen"),
                          out_dir = NULL) {
        if (is.null(out_dir)) {
            stop(paste0("====> Error: out_dir is not provided, no output directory."))
        }

        qc_type <- match.arg(qc_type)

        if (qc_type == "screen") {
            cat("Outputing bad sequences filtered out by clustering...", "\n", sep = "")
            bad_seqs <- data.table()
            for (i in 1:length(object@bad_seqs_bycluster)) {
                bad_tmp <- object@bad_seqs_bycluster[[i]]
                colnames(bad_tmp)[1] <- "seq"
                lib_tmp <- object@samples[[i]]@vep_anno
                res_tmp <- bad_tmp[lib_tmp, on = .(seq), nomatch = 0][, c("unique_oligo_name", "seq", "count")]
                colnames(res_tmp)[3] <- names(object@bad_seqs_bycluster)[i]

                if (nrow(bad_seqs) == 0) {
                    bad_seqs <- res_tmp
                } else {
                    bad_seqs <- merge(bad_seqs, res_tmp, by = c("unique_oligo_name", "seq"), all = TRUE)
                }
            }

            write.table(bad_seqs,
                        file = paste0(out_dir, "/", "failed_variants_by_cluster.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)

            cat("Outputing bad sequences filtered out by sequencing depth...", "\n", sep = "")
            bad_seqs <- data.table()
            for (i in 1:length(object@bad_seqs_bydepth)) {
                bad_tmp <- object@bad_seqs_bydepth[[i]]
                colnames(bad_tmp)[1] <- "seq"
                lib_tmp <- object@samples[[i]]@vep_anno
                res_tmp <- bad_tmp[lib_tmp, on = .(seq), nomatch = 0][, c("unique_oligo_name", "seq", "count")]
                colnames(res_tmp)[3] <- names(object@bad_seqs_bydepth)[i]

                if (nrow(bad_seqs) == 0) {
                    bad_seqs <- res_tmp
                } else {
                    bad_seqs <- merge(bad_seqs, res_tmp, by = c("unique_oligo_name", "seq"), all = TRUE)
                }
            }

            write.table(bad_seqs,
                        file = paste0(out_dir, "/", "failed_variants_by_depth.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)

            cat("Outputing bad sequences filtered out by library mapping...", "\n", sep = "")
            bad_seqs <- data.table()
            for (i in 1:length(object@bad_seqs_bylib)) {
                bad_tmp <- object@bad_seqs_bylib[[i]]
                colnames(bad_tmp)[1] <- "seq"
                lib_tmp <- object@samples[[i]]@vep_anno
                res_tmp <- bad_tmp[lib_tmp, on = .(seq), nomatch = 0][, c("unique_oligo_name", "seq", "count")]
                colnames(res_tmp)[3] <- names(object@bad_seqs_bylib)[i]

                if (nrow(bad_seqs) == 0) {
                    bad_seqs <- res_tmp
                } else {
                    bad_seqs <- merge(bad_seqs, res_tmp, by = c("unique_oligo_name", "seq"), all = TRUE)
                }
            }

            write.table(bad_seqs,
                        file = paste0(out_dir, "/", "failed_variants_by_mapping.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)
        }

        cat("Outputing missing variants in the library...", "\n", sep = "")
        bad_seqs <- data.table()
        for (i in 1:length(object@samples)) {
            bad_tmp <- object@samples[[i]]@libcounts[count == 0]

            if (nrow(bad_seqs) == 0) {
                bad_seqs <- bad_tmp
            } else {
                bad_seqs <- merge(bad_seqs, bad_tmp, all = TRUE)
            }
        }

        write.table(bad_seqs,
                    file = paste0(out_dir, "/", "missing_variants_in_library.tsv"),
                    quote = FALSE,
                    sep = "\t",
                    row.names = FALSE,
                    col.names = TRUE)
    }
)

#' initialize function
setGeneric("qcout_samqc_readlens", function(object, ...) {
  standardGeneric("qcout_samqc_readlens")
})

#' create output file of total reads stats
#'
#' @export
#' @name qcout_samqc_readlens
#' @param object    sampleQC object
#' @param len_bins  the bins of length distribution
#' @param out_dir   the output directory
setMethod(
    "qcout_samqc_readlens",
    signature = "sampleQC",
    definition = function(object,
                          len_bins = seq(0, 300, 50),
                          out_dir = NULL) {
        cols <- c("Group",
                  "Sample",
                  "Sample Info",
                  "Sample Exon",
                  "Total Reads",
                  "% 0 ~ 50",
                  "% 50 ~ 100",
                  "% 100 ~ 150",
                  "% 150 ~ 200",
                  "% 200 ~ 250",
                  "% 250 ~ 300",
                  "Pass Threshold (%)",
                  "Pass")
        df_outs <- data.frame(matrix(NA, nrow(object@stats), length(cols)))
        colnames(df_outs) <- cols

        df_outs[, 1] <- object@samples[[1]]@libname
        df_outs[, 2] <- rownames(object@stats)
        df_outs[, 3] <- object@samples_meta$sample_info
        df_outs[, 4] <- paste(object@samples_meta$transcript_id, object@samples_meta$exon_num, sep = ":")
        df_outs[, 5] <- object@stats$total_reads

        bin_per <- data.frame()
        for (i in 1:length(object@lengths)) {
            tmp_lens <- object@lengths[[i]]$length
            h <- hist(tmp_lens, breaks = len_bins, plot = FALSE)
            bin_per <- rbind(bin_per, round(h$counts / nrow(object@samples[[i]]@allcounts) * 100, 1))
        }

        df_outs[, 6] <- bin_per[, 1]
        df_outs[, 7] <- bin_per[, 2]
        df_outs[, 8] <- bin_per[, 3]
        df_outs[, 9] <- bin_per[, 4]
        df_outs[, 10] <- bin_per[, 5]
        df_outs[, 11] <- bin_per[, 6]
        df_outs[, 12] <- 90
        df_outs[, 13] <- (df_outs[, 10] + df_outs[, 11]) > df_outs[, 12]

        df_outs <- df_outs[match(mixedsort(df_outs$Sample), df_outs$Sample), ]

        if (length(out_dir) == 0) {
            reactable(df_outs, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,
                      theme = reactableTheme(
                          style = list(fontFamily = "-apple-system", fontSize = "0.75rem")),
                      columns = list(
                          "Group" = colDef(minWidth = 100),
                          "Sample" = colDef(minWidth = 100),
                          "Total Reads" = colDef(format = colFormat(separators = TRUE)),
                          "Pass" = colDef(cell = function(value) {
                                                   if (value) "\u2705" else "\u274c" })),
                      rowStyle = function(index) { if (!(df_outs[index, "Pass"])) { list(background = t_col("tomato", 0.2)) } }
                     )
        } else {
            write.table(df_outs,
                        file = paste0(out_dir, "/", "sample_qc_read_length.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)
        }
    }
)

#' initialize function
setGeneric("qcout_samqc_missing", function(object, ...) {
  standardGeneric("qcout_samqc_missing")
})

#' create output file of bad seqs which fail filtering
#'
#' @export
#' @name qcout_samqc_missing
#' @param object   sampleQC object
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_missing",
    signature = "sampleQC",
    definition = function(object,
                          out_dir = NULL) {
        cols <- c("Group",
                  "Sample",
                  "Sample Info",
                  "Sample Exon",
                  "Library Sequences",
                  "Missing Library Sequences",
                  "% Missing Library Sequences",
                  "Pass Threshold (%)",
                  "Pass")
        df_outs <- data.frame(matrix(NA, nrow(object@stats), length(cols)))
        colnames(df_outs) <- cols

        df_outs[, 1] <- object@samples[[1]]@libname
        df_outs[, 2] <- rownames(object@stats)
        df_outs[, 3] <- object@samples_meta$sample_info
        df_outs[, 4] <- paste(object@samples_meta$transcript_id, object@samples_meta$exon_num, sep = ":")
        df_outs[, 5] <- object@stats$library_seqs
        df_outs[, 6] <- object@stats$missing_meta_seqs
        tmp_out <- object@stats$per_missing_meta_seqs * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 7] <- tmp_out
        tmp_out <- object@cutoffs$per_missing_variants * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 8] <- tmp_out
        df_outs[, 9] <- object@stats$qcpass_missing_per

        df_outs <- df_outs[match(mixedsort(df_outs$Sample), df_outs$Sample), ]

        if (length(out_dir) == 0) {
            reactable(df_outs, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,
                      theme = reactableTheme(
                          style = list(fontFamily = "-apple-system", fontSize = "0.75rem")),
                      columns = list(
                          "Group" = colDef(minWidth = 100),
                          "Sample" = colDef(minWidth = 100),
                          "Library Sequences" = colDef(format = colFormat(separators = TRUE)),
                          "Missing Library Sequences" = colDef(format = colFormat(separators = TRUE)),
                          "% Missing Library Sequences" = colDef(format = colFormat(separators = TRUE),
                                                                 style = function(value) {
                                                                             if (value > object@cutoffs$per_missing_variants * 100) {
                                                                                 color <- "red"
                                                                                 fweight <- "bold"
                                                                             } else {
                                                                                 color <- "forestgreen"
                                                                                 fweight <- "plain"
                                                                             }
                                                                             list(color = color, fontWeight = fweight)}),
                          "Pass" = colDef(cell = function(value) {
                                                   if (value) "\u2705" else "\u274c" })),
                      rowStyle = function(index) { if (!(df_outs[index, "Pass"])) { list(background = t_col("tomato", 0.2)) } }
                     )
        } else {
            write.table(df_outs,
                        file = paste0(out_dir, "/", "sample_qc_stats_missing.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)
        }
    }
)

#' initialize function
setGeneric("qcout_samqc_total", function(object, ...) {
  standardGeneric("qcout_samqc_total")
})

#' create output file of total reads stats
#'
#' @export
#' @name qcout_samqc_total
#' @param object   sampleQC object
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_total",
    signature = "sampleQC",
    definition = function(object,
                          out_dir = NULL) {
        cols <- c("Group",
                  "Sample",
                  "Sample Info",
                  "Sample Exon",
                  "Accepted Reads",
                  "% Accepted Reads",
                  "Excluded Reads",
                  "% Excluded Reads",
                  "Total Reads",
                  "Pass Threshold",
                  "Pass")
        df_outs <- data.frame(matrix(NA, nrow(object@stats), length(cols)))
        colnames(df_outs) <- cols

        df_outs[, 1] <- object@samples[[1]]@libname
        df_outs[, 2] <- rownames(object@stats)
        df_outs[, 3] <- object@samples_meta$sample_info
        df_outs[, 4] <- paste(object@samples_meta$transcript_id, object@samples_meta$exon_num, sep = ":")
        df_outs[, 5] <- object@stats$accepted_reads
        tmp_out <- object@stats$accepted_reads / object@stats$total_reads * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 6] <- tmp_out
        df_outs[, 7] <- object@stats$excluded_reads
        tmp_out <- object@stats$excluded_reads / object@stats$total_reads * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 8] <- tmp_out
        df_outs[, 9] <- object@stats$total_reads
        df_outs[, 10] <- object@cutoffs$num_total_reads
        df_outs[, 11] <- object@stats$qcpass_accepted_reads

        df_outs <- df_outs[match(mixedsort(df_outs$Sample), df_outs$Sample), ]

        if (length(out_dir) == 0) {
            reactable(df_outs, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,
                      theme = reactableTheme(
                          style = list(fontFamily = "-apple-system", fontSize = "0.75rem")),
                      columns = list(
                          "Group" = colDef(minWidth = 100),
                          "Sample" = colDef(minWidth = 100),
                          "Accepted Reads" = colDef(format = colFormat(separators = TRUE)),
                          "Excluded Reads" = colDef(format = colFormat(separators = TRUE)),
                          "Total Reads" = colDef(format = colFormat(separators = TRUE),
                                                 style = function(value) {
                                                             if (value < object@cutoffs$num_total_reads) {
                                                                 color <- "red"
                                                                 fweight <- "bold"
                                                             } else {
                                                                 color <- "forestgreen"
                                                                 fweight <- "plain"
                                                             }
                                                             list(color = color, fontWeight = fweight)}),
                          "Pass Threshold" = colDef(format = colFormat(separators = TRUE)),
                          "Pass" = colDef(cell = function(value) {
                                                   if (value) "\u2705" else "\u274c" })),
                      rowStyle = function(index) { if (!(df_outs[index, "Pass"])) { list(background = t_col("tomato", 0.2)) } }
                     )
        } else {
            write.table(df_outs,
                        file = paste0(out_dir, "/", "sample_qc_stats_total.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)
        }
    }
)

#' initialize function
setGeneric("qcout_samqc_accepted", function(object, ...) {
  standardGeneric("qcout_samqc_accepted")
})

#' create output file of library reads stats
#'
#' @export
#' @name qcout_samqc_accepted
#' @param object   sampleQC object
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_accepted",
    signature = "sampleQC",
    definition = function(object,
                          out_dir = NULL) {
        cols <- c("Group",
                  "Sample",
                  "Sample Info",
                  "Sample Exon",
                  "% Library Reads",
                  "% Reference Reads",
                  "% PAM Reads",
                  "% Unmapped Reads",
                  "Pass Threshold (%)",
                  "Pass")
        df_outs <- data.frame(matrix(NA, nrow(object@stats), length(cols)))
        colnames(df_outs) <- cols

        df_outs[, 1] <- object@samples[[1]]@libname
        df_outs[, 2] <- rownames(object@stats)
        df_outs[, 3] <- object@samples_meta$sample_info
        df_outs[, 4] <- paste(object@samples_meta$transcript_id, object@samples_meta$exon_num, sep = ":")
        tmp_out <- object@stats$per_library_reads * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 5] <- tmp_out
        tmp_out <- object@stats$per_ref_reads * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 6] <- tmp_out
        tmp_out <- object@stats$per_pam_reads * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 7] <- tmp_out
        tmp_out <- object@stats$per_unmapped_reads * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 8] <- tmp_out
        df_outs[, 9] <- object@cutoffs$per_library_reads * 100
        df_outs[, 10] <- object@stats$qcpass_library_per

        df_outs <- df_outs[match(mixedsort(df_outs$Sample), df_outs$Sample), ]

        if (length(out_dir) == 0) {
            reactable(df_outs, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,
                      theme = reactableTheme(
                          style = list(fontFamily = "-apple-system", fontSize = "0.75rem")),
                      columns = list(
                          "Group" = colDef(minWidth = 100),
                          "Sample" = colDef(minWidth = 100),
                          "% Library Reads" = colDef(style = function(value) {
                                                                 if (value < object@cutoffs$per_library_reads * 100) {
                                                                    color <- "red"
                                                                    fweight <- "bold"
                                                                } else {
                                                                    color <- "forestgreen"
                                                                    fweight <- "plain"
                                                                }
                                                                list(color = color, fontWeight = fweight)}),
                          "Pass" = colDef(cell = function(value) {
                                                   if (value) "\u2705" else "\u274c" })),
                      rowStyle = function(index) { if (!(df_outs[index, "Pass"])) { list(background = t_col("tomato", 0.2)) } }
                     )
        } else {
            write.table(df_outs,
                        file = paste0(out_dir, "/", "sample_qc_stats_accepted.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)
        }
    }
)

#' initialize function
setGeneric("qcout_samqc_libcov", function(object, ...) {
  standardGeneric("qcout_samqc_libcov")
})

#' create output file of library coverage
#'
#' @export
#' @name qcout_samqc_libcov
#' @param object   sampleQC object
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_libcov",
    signature = "sampleQC",
    definition = function(object,
                          out_dir = NULL) {
        cols <- c("Group",
                  "Sample",
                  "Sample Info",
                  "Sample Exon",
                  "Total Library Reads",
                  "Total Library Sequences",
                  "Library Coverage",
                  "Median Coverage",
                  "Pass Threshold",
                  "Pass")
        df_outs <- data.frame(matrix(NA, nrow(object@stats), length(cols)))
        colnames(df_outs) <- cols

        df_outs[, 1] <- object@samples[[1]]@libname
        df_outs[, 2] <- rownames(object@stats)
        df_outs[, 3] <- object@samples_meta$sample_info
        df_outs[, 4] <- paste(object@samples_meta$transcript_id, object@samples_meta$exon_num, sep = ":")
        df_outs[, 5] <- object@stats$library_reads
        df_outs[, 6] <- object@stats$library_seqs
        tmp_out <- object@stats$library_cov
        tmp_out <- sapply(tmp_out, function(x) round(x, 0))
        df_outs[, 7] <- tmp_out
        df_outs[, 8] <- object@stats$median_cov
        df_outs[, 9] <- object@cutoffs$library_cov
        df_outs[, 10] <- object@stats$qcpass_library_cov

        df_outs <- df_outs[match(mixedsort(df_outs$Sample), df_outs$Sample), ]

        if (length(out_dir) == 0) {
            reactable(df_outs, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,
                      theme = reactableTheme(
                          style = list(fontFamily = "-apple-system", fontSize = "0.75rem")),
                      columns = list(
                          "Group" = colDef(minWidth = 100),
                          "Sample" = colDef(minWidth = 100),
                          "Total Library Reads" = colDef(format = colFormat(separators = TRUE)),
                          "Total Library Sequences" = colDef(format = colFormat(separators = TRUE)),
                          "Library Coverage" = colDef(format = colFormat(separators = TRUE),
                                                      style = function(value) {
                                                                  if (value < object@cutoffs$library_cov) {
                                                                      color <- "red"
                                                                      fweight <- "bold"
                                                                  } else {
                                                                      color <- "forestgreen"
                                                                      fweight <- "plain"
                                                                  }
                                                                  list(color = color, fontWeight = fweight)}),
                          "Pass" = colDef(cell = function(value) {
                                                   if (value) "\u2705" else "\u274c" })),
                      rowStyle = function(index) { if (!(df_outs[index, "Pass"])) { list(background = t_col("tomato", 0.2)) } }
                     )
        } else {
            write.table(df_outs,
                        file = paste0(out_dir, "/", "sample_qc_stats_coverage.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)
        }
    }
)

#' initialize function
setGeneric("qcout_samqc_pos_cov", function(object, ...) {
  standardGeneric("qcout_samqc_pos_cov")
})

#' create output file of lof percentages
#'
#' @export
#' @name qcout_samqc_pos_cov
#' @param object   sampleQC object
#' @param qc_type  screen or plasmid
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_pos_cov",
    signature = "sampleQC",
    definition = function(object,
                          qc_type = c("plasmid", "screen"),
                          out_dir = NULL) {
        qc_type <- match.arg(qc_type)

        library_dependent_counts <- data.table()
        if (qc_type == "screen") {
            library_dependent_counts <- merge_list_to_dt(object@library_counts, "sequence", "count")
            library_dependent_counts[, sequence := NULL]
        } else {
            list_counts <- list()
            for (i in 1:length(object@samples)) {
                sample <- object@samples[[i]]@sample
                list_counts[[sample]] <- object@samples[[i]]@libcounts[, c("id", "count")]
            }

            library_dependent_counts <- merge_list_to_dt(list_counts, "id", "count")
            library_dependent_counts[, id := NULL]
        }

        write.table(library_dependent_counts,
                    file = paste0(out_dir, "/", "sample_qc_stats_pos_counts.tsv"),
                    quote = FALSE,
                    sep = "\t",
                    row.names = FALSE,
                    col.names = TRUE)

        cols <- c("Group",
                  "Sample",
                  "Sample Info",
                  "Sample Exon",
                  "Chromosome",
                  "Strand",
                  "Genomic Start",
                  "Genomic End",
                  "% Low Abundance",
                  "Low Abundance cutoff",
                  "Pass Threshold (%)",
                  "Pass")
        df_outs <- data.frame(matrix(NA, nrow(object@stats), length(cols)))
        colnames(df_outs) <- cols

        df_outs[, 1] <- object@samples[[1]]@libname
        df_outs[, 2] <- rownames(object@stats)
        df_outs[, 3] <- object@samples_meta$sample_info
        df_outs[, 4] <- paste(object@samples_meta$transcript_id, object@samples_meta$exon_num, sep = ":")
        df_outs[, 5] <- sapply(object@library_counts_chr, function (x) x[[1]])
        df_outs[, 6] <- sapply(object@library_counts_chr, function (x) x[[2]])
        df_outs[, 7] <- sapply(object@library_counts_chr, function (x) x[[3]])
        df_outs[, 8] <- sapply(object@library_counts_chr, function (x) x[[4]])

        low_per <- vector()
        if (qc_type == "screen") {
            for (s in object@samples) {
                tmp_num <- sum(object@library_counts_pos[[s@sample]]$count < object@cutoffs$seq_low_count)
                low_per <- append(low_per, round(tmp_num / nrow(object@library_counts_pos[[s@sample]]) * 100, 2))
            }
        } else {
            for (s in object@samples) {
                tmp_num <- sum(s@libcounts$count < object@cutoffs$seq_low_count)
                low_per <- append(low_per, round(tmp_num / nrow(object@library_counts_pos[[s@sample]]) * 100, 2))
            }
        }
        df_outs[, 9] <- low_per

        df_outs[, 10] <- object@cutoffs$seq_low_count
        df_outs[, 11] <- (1 - object@cutoffs$low_abundance_lib_per) * 100
        df_outs[, 12] <- df_outs[, 9] < (1 - object@cutoffs$low_abundance_lib_per) * 100

        df_outs <- df_outs[match(mixedsort(df_outs$Sample), df_outs$Sample), ]

        if (length(out_dir) == 0) {
            reactable(df_outs, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,
                      theme = reactableTheme(
                          style = list(fontFamily = "-apple-system", fontSize = "0.75rem")),
                      columns = list(
                          "Group" = colDef(minWidth = 100),
                          "Sample" = colDef(minWidth = 100),
                          "Genomic Start" = colDef(format = colFormat(separators = TRUE)),
                          "Genomic End" = colDef(format = colFormat(separators = TRUE)),
                          "% Low Abundance" = colDef(minWidth = 200,
                                                           style = function(value) {
                                                                       if (value > (1 - object@cutoffs$low_abundance_lib_per) * 100) {
                                                                           color <- "red"
                                                                           fweight <- "bold"
                                                                       } else {
                                                                           color <- "forestgreen"
                                                                           fweight <- "plain"
                                                                       }
                                                                       list(color = color, fontWeight = fweight)}),
                          "Low Abundance cutoff" = colDef(minWidth = 200),
                          "Pass" = colDef(cell = function(value) {
                                                   if (value) "\u2705" else "\u274c" })),
                      rowStyle = function(index) { if (!(df_outs[index, "Pass"])) { list(background = t_col("tomato", 0.2)) } }
                     )
        } else {
            write.table(df_outs,
                        file = paste0(out_dir, "/", "sample_qc_stats_pos_coverage.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)
        }
    }
)

#' initialize function
setGeneric("qcout_samqc_pos_anno", function(object, ...) {
  standardGeneric("qcout_samqc_pos_anno")
})

#' create output file of lof percentages
#'
#' @export
#' @name qcout_samqc_pos_anno
#' @param object   sampleQC object
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_pos_anno",
    signature = "sampleQC",
    definition = function(object,
                          out_dir = NULL) {
        cols <- c("Group",
                  "Sample",
                  "Sample Info",
                  "Sample Exon",
                  "Chromosome",
                  "Strand",
                  "Genomic Start",
                  "Genomic End",
                  "% Low Abundance (LOF)",
                  "% Low Abundance (Others)",
                  "% Low Abundance (ALL)",
                  "% Low Abundance cutoff",
                  "Pass Threshold",
                  "Pass")
        df_outs <- data.frame(matrix(NA, nrow(object@stats), length(cols)))
        colnames(df_outs) <- cols

        df_outs[, 1] <- object@samples[[1]]@libname
        df_outs[, 2] <- rownames(object@stats)
        df_outs[, 3] <- object@samples_meta$sample_info
        df_outs[, 4] <- paste(object@samples_meta$transcript_id, object@samples_meta$exon_num, sep = ":")
        df_outs[, 5] <- sapply(object@library_counts_chr, function (x) x[[1]])
        df_outs[, 6] <- sapply(object@library_counts_chr, function (x) x[[2]])
        df_outs[, 7] <- sapply(object@library_counts_chr, function (x) x[[3]])
        df_outs[, 8] <- sapply(object@library_counts_chr, function (x) x[[4]])

        libcounts_pos <- as.data.frame(object@library_counts_pos_anno)
        libcounts_pos <- libcounts_pos[, c(rownames(object@stats), "consequence")]
        libcounts_pos$consequence <- ifelse(libcounts_pos$consequence == "LOF", "LOF", "Others")
        libcounts_pos[, rownames(object@stats)] <- t(t(libcounts_pos[, rownames(object@stats)]) / object@stats$accepted_reads * 100)

        # what about NA?
        #libcounts_pos[is.na(libcounts_pos)] <- 0

        lof_counts <- libcounts_pos[libcounts_pos$consequence == "LOF", rownames(object@stats)]
        # the number of seqs with low abundance
        lof_low_num <- colSums(lof_counts < object@cutoffs$low_abundance_per * 100, na.rm = TRUE)
        # the percentage of seqs with low abundance
        lof_low_per <- lof_low_num / nrow(libcounts_pos) * 100
        lof_low_per <- round(lof_low_per, 1)

        others_counts <- libcounts_pos[libcounts_pos$consequence == "Others", rownames(object@stats)]
        # the number of seqs with low abundance
        others_low_num <- colSums(others_counts < object@cutoffs$low_abundance_per * 100, na.rm = TRUE)
        # the percentage of seqs with low abundance
        others_low_per <- others_low_num / nrow(libcounts_pos) * 100
        others_low_per <- round(others_low_per, 1)

        df_outs[, 9] <- lof_low_per
        df_outs[, 10] <- others_low_per
        df_outs[, 11] <- lof_low_per + others_low_per
        df_outs[, 12] <- object@cutoffs$low_abundance_per * 100
        df_outs[, 13] <- (1 - object@cutoffs$low_abundance_lib_per) * 100
        df_outs[, 14] <- df_outs[, 11] < (1 - object@cutoffs$low_abundance_lib_per) * 100

        df_outs <- df_outs[match(mixedsort(df_outs$Sample), df_outs$Sample), ]

        if (length(out_dir) == 0) {
            reactable(df_outs, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,
                      theme = reactableTheme(
                          style = list(fontFamily = "-apple-system", fontSize = "0.75rem")),
                      columns = list(
                          "Group" = colDef(minWidth = 100),
                          "Sample" = colDef(minWidth = 100),
                          "Genomic Start" = colDef(format = colFormat(separators = TRUE)),
                          "Genomic End" = colDef(format = colFormat(separators = TRUE)),
                          "% Low Abundance (LOF)" = colDef(minWidth = 200),
                          "% Low Abundance (Others)" = colDef(minWidth = 200),
                          "% Low Abundance (ALL)" = colDef(minWidth = 200,
                                                           style = function(value) {
                                                                       if (value > (1 - object@cutoffs$low_abundance_lib_per) * 100) {
                                                                           color <- "red"
                                                                           fweight <- "bold"
                                                                       } else {
                                                                           color <- "forestgreen"
                                                                           fweight <- "plain"
                                                                       }
                                                                       list(color = color, fontWeight = fweight)}),
                          "% Low Abundance cutoff" = colDef(minWidth = 200),
                          "Pass" = colDef(cell = function(value) {
                                                   if (value) "\u2705" else "\u274c" })),
                      rowStyle = function(index) { if (!(df_outs[index, "Pass"])) { list(background = t_col("tomato", 0.2)) } }
                     )
        } else {
            write.table(df_outs,
                        file = paste0(out_dir, "/", "sample_qc_stats_pos_percentage.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)
        }
    }
)

#' initialize function
setGeneric("qcout_samqc_results", function(object, ...) {
  standardGeneric("qcout_samqc_results")
})

#' create all the output files
#'
#' @export
#' @name qcout_samqc_results
#' @param object   sampleQC object
#' @param qc_type  qc type
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_results",
    signature = "sampleQC",
    definition = function(object,
                          qc_type = c("plasmid", "screen"),
                          out_dir = NULL) {
        qc_type <- match.arg(qc_type)

        cols <- c("Group",
                  "Sample",
                  "Sample Info",
                  "Sample Exon",
                  "Gini Coefficient",
                  "Gini Coefficient Pass",
                  "Total Reads",
                  "% Missing Variants",
                  "Accepted Reads",
                  "% Mapping Reads",
                  "% Reference Reads",
                  "% Library Reads",
                  "Library Coverage",
                  "% R1 Adaptor",
                  "% R2 Adaptor")
        df_outs <- data.frame(matrix(NA, nrow(object@stats), length(cols)))
        colnames(df_outs) <- cols

        df_outs[, 1] <- object@samples[[1]]@libname
        df_outs[, 2] <- rownames(object@stats)
        df_outs[, 3] <- object@samples_meta$sample_info
        df_outs[, 4] <- paste(object@samples_meta$transcript_id, object@samples_meta$exon_num, sep = ":")
        df_outs[, 5] <- object@stats$gini_coeff_before_qc
        df_outs[, 6] <- ifelse(df_outs[, 3] < maveqc_config$gini_coeff, TRUE, FALSE)
        df_outs[, 7] <- object@stats$qcpass_total_reads
        df_outs[, 8] <- object@stats$qcpass_missing_per
        df_outs[, 9] <- object@stats$qcpass_accepted_reads
        df_outs[, 10] <- object@stats$qcpass_mapping_per
        df_outs[, 11] <- object@stats$qcpass_ref_per
        df_outs[, 12] <- object@stats$qcpass_library_per
        df_outs[, 13] <- object@stats$qcpass_library_cov
        df_outs[, 14] <- object@stats$per_r1_adaptor
        df_outs[, 15] <- object@stats$per_r2_adaptor

        df_outs <- df_outs[match(mixedsort(df_outs$Sample), df_outs$Sample), ]

        if (qc_type == "screen") {
            df_outs <- df_outs
        }

        if (length(out_dir) == 0) {
            reactable(df_outs, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,
                      theme = reactableTheme(
                          style = list(fontFamily = "-apple-system", fontSize = "0.75rem")),
                      columns = list(
                          "Group" = colDef(minWidth = 100),
                          "Sample" = colDef(minWidth = 100),
                          "Total Reads" = colDef(cell = function(value) { if (value) "\u2705" else "\u274c" }),
                          "% Missing Variants" = colDef(cell = function(value) { if (value) "\u2705" else "\u274c" }),
                          "Accepted Reads" = colDef(cell = function(value) { if (value) "\u2705" else "\u274c" }),
                          "% Mapping Reads" = colDef(cell = function(value) { if (value) "\u2705" else "\u274c" }),
                          "% Reference Reads" = colDef(cell = function(value) { if (value) "\u2705" else "\u274c" }),
                          "% Library Reads" = colDef(cell = function(value) { if (value) "\u2705" else "\u274c" }),
                          "Library Coverage" = colDef(cell = function(value) { if (value) "\u2705" else "\u274c" }))
                     )
        } else {
            write.table(df_outs,
                        file = paste0(out_dir, "/", "sample_qc_results.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)
        }
    }
)


#####################################################################################################################################################

#' initialize function
setGeneric("qcout_expqc_all", function(object, ...) {
  standardGeneric("qcout_expqc_all")
})

#' create all the output files
#'
#' @export
#' @name qcout_expqc_all
#' @param object   experimentQC object
#' @param out_dir  the output directory
setMethod(
    "qcout_expqc_all",
    signature = "experimentQC",
    definition = function(object,
                          out_dir = NULL) {
        if (is.null(out_dir)) {
            stop(paste0("====> Error: out_dir is not provided, no output directory."))
        }

        qcout_expqc_corr(object = object, out_dir = out_dir)
        qcout_expqc_deseq(object = object, eqc_type = "all", out_dir = out_dir)
    }
)

#' initialize function
setGeneric("qcout_expqc_corr", function(object, ...) {
  standardGeneric("qcout_expqc_corr")
})

#' create output file of clustering and correlation results
#'
#' @export
#' @name qcout_expqc_corr
#' @param object   experimentQC object
#' @param out_dir  the output directory
setMethod(
    "qcout_expqc_corr",
    signature = "experimentQC",
    definition = function(object,
                          out_dir = NULL) {
        df_outs <- as.data.frame(object@lib_corr_res)

        num_clusters <- length(unique(object@coldata$condition))
        name_clusters <- as.vector(unique(object@coldata$condition))
        sample_clusters <- cutree(object@lib_hclust_res, num_clusters)

        df_outs <- cbind(object@coldata[rownames(df_outs), ], df_outs)
        colnames(df_outs)[1] <- "Replicate"
        colnames(df_outs)[2] <- "Condition"
        df_outs <- cbind(sample_clusters[rownames(df_outs)], df_outs)
        colnames(df_outs)[1] <- "Cluster"
        df_outs <- cbind(rownames(df_outs), df_outs)
        colnames(df_outs)[1] <- "Sample"

        df_outs$Pass <- NA
        for (i in 1:length(name_clusters)) {
            tmp_clusters <- df_outs[df_outs$Condition == name_clusters[i], ]$Cluster
            if (length(unique(tmp_clusters)) == 1) {
                df_outs[df_outs$Condition == name_clusters[i], ]$Pass <- TRUE
            } else {
                df_outs[df_outs$Condition == name_clusters[i], ]$Pass <- FALSE
            }
        }

        df_outs_tmp <- unique(df_outs[df_outs$Pass, c("Condition", "Cluster")])
        for (i in 1:nrow(df_outs_tmp)) {
            pass_check <- TRUE
            for (j in 1:nrow(df_outs_tmp)) {
                if (i != j) {
                    if (df_outs_tmp[i, ]$Cluster == df_outs_tmp[j, ]$Cluster) {
                        pass_check <- FALSE
                    }
                }
            }
            df_outs[df_outs$Condition == df_outs_tmp[i, ]$Condition, ]$Pass <- pass_check
        }

        df_outs <- df_outs[match(mixedsort(df_outs$Sample), df_outs$Sample), ]

        is.num <- sapply(df_outs, is.numeric)
        df_outs[is.num] <- lapply(df_outs[is.num], round, 2)

        if (length(out_dir) == 0) {
            reactable(df_outs, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,
                      theme = reactableTheme(
                          style = list(fontFamily = "-apple-system", fontSize = "0.75rem")),
                      columns = list(
                          "Sample" = colDef(minWidth = 100),
                          "Cluster" = colDef(minWidth = 100),
                          "Replicate" = colDef(minWidth = 100),
                          "Condition" = colDef(minWidth = 100),
                          "Pass" = colDef(cell = function(value) {
                                                   if (value) "\u2705" else "\u274c" })),
                      rowStyle = function(index) { if (!(df_outs[index, "Pass"])) { list(background = t_col("tomato", 0.2)) } }
                     )
        } else {
            write.table(df_outs,
                        file = paste0(out_dir, "/", "experiment_qc_corr.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)
        }
    }
)

#' initialize function
setGeneric("qcout_expqc_deseq", function(object, ...) {
  standardGeneric("qcout_expqc_deseq")
})

#' create all the output files
#'
#' @export
#' @name qcout_expqc_deseq
#' @param object   experimentQC object
#' @param eqc_type library counts or all counts
#' @param out_dir  the output directory
setMethod(
    "qcout_expqc_deseq",
    signature = "experimentQC",
    definition = function(object,
                          eqc_type = c("lib", "all"),
                          out_dir = NULL) {
        if (is.null(out_dir)) {
            stop(paste0("====> Error: out_dir is not provided, no output directory."))
        }

        eqc_type <- match.arg(eqc_type)

        for (i in 1:length(object@comparisons)) {
            if (eqc_type == "lib") {
                df_outs <- object@lib_deseq_res_anno[[i]]
            } else {
                df_outs <- object@all_deseq_res_anno_adj[[i]]
            }

            tmpcons <- sort(unique(df_outs$consequence))
            cols <- c("consequence",
                      "number of depleted",
                      "number of no impact",
                      "number of enriched",
                      "lfc range",
                      "total number",
                      "total number of outside range")
            df_outs_summary <- data.frame(matrix(NA, length(tmpcons), length(cols)))
            colnames(df_outs_summary) <- cols

            df_outs_summary[, 1] <- tmpcons
            df_outs_summary[, 5] <- paste0(maveqc_config$expqc_lfc_min, " ~ ", maveqc_config$expqc_lfc_max)
            for (j in 1:nrow(df_outs_summary)) {
                df_tmp <- df_outs[consequence == df_outs_summary$consequence[j]]

                df_outs_summary[j, 2] <- nrow(df_tmp[stat == "depleted"])
                df_outs_summary[j, 3] <- nrow(df_tmp[stat == "no impact"])
                df_outs_summary[j, 4] <- nrow(df_tmp[stat == "enriched"])

                df_outs_summary[j, 6] <- nrow(df_tmp)
                df_outs_summary[j, 7] <- nrow(df_tmp[adj_log2FoldChange < maveqc_config$expqc_lfc_min | adj_log2FoldChange > maveqc_config$expqc_lfc_max])
            }

            dt_outs_summary <- as.data.table(df_outs_summary)
            dt1 <- setorder(dt_outs_summary[consequence %in% c("Synonymous_Variant", "LOF")], -consequence)
            dt2 <- dt_outs_summary[consequence %nin% c("Synonymous_Variant", "LOF")]
            dt_outs_summary <- rbind(dt1, dt2)

            write.table(df_outs,
                        file = paste0(out_dir, "/", "experiment_qc_deseq_fc.", object@comparisons[[i]], ".", eqc_type, ".tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)

            write.table(dt_outs_summary,
                        file = paste0(out_dir, "/", "experiment_qc_deseq_fc.", object@comparisons[[i]], ".", eqc_type, "_sum.tsv"),
                        quote = FALSE,
                        sep = "\t",
                        row.names = FALSE,
                        col.names = TRUE)
        }
    }
)