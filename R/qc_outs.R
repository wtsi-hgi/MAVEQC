#' initialize function
setGeneric("qcout_samqc_all", function(object, ...) {
  standardGeneric("qcout_samqc_all")
})

#' create all the output files
#'
#' @export
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

        if (qc_type == "plasmid") {
            qcout_samqc_cutoffs(object = object, out_dir = out_dir)
            qcout_samqc_readlens(object = object, out_dir = out_dir)
            qcout_samqc_total(object = object, out_dir = out_dir)
            qcout_samqc_accepted(object = object, out_dir = out_dir)
            qcout_samqc_libcov(object = object, out_dir = out_dir)
            qcout_samqc_pos_cov(object = object, qc_type = qc_type, out_dir = out_dir)
        } else {
            qcout_samqc_cutoffs(object = object, out_dir = out_dir)
            qcout_samqc_readlens(object = object, out_dir = out_dir)
            qcout_samqc_total(object = object, out_dir = out_dir)
            qcout_samqc_accepted(object = object, out_dir = out_dir)
            qcout_samqc_libcov(object = object, out_dir = out_dir)
            qcout_samqc_pos_cov(object = object, qc_type = qc_type, out_dir = out_dir)
            qcout_samqc_pos_anno(object = object, out_dir = out_dir)
        }
    }
)

#' initialize function
setGeneric("qcout_samqc_cutoffs", function(object, ...) {
  standardGeneric("qcout_samqc_cutoffs")
})

#' create output file of bad seqs which fail filtering
#'
#' @export
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

#' create output file of bad seqs which fail filtering
#'
#' @export
#' @param object   sampleQC object
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_badseqs",
    signature = "sampleQC",
    definition = function(object,
                          out_dir = NULL) {
        if (is.null(out_dir)) {
            stop(paste0("====> Error: out_dir is not provided, no output directory."))
        }

        cat("Outputing bad sequences filtered out by clustering...", "\n", sep = "")
        write.table(merge_list_to_df(object@bad_seqs_bycluster),
                    file = paste0(out_dir, "/", "failed_sequences_by_cluster.txt"),
                    quote = FALSE,
                    sep = "\t",
                    row.names = FALSE,
                    col.names = TRUE)

        cat("Outputing bad sequences filtered out by sequencing depth...", "\n", sep = "")
        write.table(object@bad_seqs_bydepth,
                    file = paste0(out_dir, "/", "failed_sequences_by_depth.txt"),
                    quote = FALSE,
                    sep = "\t",
                    row.names = FALSE,
                    col.names = TRUE)

        cat("Outputing bad sequences filtered out by library mapping...", "\n", sep = "")
        write.table(object@bad_seqs_bylib,
                    file = paste0(out_dir, "/", "failed_sequences_by_mapping.txt"),
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
                  "Total Reads",
                  "% 0 ~ 50",
                  "% 50 ~ 100",
                  "% 100 ~ 150",
                  "% 150 ~ 200",
                  "% 200 ~ 250",
                  "% 250 ~ 300",
                  "Pass Threshold",
                  "Pass")
        df_outs <- data.frame(matrix(NA, nrow(object@stats), length(cols)))
        colnames(df_outs) <- cols

        df_outs[, 1] <- object@samples[[1]]@libname
        df_outs[, 2] <- rownames(object@stats)
        df_outs[, 3] <- object@stats$total_reads

        bin_per <- data.frame()
        for (i in 1:length(object@lengths)) {
            tmp_lens <- object@lengths[[i]]$length
            h <- hist(tmp_lens, breaks = len_bins, plot = FALSE)
            bin_per <- rbind(bin_per, round(h$counts / nrow(object@samples[[i]]@allcounts) * 100, 1))
        }

        df_outs[, 4] <- bin_per[, 1]
        df_outs[, 5] <- bin_per[, 2]
        df_outs[, 6] <- bin_per[, 3]
        df_outs[, 7] <- bin_per[, 4]
        df_outs[, 8] <- bin_per[, 5]
        df_outs[, 9] <- bin_per[, 6]
        df_outs[, 10] <- 90
        df_outs[, 11] <- (df_outs[, 8] + df_outs[, 9]) > df_outs[, 10]

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
                                                   if (value) "\u2705" else "\u274c" }))
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
setGeneric("qcout_samqc_total", function(object, ...) {
  standardGeneric("qcout_samqc_total")
})

#' create output file of total reads stats
#'
#' @export
#' @param object   sampleQC object
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_total",
    signature = "sampleQC",
    definition = function(object,
                          out_dir = NULL) {
        cols <- c("Group",
                  "Sample",
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
        df_outs[, 3] <- object@stats$accepted_reads
        tmp_out <- object@stats$accepted_reads / object@stats$total_reads * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 4] <- tmp_out
        df_outs[, 5] <- object@stats$excluded_reads
        tmp_out <- object@stats$excluded_reads / object@stats$total_reads * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 6] <- tmp_out
        df_outs[, 7] <- object@stats$total_reads
        df_outs[, 8] <- object@cutoffs$total_reads
        df_outs[, 9] <- object@stats$qcpass_accepted_reads

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
                                                             if (value < object@cutoffs$total_reads) {
                                                                 color <- "red"
                                                                 fweight <- "bold"
                                                             } else {
                                                                 color <- "forestgreen"
                                                                 fweight <- "plain"
                                                             }
                                                             list(color = color, fontWeight = fweight)}),
                          "Pass Threshold" = colDef(format = colFormat(separators = TRUE)),
                          "Pass" = colDef(cell = function(value) {
                                                   if (value) "\u2705" else "\u274c" }))
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
#' @param object   sampleQC object
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_accepted",
    signature = "sampleQC",
    definition = function(object,
                          out_dir = NULL) {
        cols <- c("Group",
                  "Sample",
                  "% Library Reads",
                  "% Reference Reads",
                  "% PAM Reads",
                  "% Unmapped Reads",
                  "Pass Threshold",
                  "Pass")
        df_outs <- data.frame(matrix(NA, nrow(object@stats), length(cols)))
        colnames(df_outs) <- cols

        df_outs[, 1] <- object@samples[[1]]@libname
        df_outs[, 2] <- rownames(object@stats)
        tmp_out <- object@stats$per_library_reads * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 3] <- tmp_out
        tmp_out <- object@stats$per_ref_reads * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 4] <- tmp_out
        tmp_out <- object@stats$per_pam_reads * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 5] <- tmp_out
        tmp_out <- object@stats$per_unmapped_reads * 100
        tmp_out <- sapply(tmp_out, function(x) round(x, 1))
        df_outs[, 6] <- tmp_out
        df_outs[, 7] <- object@cutoffs$library_percent * 100
        df_outs[, 8] <- object@stats$qcpass_library_per

        df_outs <- df_outs[match(mixedsort(df_outs$Sample), df_outs$Sample), ]

        if (length(out_dir) == 0) {
            reactable(df_outs, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,
                      theme = reactableTheme(
                          style = list(fontFamily = "-apple-system", fontSize = "0.75rem")),
                      columns = list(
                          "Group" = colDef(minWidth = 100),
                          "Sample" = colDef(minWidth = 100),
                          "% Library Reads" = colDef(style = function(value) {
                                                                 if (value < object@cutoffs$library_percent * 100) {
                                                                    color <- "red"
                                                                    fweight <- "bold"
                                                                } else {
                                                                    color <- "forestgreen"
                                                                    fweight <- "plain"
                                                                }
                                                                list(color = color, fontWeight = fweight)}),
                          "Pass" = colDef(cell = function(value) {
                                                   if (value) "\u2705" else "\u274c" }))
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
#' @param object   sampleQC object
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_libcov",
    signature = "sampleQC",
    definition = function(object,
                          out_dir = NULL) {
        cols <- c("Group",
                  "Sample",
                  "Total Library Reads",
                  "Total Template Oligo Sequences",
                  "Library Coverage",
                  "Pass Threshold",
                  "Pass")
        df_outs <- data.frame(matrix(NA, nrow(object@stats), length(cols)))
        colnames(df_outs) <- cols

        df_outs[, 1] <- object@samples[[1]]@libname
        df_outs[, 2] <- rownames(object@stats)
        df_outs[, 3] <- object@stats$library_reads
        df_outs[, 4] <- object@stats$library_seqs
        tmp_out <- object@stats$library_cov
        tmp_out <- sapply(tmp_out, function(x) round(x, 0))
        df_outs[, 5] <- tmp_out
        df_outs[, 6] <- object@cutoffs$library_cov
        df_outs[, 7] <- object@stats$qcpass_library_cov

        df_outs <- df_outs[match(mixedsort(df_outs$Sample), df_outs$Sample), ]

        if (length(out_dir) == 0) {
            reactable(df_outs, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,
                      theme = reactableTheme(
                          style = list(fontFamily = "-apple-system", fontSize = "0.75rem")),
                      columns = list(
                          "Group" = colDef(minWidth = 100),
                          "Sample" = colDef(minWidth = 100),
                          "Total Library Reads" = colDef(format = colFormat(separators = TRUE)),
                          "Total Template Oligo Sequences" = colDef(format = colFormat(separators = TRUE)),
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
                                                   if (value) "\u2705" else "\u274c" }))
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

        cols <- c("Group",
                  "Sample",
                  "Chromosome",
                  "Strand",
                  "Genomic Start",
                  "Genomic End",
                  "% Low Abundance",
                  "Low Abundance cutoff",
                  "Pass Threshold",
                  "Pass")
        df_outs <- data.frame(matrix(NA, nrow(object@stats), length(cols)))
        colnames(df_outs) <- cols

        df_outs[, 1] <- object@samples[[1]]@libname
        df_outs[, 2] <- rownames(object@stats)
        df_outs[, 3] <- sapply(object@library_counts_chr, function (x) x[[1]])
        df_outs[, 4] <- sapply(object@library_counts_chr, function (x) x[[2]])
        df_outs[, 5] <- sapply(object@library_counts_chr, function (x) x[[3]])
        df_outs[, 6] <- sapply(object@library_counts_chr, function (x) x[[4]])

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
        df_outs[, 7] <- low_per

        df_outs[, 8] <- object@cutoffs$seq_low_count
        df_outs[, 9] <- (1 - object@cutoffs$low_abundance_lib_per) * 100
        df_outs[, 10] <- df_outs[, 7] < (1 - object@cutoffs$low_abundance_lib_per) * 100

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
                                                   if (value) "\u2705" else "\u274c" }))
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
#' @param object   sampleQC object
#' @param out_dir  the output directory
setMethod(
    "qcout_samqc_pos_anno",
    signature = "sampleQC",
    definition = function(object,
                          out_dir = NULL) {
        cols <- c("Group",
                  "Sample",
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
        df_outs[, 3] <- sapply(object@library_counts_chr, function (x) x[[1]])
        df_outs[, 4] <- sapply(object@library_counts_chr, function (x) x[[2]])
        df_outs[, 5] <- sapply(object@library_counts_chr, function (x) x[[3]])
        df_outs[, 6] <- sapply(object@library_counts_chr, function (x) x[[4]])

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

        df_outs[, 7] <- lof_low_per
        df_outs[, 8] <- others_low_per
        df_outs[, 9] <- lof_low_per + others_low_per
        df_outs[, 10] <- object@cutoffs$low_abundance_per * 100
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
                                                   if (value) "\u2705" else "\u274c" }))
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
