#' initialize function
setGeneric("run_sample_qc", function(object, ...) {
  standardGeneric("run_sample_qc")
})

#' run sample QC for the list of samples
#'
#' @export
#' @param object                 sampleQC object
#' @param qc_type                plasmid or screen
#' @param cutoff_total           qc cutoff of the total reads
#' @param cutoff_missing_per     qc cutoff of the missing variant percentage
#' @param cutoff_low_count       count cutoff of library reads
#' @param cutoff_low_sample_per  sample percentage cutoff of library reads
#' @param cutoff_accepted        qc cutoff of the total accepted reads
#' @param cutoff_mapping_per     qc cutoff of mapping percentage (ref + pam + library)
#' @param cutoff_ref_per         qc cutoff of reference percentage
#' @param cutoff_library_per     qc cutoff of library reads percentage
#' @param cutoff_library_cov     qc cutoff of library coverage
#' @param cutoff_low_per         qc cutoff of low abundance percentage for LOF plot
#' @param cutoff_low_lib_per     qc cutoff of the percentage of library sequences with low abundance for LOF plot
#' @return object
setMethod(
    "run_sample_qc",
    signature = "sampleQC",
    definition = function(object,
                          qc_type = c("plasmid", "screen"),
                          cutoff_total = maveqc_config$sqc_total,
                          cutoff_missing_per = maveqc_config$sqc_missing,
                          cutoff_low_count = maveqc_config$sqc_low_count,
                          cutoff_low_sample_per = maveqc_config$sqc_low_sample_per,
                          cutoff_accepted = maveqc_config$sqc_accepted,
                          cutoff_mapping_per = maveqc_config$sqc_mapping_per,
                          cutoff_ref_per = maveqc_config$sqc_ref_per,
                          cutoff_library_per = maveqc_config$sqc_library_per,
                          cutoff_library_cov = maveqc_config$sqc_library_cov,
                          cutoff_low_per = maveqc_config$sqc_low_per,
                          cutoff_low_lib_per = maveqc_config$sqc_low_lib_per) {
        #----------#
        # checking #
        #----------#
        if (length(object@samples) == 0) {
            stop(paste0("====> Error: no sample found in the sampleQC object!"))
        }

        qc_type <- match.arg(qc_type)
        if (qc_type == "screen") {
            if (length(object@samples_ref) == 0) {
                stop(paste0("====> Error: samples_ref is empty! Screen QC must have reference samples."))
            }
        }

        cols <- c("num_total_reads",
                  "per_missing_variants",
                  "seq_low_count",
                  "seq_low_sample_per",
                  "num_accepted_reads",
                  "per_mapping_reads",
                  "per_ref_reads",
                  "per_library_reads",
                  "library_cov",
                  "low_abundance_per",
                  "low_abundance_lib_per")
        df_cutoffs <- data.frame(matrix(NA, 1, length(cols)))
        colnames(df_cutoffs) <- cols

        df_cutoffs$num_total_reads <- cutoff_total
        df_cutoffs$per_missing_variants <- cutoff_missing_per
        df_cutoffs$seq_low_count <- cutoff_low_count
        df_cutoffs$seq_low_sample_per <- cutoff_low_sample_per
        df_cutoffs$num_accepted_reads <- cutoff_accepted
        df_cutoffs$per_mapping_reads <- cutoff_mapping_per
        df_cutoffs$per_ref_reads <- cutoff_ref_per
        df_cutoffs$per_library_reads <- cutoff_library_per
        df_cutoffs$library_cov <- cutoff_library_cov
        df_cutoffs$low_abundance_per <- cutoff_low_per
        df_cutoffs$low_abundance_lib_per <- cutoff_low_lib_per

        object@cutoffs <- df_cutoffs

        #-------------------------------------------#
        # 1. Filtering by the total number of reads #
        #-------------------------------------------#
        cat("Filtering by the total number of reads...", "\n", sep = "")

        sample_names <- character()
        for (s in object@samples) {
            sample_names <- append(sample_names, s@sample)

            object@stats[s@sample, ]$per_r1_adaptor <- s@per_r1_adaptor
            object@stats[s@sample, ]$per_r2_adaptor <- s@per_r2_adaptor

            object@stats[s@sample, ]$total_reads <- s@allstats$total_counts
            object@stats[s@sample, ]$ref_reads <- s@allstats_qc$num_ref_reads
            object@stats[s@sample, ]$pam_reads <- s@allstats_qc$num_pam_reads

            object@stats[s@sample, ]$gini_coeff_before_qc <- s@libstats_qc$gini_coeff
        }

        #---------------------------------------#
        # 2. Filtering by low counts            #
        #    a) k-means clustering on screen QC #
        #    a) hard cutoff on Plasmid QC       #
        #---------------------------------------#
        cat("Filtering by low counts...", "\n", sep = "")

        if (qc_type == "screen") {
            cat("    |--> Creating k-means clusters...", "\n", sep = "")

            ref_counts <- data.table()
            for (s in object@samples_ref) {
                tmp_counts <- s@allcounts[, c("sequence", "count")]
                tmp_counts <- as.data.table(tmp_counts)

                if (nrow(ref_counts) == 0) {
                    ref_counts <- tmp_counts
                    colnames(ref_counts) <- c("sequence", s@sample)
                } else {
                    tmp_cols <- colnames(ref_counts)
                    ref_counts <- merge(ref_counts, tmp_counts, by = "sequence", all = TRUE)
                    colnames(ref_counts) <- c(tmp_cols, s@sample)
                }
            }
            ref_counts[, count := rowSums(.SD, na.rm = TRUE), .SDcols = 2:ncol(ref_counts)]
            ref_counts[, count_log2 := log2(count + 1)]

            kmeans_res <- Ckmeans.1d.dp(ref_counts$count_log2, k = 2, y = 1)
            ref_counts$cluster <- kmeans_res$cluster

            cat("    |--> Filtering using clusters...", "\n", sep = "")

            # filtering sequences on input samples by filtered set
            for (s in object@samples) {
                cat("        |--> Filtering on ", s@sample, "\n", sep = "")

                unfiltered_counts <- s@allcounts[, c("sequence", "count")]
                object@seq_clusters[[s@sample]] <- ref_counts

                # considering missing seqs
                object@accepted_counts[[s@sample]] <- unfiltered_counts[ref_counts[cluster == 2], on = .(sequence), nomatch = 0]

                object@bad_seqs_bycluster[[s@sample]] <- unfiltered_counts[!ref_counts[cluster == 2], on = .(sequence)]
            }
        } else {
            cat("    |--> Creating k-means clusters...", "\n", sep = "")

            for (s in object@samples) {
                cat("        |--> Filtering on ", s@sample, "\n", sep = "")

                tmp_counts <- s@allcounts[, c("sequence", "count")]
                tmp_counts <- as.data.table(tmp_counts)
                tmp_counts[, count_log2 := log2(count + 1)]

                kmeans_res <- Ckmeans.1d.dp(tmp_counts$count_log2, k = 2, y = 1)
                tmp_counts$cluster <- kmeans_res$cluster

                object@seq_clusters[[s@sample]] <- tmp_counts

                object@accepted_counts[[s@sample]] <- tmp_counts[cluster == 2]

                object@bad_seqs_bycluster[[s@sample]] <- tmp_counts[cluster == 1, c("sequence", "count")]
            }
        }

        #-------------------------------------#
        # 3. Filtering by depth and samples   #
        #    a) count >= X                    #
        #    b) in >= X% of samples           #
        #-------------------------------------#

        # if plasmid qc, don't apply percentage filtering as samples have different seqs
        if (qc_type == "screen") {
            cat("Filtering by depth and percentage in samples...", "\n", sep = "")

            # note library independent counts may have different sequences
            accepted_counts <- merge_list_to_dt(object@accepted_counts, "sequence", "count")
            accepted_counts[, sample_number := rowSums(.SD >= cutoff_low_count, na.rm = TRUE), .SDcols = 2:ncol(accepted_counts)]
            accepted_counts[, sample_percentage := sample_number / length(sample_names)]

            for (s in object@samples) {
                cols <- c("sequence", s@sample)
                object@accepted_counts[[s@sample]] <- na.omit(accepted_counts[sample_percentage >= cutoff_low_sample_per, ..cols], cols = s@sample)
                colnames(object@accepted_counts[[s@sample]]) <- c("sequence", "count")

                cols <- c("sequence", s@sample)
                object@bad_seqs_bydepth[[s@sample]] <- na.omit(accepted_counts[sample_percentage < cutoff_low_sample_per, ..cols], cols = s@sample)
                colnames(object@bad_seqs_bydepth[[s@sample]]) <- c("sequence", "count")
            }
        } else {
            cat("Filtering by depth...", "\n", sep = "")

            for (s in object@samples) {
                tmp_counts <- object@accepted_counts[[s@sample]]
                cols <- c("sequence", "count")

                object@accepted_counts[[s@sample]] <- tmp_counts[count >= cutoff_low_count, ..cols]

                object@bad_seqs_bydepth[[s@sample]] <- tmp_counts[count < cutoff_low_count, ..cols]
            }
        }

        #--------------------------------------#
        # 4. Filtering by library mapping      #
        #    a) reads mapped to VaLiAnT output #
        #--------------------------------------#
        cat("Filtering by library mapping...", "\n", sep = "")

        for (s in object@samples) {
            tmp_counts <- object@accepted_counts[[s@sample]]

            # meta_mseqs without ref and pam by format_count
            # and using library dependent sequences instead of meta sequences
            # accepted_counts have ref and pam
            object@library_counts[[s@sample]] <- object@accepted_counts[[s@sample]][sequence %in% s@meta_mseqs]
            object@unmapped_counts[[s@sample]] <- object@accepted_counts[[s@sample]][sequence %nin% c(s@meta_mseqs, s@refseq, s@pamseq)]

            object@bad_seqs_bylib[[s@sample]] <- object@unmapped_counts[[s@sample]]
        }

        #--------------------------------------#
        # 5. Filtering by library coverage     #
        #    a) library reads / oligos in meta #
        #--------------------------------------#
        cat("Filtering by library coverage...", "\n", sep = "")

        for (s in object@samples) {
            object@stats[s@sample, ]$accepted_reads <- sum(object@accepted_counts[[s@sample]]$count, na.rm = TRUE)
            object@stats[s@sample, ]$excluded_reads <- object@stats[s@sample, ]$total_reads - object@stats[s@sample, ]$accepted_reads
            object@stats[s@sample, ]$library_reads <- sum(object@library_counts[[s@sample]]$count, na.rm = TRUE)
            object@stats[s@sample, ]$unmapped_reads <- sum(object@unmapped_counts[[s@sample]]$count, na.rm = TRUE)

            object@stats[s@sample, ]$per_library_reads <- object@stats[s@sample, ]$library_reads / object@stats[s@sample, ]$accepted_reads
            object@stats[s@sample, ]$per_unmapped_reads <- object@stats[s@sample, ]$unmapped_reads / object@stats[s@sample, ]$accepted_reads
            object@stats[s@sample, ]$per_ref_reads <- object@stats[s@sample, ]$ref_reads / object@stats[s@sample, ]$accepted_reads
            object@stats[s@sample, ]$per_pam_reads <- object@stats[s@sample, ]$pam_reads / object@stats[s@sample, ]$accepted_reads

            object@stats[s@sample, ]$missing_meta_seqs <- length(s@missing_meta_seqs)
            object@stats[s@sample, ]$per_missing_meta_seqs <- length(s@missing_meta_seqs) / length(s@meta_mseqs)

            object@stats[s@sample, ]$library_seqs <- length(s@meta_mseqs)
            object@stats[s@sample, ]$median_cov <- median(object@library_counts[[s@sample]]$count)
            object@stats[s@sample, ]$library_cov <- as.integer(object@stats[s@sample, ]$library_reads / length(s@meta_mseqs))
        }

        #--------------------- -----------------#
        # 6. Sorting library counts by position #
        #---------------------------------------#
        cat("Sorting library counts by position...", "\n", sep = "")

        # main issue:
        # meta seqs have adaptors
        # oligo names are not unique
        # a seq has a consequence, but has many names, don't know which name is right

        for (s in object@samples) {
            cat("    |--> Sorting on ", s@sample, "\n", sep = "")

            tmp_meta <- s@valiant_meta[, c("oligo_name", "mut_position")]

            # fecth oligo name using library dependent counts
            tmp_map <- s@libcounts[, c("name", "sequence")]

            libcounts_pos <- object@library_counts[[s@sample]]
            libcounts_pos[tmp_map, oligo_name := i.name, on = .(sequence)]
            libcounts_pos[tmp_meta, position := i.mut_position, on = .(oligo_name)]
            setorder(libcounts_pos, cols = "position")

            # in case some sequences don't match in the library
            libcounts_pos <- na.omit(libcounts_pos, cols = "position")

            object@library_counts_pos[[s@sample]] <- libcounts_pos
            object@library_counts_chr[[s@sample]] <- c(unique(s@valiant_meta$ref_chr),
                                                       unique(s@valiant_meta$ref_strand),
                                                       unique(s@valiant_meta$ref_start),
                                                       unique(s@valiant_meta$ref_end))
        }

        #------------------------#
        # 7. Gini coeff after qc #
        #------------------------#
        cat("Calculating gini coefficiency...", "\n", sep = "")

        for (s in object@samples) {
            gini_coeff <- cal_gini(object@library_counts[[s@sample]]$count, corr = FALSE, na.rm = TRUE)
            object@stats[s@sample, ]$gini_coeff_after_qc <- round(gini_coeff, 3)
        }

        #------------------#
        # 8. QC results    #
        #------------------#
        object@stats$per_missing_meta_seqs <- unlist(lapply(object@stats$per_missing_meta_seqs, function(x) round(x, 4)))
        object@stats$per_library_reads <- unlist(lapply(object@stats$per_library_reads, function(x) round(x, 4)))
        object@stats$per_unmapped_reads <- unlist(lapply(object@stats$per_unmapped_reads, function(x) round(x, 4)))
        object@stats$per_ref_reads <- unlist(lapply(object@stats$per_ref_reads, function(x) round(x, 4)))
        object@stats$per_pam_reads <- unlist(lapply(object@stats$per_pam_reads, function(x) round(x, 4)))
        object@stats$per_r1_adaptor <- unlist(lapply(object@stats$per_r1_adaptor, function(x) round(x, 4)))
        object@stats$per_r2_adaptor <- unlist(lapply(object@stats$per_r2_adaptor, function(x) round(x, 4)))

        object@stats$qcpass_total_reads <- unlist(lapply(object@stats$total_reads, function(x) ifelse(x >= cutoff_total, TRUE, FALSE)))
        object@stats$qcpass_missing_per <- unlist(lapply(object@stats$per_missing_meta_seqs, function(x) ifelse(x < cutoff_missing_per, TRUE, FALSE)))
        object@stats$qcpass_accepted_reads <- unlist(lapply(object@stats$accepted_reads, function(x) ifelse(x >= cutoff_accepted, TRUE, FALSE)))
        object@stats$qcpass_mapping_per <- unlist(lapply(object@stats$per_unmapped_reads, function(x) ifelse(x < (1 - cutoff_mapping_per), TRUE, FALSE)))
        object@stats$qcpass_ref_per <- unlist(lapply(object@stats$per_ref_reads, function(x) ifelse(x < cutoff_ref_per, TRUE, FALSE)))
        object@stats$qcpass_library_per <- unlist(lapply(object@stats$per_library_reads, function(x) ifelse(x >= cutoff_library_per, TRUE, FALSE)))
        object@stats$qcpass_library_cov <- unlist(lapply(object@stats$library_cov, function(x) ifelse(x >= cutoff_library_cov, TRUE, FALSE)))

        qc_lables <- c("qcpass_total_reads",
                       "qcpass_missing_per",
                       "qcpass_accepted_reads",
                       "qcpass_mapping_per",
                       "qcpass_ref_per",
                       "qcpass_library_per",
                       "qcpass_library_cov")
        object@stats$qcpass <- apply(object@stats[, qc_lables], 1, function(x) all(x))

        #------------------------#
        # 9. Filtered samples    #
        #------------------------#

        object@filtered_samples <- rownames(object@stats[object@stats$qcpass == TRUE, ])

        #-------------------------#
        # 10. map vep consequence #
        #-------------------------#

        # main issue:
        # oligo names are not unique to sequence, cannot use as reference
        # a seq has a consequence, but has many names, don't which name is right
        # so, vep must have right seq, otherwise cannot determine the consequence

        # if plasmid qc, don't apply
        if (qc_type == "screen") {
            cat("Mapping consequencing annotation...", "\n", sep = "")

            # assuming all the samples have the sample library sequences and corresponding consequences
            # using sequences in vep annotation to identify consequences
            # assuming unique_oligo_name in vep_anno is unique
            vep_anno <- object@samples[[1]]@vep_anno[, c("unique_oligo_name", "seq", "summary_plot")]
            colnames(vep_anno) <- c("oligo_name", "sequence", "consequence")

            # merge counts, but use data table in case rownames of data frame is not unique
            library_counts_anno <- merge_list_to_dt(object@library_counts, "sequence", "count")
            cols <- colnames(library_counts_anno)
            library_counts_anno[vep_anno, c("oligo_name", "consequence") := .(oligo_name, consequence), on = .(sequence)]
            setcolorder(library_counts_anno, c("oligo_name", "consequence", cols))
            object@library_counts_anno <- library_counts_anno

            # merge all the sorted library counts
            library_counts_pos_anno <- data.table()
            for (s in object@samples) {
                tmp_pos <- object@library_counts_pos[[s@sample]][, c("sequence", "position", "count")]
                colnames(tmp_pos) <- c("sequence", "position", s@sample)

                if (nrow(library_counts_pos_anno) == 0) {
                    library_counts_pos_anno <- tmp_pos
                } else {
                    library_counts_pos_anno <- merge(library_counts_pos_anno, tmp_pos, by = c("sequence", "position"), all = TRUE)
                }
            }

            cols <- colnames(library_counts_pos_anno)[-2]
            library_counts_pos_anno[vep_anno, c("oligo_name", "consequence") := .(oligo_name, consequence), on = .(sequence)]
            setcolorder(library_counts_pos_anno, c("oligo_name", "consequence", "position", cols))
            object@library_counts_pos_anno <- library_counts_pos_anno
        }

        return(object)
    }
)
