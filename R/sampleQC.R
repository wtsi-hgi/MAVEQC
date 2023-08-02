#' initialize function
setGeneric("run_sample_qc", function(object, ...) {
  standardGeneric("run_sample_qc")
})

#' run sample QC for the list of samples
#'
#' @export
#' @param object                 sampleQC object
#' @param qc_type                plasmid or screen
#' @param cutoff_low_count       count cutoff of library reads
#' @param cutoff_low_sample_per  sample percentage cutoff of library reads
#' @param cutoff_filtered        qc cutoff of the total filtered reads
#' @param cutoff_mapping_per     qc cutoff of mapping percentage (ref + pam + library)
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
                          cutoff_low_count = maveqc_config$sqc_low_count,
                          cutoff_low_sample_per = maveqc_config$sqc_low_sample_per,
                          cutoff_filtered = maveqc_config$sqc_filtered,
                          cutoff_mapping_per = maveqc_config$sqc_mapping_per,
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

        cols <- c("total_reads",
                  "seq_low_count",
                  "seq_low_sample_per",
                  "library_percent",
                  "library_cov",
                  "low_abundance_per",
                  "low_abundance_lib_per")
        df_cutoffs <- data.frame(matrix(NA, 1, length(cols)))
        colnames(df_cutoffs) <- cols

        df_cutoffs$total_reads <- cutoff_filtered
        df_cutoffs$seq_low_count <- cutoff_low_count
        df_cutoffs$seq_low_sample_per <- cutoff_low_sample_per
        df_cutoffs$library_percent <- cutoff_library_per
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

                unfiltered_counts <- s@allcounts$count
                names(unfiltered_counts) <- s@allcounts$sequence

                object@seq_clusters[[s@sample]] <- ref_counts
                object@accepted_seqs[[s@sample]] <- ref_counts[cluster == 2]$sequence

                # considering missing seqs
                name_check <- names(unfiltered_counts) %in% object@accepted_seqs[[s@sample]]
                object@accepted_counts[[s@sample]] <- unfiltered_counts[name_check]

                object@bad_seqs_bycluster[[s@sample]] <- ref_counts[cluster == 1]$sequence
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
                object@accepted_seqs[[s@sample]] <- tmp_counts[cluster == 2]$sequence

                object@accepted_counts[[s@sample]] <- tmp_counts[cluster == 2]$count
                names(object@accepted_counts[[s@sample]]) <- tmp_counts[cluster == 2]$sequence

                object@bad_seqs_bycluster[[s@sample]] <- tmp_counts[cluster == 1]$sequence
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
            accepted_counts <- merge_list_to_df(object@accepted_counts)
            accepted_counts_depth <- cbind(accepted_counts, rowSums(accepted_counts >= cutoff_low_count, na.rm = TRUE))
            accepted_counts_depth <- cbind(accepted_counts_depth, accepted_counts_depth[, ncol(accepted_counts_depth)] / length(sample_names))
            colnames(accepted_counts_depth) <- c(sample_names, "sample_number", "sample_percentage")

            accepted_counts_final <- accepted_counts_depth[accepted_counts_depth$sample_percentage >= cutoff_low_sample_per, sample_names]

            for (s in object@samples) {
                tmp_counts <- accepted_counts_final[, s@sample]
                names(tmp_counts) <- rownames(accepted_counts_final)
                object@accepted_counts[[s@sample]] <- tmp_counts[!is.na(tmp_counts)]

                f_per <- accepted_counts_depth$sample_percentage < cutoff_low_sample_per
                f_na <- !is.na(accepted_counts_depth[, s@sample])
                accepted_counts_bad <- accepted_counts_depth[f_per & f_na, s@sample, drop = FALSE]
                object@bad_seqs_bydepth[[s@sample]] <- rownames(accepted_counts_bad)
            }
        } else {
            cat("Filtering by depth...", "\n", sep = "")

            for (s in object@samples) {
                tmp_counts <- object@accepted_counts[[s@sample]]
                object@accepted_counts[[s@sample]] <- tmp_counts[tmp_counts >= cutoff_low_count]

                object@bad_seqs_bydepth[[s@sample]] <- names(tmp_counts[tmp_counts < cutoff_low_count])
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
            # accepted_counts have ref and pam
            object@library_counts[[s@sample]] <- tmp_counts[names(tmp_counts) %in% s@meta_mseqs]
            object@unmapped_counts[[s@sample]] <- tmp_counts[names(tmp_counts) %nin% c(s@meta_mseqs, s@refseq, s@pamseq)]

            object@bad_seqs_bylib[[s@sample]] <- names(object@unmapped_counts[[s@sample]])
        }

        #--------------------------------------#
        # 5. Filtering by library coverage     #
        #    a) library reads / oligos in meta #
        #--------------------------------------#
        cat("Filtering by library coverage...", "\n", sep = "")

        for (s in object@samples) {
            sample_name <- s@sample
            object@stats[sample_name, ]$accepted_reads <- sum(object@accepted_counts[[sample_name]], na.rm = TRUE)
            object@stats[sample_name, ]$excluded_reads <- object@stats[sample_name, ]$total_reads - object@stats[sample_name, ]$accepted_reads
            object@stats[sample_name, ]$library_reads <- sum(object@library_counts[[sample_name]], na.rm = TRUE)
            object@stats[sample_name, ]$unmapped_reads <- sum(object@unmapped_counts[[sample_name]], na.rm = TRUE)

            object@stats[sample_name, ]$per_library_reads <- object@stats[sample_name, ]$library_reads / object@stats[sample_name, ]$accepted_reads
            object@stats[sample_name, ]$per_unmapped_reads <- object@stats[sample_name, ]$unmapped_reads / object@stats[sample_name, ]$accepted_reads
            object@stats[sample_name, ]$per_ref_reads <- object@stats[sample_name, ]$ref_reads / object@stats[sample_name, ]$accepted_reads
            object@stats[sample_name, ]$per_pam_reads <- object@stats[sample_name, ]$pam_reads / object@stats[sample_name, ]$accepted_reads

            object@stats[sample_name, ]$missing_meta_seqs <- length(s@missing_meta_seqs)

            object@stats[sample_name, ]$library_seqs <- length(s@meta_mseqs)
            object@stats[sample_name, ]$library_cov <- as.integer(object@stats[sample_name, ]$library_reads / length(s@meta_mseqs))
        }

        #--------------------- -----------------#
        # 6. Sorting library counts by position #
        #---------------------------------------#
        cat("Sorting library counts by position...", "\n", sep = "")

        # main issue:
        # meta seqs are not right
        # oligo names are not unique
        # a seq has a consequence, but has many names, don't know which name is right

        for (s in object@samples) {
            cat("    |--> Sorting on ", s@sample, "\n", sep = "")

            tmp_meta <- s@valiant_meta[, c("oligo_name", "mut_position")]
            tmp_meta <- as.data.table(tmp_meta)

            # fecth oligo name using library dependent counts
            tmp_map <- s@libcounts[, c("name", "sequence")]
            colnames(tmp_map) <- c("oligo_name", "seq")
            tmp_map <- as.data.table(tmp_map)

            tmp_lib <- object@library_counts[[s@sample]]
            libcounts_pos <- cbind(tmp_lib, rep(0, length(tmp_lib)))
            libcounts_pos <- as.data.frame(libcounts_pos)
            colnames(libcounts_pos) <- c("count", "position")
            libcounts_pos$seq <- names(tmp_lib)
            libcounts_pos <- as.data.table(libcounts_pos)

            # multiple oligo names have the same seq and the same position
            # so here is fine, but need to refine
            # oligo name may be not right, but position is good
            libcounts_pos[tmp_map, oligo_name := i.oligo_name, on = .(seq)]
            libcounts_pos[tmp_meta, position := i.mut_position, on = .(oligo_name)]
            setorder(libcounts_pos, cols = "position")

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
            gini_coeff <- cal_gini(object@library_counts[[s@sample]], corr = FALSE, na.rm = TRUE)
            object@stats[s@sample, ]$gini_coeff_after_qc <- round(gini_coeff, 3)
        }

        #------------------#
        # 8. QC results    #
        #------------------#

        object@stats$qcpass_accepted_reads <- unlist(lapply(object@stats$accepted_reads, function(x) ifelse(x >= cutoff_filtered, TRUE, FALSE)))
        object@stats$qcpass_mapping_per <- unlist(lapply(object@stats$per_unmapped_reads, function(x) ifelse(x < (1 - cutoff_mapping_per), TRUE, FALSE)))
        object@stats$qcpass_library_per <- unlist(lapply(object@stats$per_library_reads, function(x) ifelse(x >= cutoff_library_per, TRUE, FALSE)))
        object@stats$qcpass_library_cov <- unlist(lapply(object@stats$library_cov, function(x) ifelse(x >= cutoff_library_cov, TRUE, FALSE)))

        qc_lables <- c("qcpass_accepted_reads", "qcpass_mapping_per", "qcpass_library_per", "qcpass_library_cov")
        object@stats$qcpass <- apply(object@stats[, qc_lables], 1, function(x) all(x))

        #------------------------#
        # 9. Filtered samples    #
        #------------------------#

        object@filtered_samples <- rownames(object@stats[object@stats$qcpass == TRUE, ])

        #-------------------------#
        # 10. map vep consequence #
        #-------------------------#

        # main issue:
        # meta_consequence seqs are reverse complement, fix it in create_sge_object
        # but it may change later
        # oligo names are not unique to sequence, cannot use as reference
        # a seq has a consequence, but has many names, don't which name is right
        # so, vep must have right seq, otherwise cannot determine the consequence

        # if plasmid qc, don't apply
        if (qc_type == "screen") {
            cat("Mapping consequencing annotation...", "\n", sep = "")

            # assuming all the samples have the sample library sequences and corresponding consequences
            # using sequences in vep annotation to identify consequences
            vep_anno <- object@samples[[1]]@vep_anno[, c("unique_oligo_name", "seq", "summary_plot")]
            colnames(vep_anno) <- c("oligo_name", "seq", "consequence")
            vep_anno <- as.data.table(vep_anno)

            # merge counts, but use data table in case rownames of data frame is not unique
            library_counts_anno <- merge_list_to_df(object@library_counts)
            library_counts_anno[is.na(library_counts_anno)] <- 0

            library_counts_anno$seq <- rownames(library_counts_anno)
            library_counts_anno <- as.data.table(library_counts_anno)

            library_counts_anno[vep_anno, consequence := i.consequence, on = .(seq)]
            object@library_counts_anno <- library_counts_anno

            # merge all the sorted library counts
            library_counts_pos_anno <- data.table()
            for (s in object@samples) {
                tmp_pos <- object@library_counts_pos[[s@sample]]
                tmp_pos <- tmp_pos[, c("seq", "position", "count")]
                colnames(tmp_pos) <- c("seq", "position", s@sample)

                if (nrow(library_counts_pos_anno) == 0) {
                    library_counts_pos_anno <- tmp_pos
                } else {
                    library_counts_pos_anno <- merge(library_counts_pos_anno, tmp_pos, by = c("seq", "position"), all = TRUE)
                }
            }

            library_counts_pos_anno[vep_anno, consequence := i.consequence, on = .(seq)]
            object@library_counts_pos_anno <- library_counts_pos_anno
        }

        return(object)
    }
)
