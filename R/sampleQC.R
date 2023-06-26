#' initialize function
setGeneric("run_sample_qc", function(object, ...) {
  standardGeneric("run_sample_qc")
})

#' run sample QC for the list of samples
#'
#' @export
#' @param object          sampleQC object
#' @param qc_type         plasmid or screen
#' @param effect_count    count cutoff of effective reads
#' @param effect_per      sample percentage cutoff of effective reads
#' @param cutoff_filtered qc cutoff of the total filtered reads
#' @param cutoff_mapping  qc cutoff of mapping percentage (ref + pam + effect)
#' @param cutoff_effect   qc cutoff of effective reads percentage
#' @param cutoff_cov      qc cutoff of effective coverage
#' @return object
setMethod(
    "run_sample_qc",
    signature = "sampleQC",
    definition = function(object,
                          qc_type,
                          effect_count = 5,
                          effect_per = 0.25,
                          cutoff_filtered = 1000000,
                          cutoff_mapping = 0.6,
                          cutoff_effect = 0.4,
                          cutoff_cov = 100) {
        #----------#
        # checking #
        #----------#
        if (length(object@samples) == 0) {
            stop(paste0("====> Error: no sample found in the sampleQC object!"))
        }

        if (length(qc_type) == 0) {
            stop(paste0("====> Error: please provide QC type."))
        } else {
            if (qc_type %in% c("plasmid", "screen")) {
                if (qc_type == "screen") {
                    if (length(object@samples_ref) == 0) {
                        stop(paste0("====> Error: samples_ref is empty! Screen QC must have reference samples."))
                    }
                }
            } else {
                stop(paste0("====> Error: wrong QC type! Please use plasmid or screen."))
            }
        }

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

            # merging reference counts, data.table() to speed up
            ref_counts <- data.table()
            for (s in object@samples_ref) {
                tmp_counts <- s@allcounts[, c("sequence", "count")]
                tmp_counts <- as.data.table(tmp_counts)

                if (nrow(ref_counts) == 0) {
                    ref_counts <- tmp_counts
                } else {
                    ref_counts <- merge(ref_counts, tmp_counts, by = "sequence", all = TRUE)
                    tmp_cols <- vector()
                    for (i in 1:(ncol(ref_counts) - 1)) {
                        tmp_cols <- c(tmp_cols, paste0("ref", i))
                    }
                    colnames(ref_counts) <- c("sequence", tmp_cols)
                }
            }
            ref_counts <- as.data.frame(ref_counts)
            rownames(ref_counts) <- ref_counts$sequence
            ref_counts <- subset(ref_counts, select = -sequence)

            ref_counts_merged <- rowSums(ref_counts, na.rm = TRUE)
            ref_counts_merged_log2 <- log2(ref_counts_merged + 1)

            # create filtered set of sequences by k-means clustering
            kmeans_res <- Ckmeans.1d.dp(ref_counts_merged_log2, k = 2, y = 1)
            ref_clusters <- cbind(ref_counts_merged, ref_counts_merged_log2, kmeans_res$cluster)
            colnames(ref_clusters) <- c("count", "count_log2", "cluster")
            ref_clusters <- data.frame(ref_clusters)

            object@seq_clusters[["ref"]] <- ref_clusters
            object@filtered_seqs <- rownames(ref_clusters[ref_clusters$cluster == 2, ])
            object@bad_seqs$by_cluster <- rownames(ref_clusters[ref_clusters$cluster == 1, ])

            cat("    |--> Filtering using clusters...", "\n", sep = "")

            # filtering sequences on input samples by filtered set
            unfiltered_counts <- data.table()
            for (s in object@samples) {
                cat("        |--> Filtering on ", s@sample, "\n", sep = "")

                tmp_counts <- s@allcounts[, c("sequence", "count")]
                tmp_counts <- as.data.table(tmp_counts)

                if (nrow(unfiltered_counts) == 0) {
                    unfiltered_counts <- tmp_counts
                } else {
                    unfiltered_counts <- merge(unfiltered_counts, tmp_counts, by = "sequence", all = TRUE)
                    tmp_cols <- vector()
                    for (i in 1:(ncol(unfiltered_counts) - 1)) {
                        tmp_cols <- c(tmp_cols, paste0("un", i))
                    }
                    colnames(unfiltered_counts) <- c("sequence", tmp_cols)
                }
            }
            unfiltered_counts <- as.data.frame(unfiltered_counts)
            rownames(unfiltered_counts) <- unfiltered_counts$sequence
            unfiltered_counts <- subset(unfiltered_counts, select = -sequence)
            colnames(unfiltered_counts) <- sample_names

            filtered_counts <- unfiltered_counts[object@filtered_seqs, ]
        } else {
            filtered_counts <- data.table()
            for (s in object@samples) {
                tmp_counts <- s@allcounts$count
                names(tmp_counts) <- s@allcounts$sequence

                tmp_counts_log2 <- log2(tmp_counts + 1)
                kmeans_res <- Ckmeans.1d.dp(tmp_counts_log2, k = 2, y = 1)
                tmp_clusters <- cbind(tmp_counts, tmp_counts_log2, kmeans_res$cluster)
                colnames(tmp_clusters) <- c("count", "count_log2", "cluster")
                tmp_clusters <- data.frame(tmp_clusters)

                object@seq_clusters[[s@sample]] <- tmp_clusters

                tmp_counts_filtered <- tmp_clusters[tmp_clusters$cluster == 2, "count", drop = FALSE]
                tmp_counts_filtered$sequence <- rownames(tmp_counts_filtered)
                tmp_counts_filtered <- as.data.table(tmp_counts_filtered)

                if (nrow(filtered_counts) == 0) {
                    filtered_counts <- tmp_counts_filtered
                } else {
                    filtered_counts <- merge(filtered_counts, tmp_counts_filtered, by = "sequence", all = TRUE)
                    tmp_cols <- vector()
                    for (i in 1:(ncol(filtered_counts) - 1)) {
                        tmp_cols <- c(tmp_cols, paste0("fil", i))
                    }
                    colnames(filtered_counts) <- c("sequence", tmp_cols)
                }
            }
            filtered_counts <- as.data.frame(filtered_counts)
            rownames(filtered_counts) <- filtered_counts$sequence
            filtered_counts <- subset(filtered_counts, select = -sequence)
            colnames(filtered_counts) <- sample_names

            tmp_seqs <- data.frame()
            for (s in object@seq_clusters) {
                if (nrow(tmp_seqs) == 0) {
                    tmp_seqs <- rownames(s[s$cluster == 1, ])
                } else {
                    tmp_seqs <- cbind_fill(tmp_seqs, rownames(s[s$cluster == 1, ]))
                }
            }
            colnames(tmp_seqs) <- sample_names
            object@bad_seqs$by_cluster <- tmp_seqs
        }

        #-------------------------------------#
        # 3. Filtering by depth and samples   #
        #    a) count >= X                    #
        #    b) in >= X% of samples           #
        #-------------------------------------#
        cat("Filtering by depth and samples...", "\n", sep = "")

        filtered_counts_final <- cbind(filtered_counts, rowSums(filtered_counts >= effect_count, na.rm = TRUE))
        filtered_counts_final <- cbind(filtered_counts_final, filtered_counts_final[, ncol(filtered_counts_final)] / length(sample_names))
        colnames(filtered_counts_final) <- c(sample_names, "sample_number", "sample_percentage")

        object@filtered_counts <- filtered_counts_final[filtered_counts_final$sample_percentage >= effect_per, sample_names]

        tmp_seqs <- rownames(filtered_counts)
        object@bad_seqs$by_count <- tmp_seqs[tmp_seqs %nin% rownames(object@filtered_counts)]

        #--------------------------------------#
        # 4. Filtering by effective mapping    #
        #    a) reads mapped to VaLiAnT output #
        #--------------------------------------#
        cat("Filtering by effective mapping...", "\n", sep = "")

        ref_seqs <- vector()
        pam_seqs <- vector()
        meta_mseqs <- vector()
        for (s in object@samples) {
            ref_seqs <- c(ref_seqs, s@refseq)
            pam_seqs <- c(pam_seqs, s@pamseq)
            meta_mseqs <- c(meta_mseqs, s@meta_mseqs)
        }
        ref_seqs <- unique(ref_seqs)
        pam_seqs <- unique(pam_seqs)
        meta_mseqs <- unique(meta_mseqs)

        # mapped to valiant but not ref/pam
        effective_counts <- object@filtered_counts[rownames(object@filtered_counts)%in%meta_mseqs, ]
        effective_counts <- effective_counts[rownames(effective_counts)%nin%c(ref_seqs, pam_seqs), ]
        object@effective_counts <- effective_counts

        # not mapped to valiant but not ref/pam
        unmapped_counts <- object@filtered_counts[rownames(object@filtered_counts)%nin%meta_mseqs, ]
        unmapped_counts <- unmapped_counts[rownames(unmapped_counts)%nin%c(ref_seqs, pam_seqs), ]

        object@bad_seqs$by_effect <- rownames(unmapped_counts)

        for (s in object@samples) {
            samplename <- s@sample
            object@stats[samplename, ]$filtered_reads <- sum(object@filtered_counts[, samplename], na.rm = TRUE)
            object@stats[samplename, ]$effective_reads <- sum(object@effective_counts[, samplename], na.rm = TRUE)
            object@stats[samplename, ]$unmapped_reads <- sum(unmapped_counts[, samplename], na.rm = TRUE)
            object@stats[samplename, ]$failed_reads <- object@stats[samplename, ]$total_reads - object@stats[samplename, ]$filtered_reads

            object@stats[samplename, ]$per_effective_reads <- object@stats[samplename, ]$effective_reads / object@stats[samplename, ]$filtered_reads
            object@stats[samplename, ]$per_unmapped_reads <- object@stats[samplename, ]$unmapped_reads / object@stats[samplename, ]$filtered_reads
            object@stats[samplename, ]$per_ref_reads <- object@stats[samplename, ]$ref_reads / object@stats[samplename, ]$filtered_reads
            object@stats[samplename, ]$per_pam_reads <- object@stats[samplename, ]$pam_reads / object@stats[samplename, ]$filtered_reads

            object@stats[samplename, ]$missing_meta_seqs <- length(s@missing_meta_seqs)
        }

        #----------------------------------------#
        # 5. Filtering by effective coverage     #
        #    a) effective reads / oligos in meta #
        #----------------------------------------#
        cat("Filtering by effective coverage...", "\n", sep = "")

        for (s in object@samples) {
            samplename <- s@sample
            object@stats[samplename, ]$effective_cov <- object@stats[samplename, ]$effective_reads / length(s@meta_mseqs)
        }
        object@stats$effective_cov <- as.integer(object@stats$effective_cov)

        #--------------------- -------------------#
        # 6. Sorting effective counts by position #
        #-----------------------------------------#
        cat("Sorting effective counts by position...", "\n", sep = "")

        effcounts_pos <- data.frame()
        for (s in object@samples) {
            obj_effcounts <- object@effective_counts[, s@sample, drop = FALSE]
            obj_effcounts <- obj_effcounts[!is.na(obj_effcounts[, 1]), 1, drop = FALSE]

            tmp_effcounts <- data.frame(matrix(NA, length(s@meta_mseqs), 1))
            colnames(tmp_effcounts) <- s@sample
            rownames(tmp_effcounts) <- s@meta_mseqs
            tmp_effcounts[, 1] <- 0

            tmp_effcounts[rownames(obj_effcounts), 1] <- obj_effcounts[, 1]
            rownames(tmp_effcounts) <- 1:length(s@meta_mseqs)

            if (nrow(effcounts_pos) == 0) {
                effcounts_pos <- tmp_effcounts
            } else {
                effcounts_pos <- cbind_fill(effcounts_pos, tmp_effcounts)
            }
        }
        # cbind_fill ignores rownames
        # assuming all the samples have the same meta sheet
        rownames(effcounts_pos) <- object@samples[[1]]@meta_mseqs
        colnames(effcounts_pos) <- sample_names

        object@effective_counts_pos <- as.data.frame(effcounts_pos)

        #------------------------#
        # 7. Gini coeff after qc #
        #------------------------#
        cat("Calculating gini coefficiency...", "\n", sep = "")

        for (s in object@samples) {
            gini_coeff <- cal_gini(object@effective_counts[, s@sample], corr = FALSE, na.rm = TRUE)
            object@stats[s@sample, ]$gini_coeff_after_qc <- round(gini_coeff, 3)
        }

        #------------------#
        # 8. QC results    #
        #------------------#
        object@stats$qcpass_filtered_reads <- unlist(lapply(object@stats$filtered_reads, function(x) ifelse(x >= cutoff_filtered, TRUE, FALSE)))
        object@stats$qcpass_mapping_per <- unlist(lapply(object@stats$per_unmapped_reads, function(x) ifelse(x < (1 - cutoff_mapping), TRUE, FALSE)))
        object@stats$qcpass_effective_per <- unlist(lapply(object@stats$per_effective_reads, function(x) ifelse(x >= cutoff_effect, TRUE, FALSE)))
        object@stats$qcpass_effective_cov <- unlist(lapply(object@stats$effective_cov, function(x) ifelse(x >= cutoff_cov, TRUE, FALSE)))

        qc_lables <- c("qcpass_filtered_reads", "qcpass_mapping_per", "qcpass_effective_per", "qcpass_effective_cov")
        object@stats$qcpass <- apply(object@stats[, qc_lables], 1, function(x) all(x))

        #------------------------#
        # 9. Filtered samples    #
        #------------------------#

        object@filtered_samples <- rownames(object@stats[object@stats$qcpass == TRUE, ])

        return(object)
    }
)

#' initialize function
setGeneric("run_sample_qc_deseq2", function(object, ...) {
  standardGeneric("run_sample_qc_deseq2")
})

#' run DESeq2 for the list of samples
#'
#' @export
#' @param object     sampleQC object
#' @param ds_coldata a data frame of col data for deseq2
#' @param ds_ref     the reference condition, eg. D4
#' @return object
setMethod(
    "run_sample_qc_deseq2",
    signature = "sampleQC",
    definition = function(object,
                          ds_coldata,
                          ds_ref) {
        #----------#
        # checking #
        #----------#
        if (length(object@filtered_samples) == 0) {
            stop(paste0("====> Error: please run run_sample_qc first!"))
        }

        if (nrow(object@samples[[1]]@vep_anno) == 0) {
            stop(paste0("====> Error: no consequence annotation found!"))
        }

        if ("condition" %nin% colnames(ds_coldata)) {
            stop(paste0("====> Error: coldata must have condition values!"))
        } else {
            ds_coldata <- as.data.frame(ds_coldata)

            ds_coldata$condition <- factor(ds_coldata$condition)
            ds_coldata$condition <- factor(ds_coldata$condition, levels = mixsort(levels(ds_coldata$condition)))

            ds_coldata$replicate <- factor(ds_coldata$replicate)
            ds_coldata$replicate <- factor(ds_coldata$replicate, levels = mixsort(levels(ds_coldata$replicate)))
        }

        #------------------------#
        # 1. map vep consequence #
        #------------------------#
        cat("Mapping consequencing annotation...", "\n", sep = "")

        # assuming all the samples have the same vep annotations
        vep_anno <- object@samples[[1]]@vep_anno
        rownames(vep_anno) <- vep_anno$Seq # may change

        # effective counts
        effective_counts_anno <- object@effective_counts[, object@filtered_samples]
        effective_counts_anno[is.na(effective_counts_anno)] <- 0

        effective_counts_anno$consequence <- vep_anno[rownames(effective_counts_anno), ]$assigned # may change
        rownames(effective_counts_anno) <- vep_anno[rownames(effective_counts_anno), ]$unique_oligo_name # may change

        object@effective_counts_anno <- effective_counts_anno

        # sorted effective counts
        effective_counts_pos_anno <- object@effective_counts_pos[, object@filtered_samples]
        effective_counts_pos_anno[is.na(effective_counts_pos_anno)] <- 0

        effective_counts_pos_anno$consequence <- vep_anno[rownames(effective_counts_pos_anno), ]$assigned # may change
        # careful: not all the meta mseqs have consequence annotation, so they don't have IDs in vep
        #rownames(effective_counts_pos_anno) <- vep_anno[rownames(effective_counts_pos_anno), ]$unique_oligo_name # may change

        object@effective_counts_pos_anno <- effective_counts_pos_anno

        #-------------------#
        # 2. DESeq2 process #
        #-------------------#

        # run control
        cat("Running control deseq2 to get size factor...", "\n", sep = "")

        syn_counts <- effective_counts_anno[effective_counts_anno$consequence == "synonymous", rownames(ds_coldata)]

        syn_ds_obj <- DESeqDataSetFromMatrix(countData = syn_counts, colData = ds_coldata, design = ~condition)
        syn_ds_obj <- syn_ds_obj[rowSums(counts(syn_ds_obj)) > 0, ]
        syn_ds_obj$condition <- factor(syn_ds_obj$condition, levels = mixsort(levels(syn_ds_obj$condition)))
        syn_ds_obj$condition <- relevel(syn_ds_obj$condition, ref = ds_ref)
        syn_ds_obj <- estimateSizeFactors(syn_ds_obj)

        # run all
        cat("Running deseq2 on all the filtered samples...", "\n", sep = "")
        deseq_counts <- effective_counts_anno[, rownames(ds_coldata)]

        ds_obj <- DESeqDataSetFromMatrix(countData = deseq_counts, colData = ds_coldata, design = ~condition)
        ds_obj <- ds_obj[rowSums(counts(ds_obj)) > 0, ]
        ds_obj$condition <- factor(ds_obj$condition, levels = mixsort(levels(ds_obj$condition)))
        ds_obj$condition <- relevel(ds_obj$condition, ref = ds_ref)
        sizeFactors(ds_obj) <- sizeFactors(syn_ds_obj)

        ds_obj <- DESeq(ds_obj, fitType = "local", quiet = TRUE)
        ds_rlog <- rlog(ds_obj)

        object@deseq_rlog <- as.data.frame(assay(ds_rlog))

        cat("Creating comparision and calculating logFC...", "\n", sep = "")
        conds <- levels(ds_coldata$condition)
        ds_contrast <- list()
        for (i in 1:length(conds)) {
            if (conds[i] != ds_ref) {
                ds_contrast <- append(ds_contrast, paste0("condition_", conds[i], "_vs_", ds_ref))
            }
        }

        ds_res <- degComps(ds_obj,
                           combs = "condition",
                           contrast = ds_contrast,
                           alpha = 0.05,
                           skip = FALSE,
                           type = "apeglm",
                           pairs = FALSE,
                           fdr = "default")

        object@deseq_res <- ds_res

        return(object)
    }
)