#' initialize function
setGeneric("run_experiment_qc", function(object, ...) {
  standardGeneric("run_experiment_qc")
})

#' run DESeq2 for the list of samples
#'
#' @export
#' @param object  experimentQC object
#' @param pcut    the padj cutoff
#' @param dcut    the depleted log2 fold change cutoff
#' @param ecut    the enriched log2 fold change cutoff
#' @param ntop    the number of top variances
#' @return object
setMethod(
    "run_experiment_qc",
    signature = "experimentQC",
    definition = function(object,
                          pcut = maveqc_config$expqc_padj,
                          dcut = maveqc_config$expqc_lfc_depleted,
                          ecut = maveqc_config$expqc_lfc_enriched,
                          ntop = maveqc_config$expqc_top_variants) {
        cat("Running DESeq2 on library counts...", "\n", sep = "")
        object <- run_experiment_qc_lib_lfc(object, pcut = pcut, dcut = dcut, ecut = ecut, ntop = ntop)

        cat("Running DESeq2 on all counts...", "\n", sep = "")
        object <- run_experiment_qc_all_lfc(object, pcut = pcut, dcut = dcut, ecut = ecut)
    }
)


#' initialize function
setGeneric("run_experiment_qc_lib_lfc", function(object, ...) {
  standardGeneric("run_experiment_qc_lib_lfc")
})

#' run DESeq2 for the list of samples
#'
#' @export
#' @param object  experimentQC object
#' @param pcut    the padj cutoff
#' @param dcut    the depleted log2 fold change cutoff
#' @param ecut    the enriched log2 fold change cutoff
#' @param ntop    the number of top variances
#' @return object
setMethod(
    "run_experiment_qc_lib_lfc",
    signature = "experimentQC",
    definition = function(object,
                          pcut,
                          dcut,
                          ecut,
                          ntop) {
        #----------------------------#
        # 1. calculating size factor #
        #----------------------------#
        # run control
        cat("    |--> Running control DESeq2 to get size factor...", "\n", sep = "")

        library_counts_anno <- as.data.frame(object@library_counts_anno)
        library_counts_anno[is.na(library_counts_anno)] <- 0
        ds_coldata <- object@coldata

        # rownames are necessary for DESeq2, otherwise error happens to assign values in function
        rownames(library_counts_anno) <- library_counts_anno$sequence

        syn_counts <- library_counts_anno[library_counts_anno$consequence == "Synonymous_Variant", rownames(ds_coldata)]

        suppressMessages(syn_ds_obj <- DESeqDataSetFromMatrix(countData = syn_counts, colData = ds_coldata, design = ~condition))
        syn_ds_obj <- syn_ds_obj[rowSums(counts(syn_ds_obj)) > 0, ]
        syn_ds_obj$condition <- factor(syn_ds_obj$condition, levels = mixedsort(levels(syn_ds_obj$condition)))
        syn_ds_obj$condition <- relevel(syn_ds_obj$condition, ref = object@ref_condition)
        syn_ds_obj <- estimateSizeFactors(syn_ds_obj)

        # run all
        cat("    |--> Applying size factor to get DESeq2 normalised counts...", "\n", sep = "")
        deseq_counts <- library_counts_anno[, rownames(ds_coldata)]

        suppressMessages(ds_obj <- DESeqDataSetFromMatrix(countData = deseq_counts, colData = ds_coldata, design = ~condition))
        ds_obj <- ds_obj[rowSums(counts(ds_obj)) > 0, ]
        ds_obj$condition <- factor(ds_obj$condition, levels = mixedsort(levels(ds_obj$condition)))
        ds_obj$condition <- relevel(ds_obj$condition, ref = object@ref_condition)
        sizeFactors(ds_obj) <- sizeFactors(syn_ds_obj)

        suppressMessages(ds_obj <- DESeq(ds_obj, quiet = TRUE))
        ds_rlog <- rlog(ds_obj)

        object@lib_deseq_rlog <- as.data.frame(assay(ds_rlog))

        #-----------------------#
        # 2. clustering and PCA #
        #-----------------------#
        cat("    |--> Clustering and PCA...", "\n", sep = "")

        sample_dist <- dist(t(object@lib_deseq_rlog), method = "euclidean")
        sample_hclust <- hclust(d = sample_dist, method = "ward.D2")

        object@lib_hclust_res <- sample_hclust
        object@lib_corr_res <- cor(scale(as.matrix(object@lib_deseq_rlog)))

        pca_input <- as.matrix(object@lib_deseq_rlog)
        rv <- rowVars(pca_input)
        select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
        pca <- prcomp(t(pca_input[select, ]), center = TRUE, scale = TRUE)

        object@lib_pca_res <- pca

        #-------------------#
        # 3. DESeq2 results #
        #-------------------#
        cat("    |--> Calculating DESeq2 LFC...", "\n", sep = "")

        suppressMessages(ds_res <- degComps(ds_obj,
                                            combs = "condition",
                                            contrast = object@comparisons,
                                            alpha = 0.05,
                                            skip = FALSE,
                                            type = "apeglm",
                                            pairs = FALSE,
                                            fdr = "default"))

        object@lib_deseq_res <- ds_res

        library_counts_pos_anno <- object@library_counts_pos_anno

        comparisions <- names(object@lib_deseq_res)
        for (i in 1:length(object@lib_deseq_res)) {
            res <- object@lib_deseq_res[[i]]$shrunken[, c("log2FoldChange", "padj")]
            res$sequence <- rownames(res)
            res <- as.data.table(res)

            res[library_counts_pos_anno, c("oligo_name", "position", "consequence") := .(oligo_name, position, consequence), on = .(sequence)]

            res$stat <- "no impact"
            res[(res$padj < pcut) & (res$log2FoldChange > ecut), ]$stat <- "enriched"
            res[(res$padj < pcut) & (res$log2FoldChange < dcut), ]$stat <- "depleted"

            res$stat <- factor(res$stat, levels = c("no impact", "enriched", "depleted"))
            setcolorder(res, c("oligo_name", "consequence", "position", "log2FoldChange", "padj", "stat", "sequence"))

            res <- na.omit(res, cols = "consequence")

            object@lib_deseq_res_anno[[comparisions[i]]] <- res
        }

        return(object)
    }
)

#' initialize function
setGeneric("run_experiment_qc_all_lfc", function(object, ...) {
  standardGeneric("run_experiment_qc_all_lfc")
})

#' run DESeq2 for the list of samples
#'
#' @export
#' @param object  experimentQC object
#' @param pcut    the padj cutoff
#' @param dcut    the depleted log2 fold change cutoff
#' @param ecut    the enriched log2 fold change cutoff
#' @return object
setMethod(
    "run_experiment_qc_all_lfc",
    signature = "experimentQC",
    definition = function(object,
                          pcut,
                          dcut,
                          ecut) {
        #----------------------------#
        # 1. calculating size factor #
        #----------------------------#
        # run control
        cat("    |--> Running normalisation using total number of counts...", "\n", sep = "")

        accepted_counts <- as.data.frame(object@accepted_counts)
        accepted_counts[is.na(accepted_counts)] <- 0
        ds_coldata <- object@coldata

        # rownames are necessary for DESeq2, otherwise error happens to assign values in function
        rownames(accepted_counts) <- accepted_counts$sequence
        accepted_counts <- accepted_counts[, rownames(ds_coldata)]

        sizeFactor <- colSums(accepted_counts) / 1000000
        normFactor <- t(replicate(nrow(accepted_counts), sizeFactor))

        suppressMessages(ds_obj <- DESeqDataSetFromMatrix(countData = accepted_counts, colData = ds_coldata, design = ~condition))
        ds_obj$condition <- factor(ds_obj$condition, levels = mixedsort(levels(ds_obj$condition)))
        ds_obj$condition <- relevel(ds_obj$condition, ref = object@ref_condition)
        normalizationFactors(ds_obj) <- normFactor[, rownames(ds_coldata)]
        suppressMessages(ds_obj <- DESeq(ds_obj, quiet = TRUE))
        ds_rlog <- rlog(ds_obj)

        object@all_deseq_rlog <- as.data.frame(assay(ds_rlog))

        #---------------------------#
        # 2. calculating DESeq2 LFC #
        #---------------------------#
        cat("    |--> Calculating DESeq2 LFC...", "\n", sep = "")

        suppressMessages(ds_res <- degComps(ds_obj,
                                            combs = "condition",
                                            contrast = object@comparisons,
                                            alpha = 0.05,
                                            skip = FALSE,
                                            type = "apeglm",
                                            pairs = FALSE,
                                            fdr = "default"))

        object@all_deseq_res <- ds_res

        library_counts_pos_anno <- object@library_counts_pos_anno

        comparisions <- names(object@all_deseq_res)
        for (i in 1:length(object@all_deseq_res)) {
            res <- object@all_deseq_res[[i]]$raw[, c("log2FoldChange", "lfcSE", "padj")]
            res$sequence <- rownames(res)
            res <- as.data.table(res)

            res[library_counts_pos_anno, c("oligo_name", "position", "consequence") := .(oligo_name, position, consequence), on = .(sequence)]

            res$stat <- "no impact"
            res[(res$padj < pcut) & (res$log2FoldChange > ecut), ]$stat <- "enriched"
            res[(res$padj < pcut) & (res$log2FoldChange < dcut), ]$stat <- "depleted"

            res$stat <- factor(res$stat, levels = c("no impact", "enriched", "depleted"))
            setcolorder(res, c("oligo_name", "consequence", "position", "log2FoldChange", "lfcSE", "padj", "stat", "sequence"))

            object@all_deseq_res_anno[[comparisions[i]]] <- res
        }

        #-----------------------------------#
        # 3. adjusting DESeq2 LFC & p value #
        #-----------------------------------#
        cat("    |--> Adjusting DESeq2 LFC & p value...", "\n", sep = "")

        control_consequences <- c("Synonymous_Variant", "Intronic_Variant")
        for (i in 1:length(object@all_deseq_res_anno)) {
            control_res <- object@all_deseq_res_anno[[i]][consequence %in% control_consequences]
            res <- object@all_deseq_res_anno[[i]]

            control_median_lfc <- median(control_res$log2FoldChange)
            res$adj_log2FoldChange <- res$log2FoldChange - control_median_lfc

            res$adj_score <- res$adj_log2FoldChange / res$lfcSE
            res$adj_pval <- pnorm(abs(res$adj_score), lower.tail = FALSE) * 2
            res$adj_bh <- p.adjust(res$adj_pval, method = "BH")

            res$stat <- "no impact"
            res[(res$adj_bh < pcut) & (res$adj_log2FoldChange > ecut), ]$stat <- "enriched"
            res[(res$adj_bh < pcut) & (res$adj_log2FoldChange < dcut), ]$stat <- "depleted"

            res$stat <- factor(res$stat, levels = c("no impact", "enriched", "depleted"))

            res <- na.omit(res, cols = "consequence")

            object@all_deseq_res_anno_adj[[comparisions[i]]] <- res
        }

        return(object)
    }
)