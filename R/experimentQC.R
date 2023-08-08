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
                          pcut = 0.05,
                          dcut = 0,
                          ecut = 0,
                          ntop = 500) {
        #----------------------------#
        # 1. calculating size factor #
        #----------------------------#
        # run control
        cat("Running control deseq2 to get size factor...", "\n", sep = "")

        library_counts_anno <- as.data.frame(object@library_counts_anno)
        library_counts_anno[is.na(library_counts_anno)] <- 0
        ds_coldata <- object@coldata

        # rownames are necessary for DESeq2, otherwise error happens to assign values in function
        rownames(library_counts_anno) <- library_counts_anno$sequence

        syn_counts <- library_counts_anno[library_counts_anno$consequence == "Synonymous_Variant", rownames(ds_coldata)]

        syn_ds_obj <- DESeqDataSetFromMatrix(countData = syn_counts, colData = ds_coldata, design = ~condition)
        syn_ds_obj <- syn_ds_obj[rowSums(counts(syn_ds_obj)) > 0, ]
        syn_ds_obj$condition <- factor(syn_ds_obj$condition, levels = mixedsort(levels(syn_ds_obj$condition)))
        syn_ds_obj$condition <- relevel(syn_ds_obj$condition, ref = object@ref_condition)
        syn_ds_obj <- estimateSizeFactors(syn_ds_obj)

        # run all
        cat("Running deseq2 on all the filtered samples...", "\n", sep = "")
        deseq_counts <- library_counts_anno[, rownames(ds_coldata)]

        ds_obj <- DESeqDataSetFromMatrix(countData = deseq_counts, colData = ds_coldata, design = ~condition)
        ds_obj <- ds_obj[rowSums(counts(ds_obj)) > 0, ]
        ds_obj$condition <- factor(ds_obj$condition, levels = mixedsort(levels(ds_obj$condition)))
        ds_obj$condition <- relevel(ds_obj$condition, ref = object@ref_condition)
        sizeFactors(ds_obj) <- sizeFactors(syn_ds_obj)

        ds_obj <- DESeq(ds_obj, fitType = "local", quiet = TRUE)
        ds_rlog <- rlog(ds_obj)

        object@deseq_rlog <- as.data.frame(assay(ds_rlog))

        #-----------------------#
        # 2. clustering and PCA #
        #-----------------------#
        sample_dist <- dist(t(object@deseq_rlog), method = "euclidean")
        sample_hclust <- hclust(d = sample_dist, method = "ward.D2")

        object@hclust_res <- sample_hclust
        object@corr_res <- cor(scale(as.matrix(object@deseq_rlog)))

        pca_input <- as.matrix(object@deseq_rlog)
        rv <- rowVars(pca_input)
        select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
        pca <- prcomp(t(pca_input[select, ]), center = TRUE, scale = TRUE)

        object@pca_res <- pca

        #-------------------#
        # 3. DESeq2 results #
        #-------------------#

        ds_res <- degComps(ds_obj,
                           combs = "condition",
                           contrast = object@comparisons,
                           alpha = 0.05,
                           skip = FALSE,
                           type = "apeglm",
                           pairs = FALSE,
                           fdr = "default")

        object@deseq_res <- ds_res

        library_counts_pos_anno <- object@library_counts_pos_anno

        comparisions <- names(object@deseq_res)
        for (i in 1:length(object@deseq_res)) {
            res <- object@deseq_res[[i]]$shrunken[, c("log2FoldChange", "padj")]
            res$sequence <- rownames(res)
            res <- as.data.table(res)

            res[library_counts_pos_anno, position := i.position, on = .(sequence)]
            res[library_counts_pos_anno, consequence := i.consequence, on = .(sequence)]

            res$stat <- "no impact"
            res[(res$padj < pcut) & (res$log2FoldChange > ecut), ]$stat <- "enriched"
            res[(res$padj < pcut) & (res$log2FoldChange < dcut), ]$stat <- "depleted"

            res$stat <- factor(res$stat, levels = c("no impact", "enriched", "depleted"))

            object@deseq_res_anno[[comparisions[i]]] <- res
        }

        return(object)
    }
)