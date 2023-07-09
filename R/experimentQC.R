
#' initialize function
setGeneric("run_experiment_qc", function(object, ...) {
  standardGeneric("run_experiment_qc")
})

#' run DESeq2 for the list of samples
#'
#' @export
#' @param object     sampleQC object
#' @return object
setMethod(
    "run_experiment_qc",
    signature = "experimentQC",
    definition = function(object) {
        #-------------------#
        # 1. DESeq2 process #
        #-------------------#

        # run control
        cat("Running control deseq2 to get size factor...", "\n", sep = "")

        library_counts_anno <- as.data.frame(object@library_counts_anno)
        ds_coldata <- object@coldata

        # rownames are necessary for DESeq2, otherwise error happens to assign values in function
        rownames(library_counts_anno) <- library_counts_anno$seq

        syn_counts <- library_counts_anno[library_counts_anno$consequence == "Synonymous_Variant", rownames(ds_coldata)]

        syn_ds_obj <- DESeqDataSetFromMatrix(countData = syn_counts, colData = ds_coldata, design = ~condition)
        syn_ds_obj <- syn_ds_obj[rowSums(counts(syn_ds_obj)) > 0, ]
        syn_ds_obj$condition <- factor(syn_ds_obj$condition, levels = mixsort(levels(syn_ds_obj$condition)))
        syn_ds_obj$condition <- relevel(syn_ds_obj$condition, ref = object@ref_condition)
        syn_ds_obj <- estimateSizeFactors(syn_ds_obj)

        # run all
        cat("Running deseq2 on all the filtered samples...", "\n", sep = "")
        deseq_counts <- library_counts_anno[, rownames(ds_coldata)]

        ds_obj <- DESeqDataSetFromMatrix(countData = deseq_counts, colData = ds_coldata, design = ~condition)
        ds_obj <- ds_obj[rowSums(counts(ds_obj)) > 0, ]
        ds_obj$condition <- factor(ds_obj$condition, levels = mixsort(levels(ds_obj$condition)))
        ds_obj$condition <- relevel(ds_obj$condition, ref = object@ref_condition)
        sizeFactors(ds_obj) <- sizeFactors(syn_ds_obj)

        ds_obj <- DESeq(ds_obj, fitType = "local", quiet = TRUE)
        ds_rlog <- rlog(ds_obj)

        object@deseq_rlog <- as.data.frame(assay(ds_rlog))

        ds_res <- degComps(ds_obj,
                           combs = "condition",
                           contrast = object@comparisons,
                           alpha = 0.05,
                           skip = FALSE,
                           type = "apeglm",
                           pairs = FALSE,
                           fdr = "default")

        object@deseq_res <- ds_res

        return(object)
    }
)