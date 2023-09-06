#' A class representing a SGE object
#'
#' @export
#' @name SGE
#' @slot sample             the sample name
#' @slot libname            library name
#' @slot libtype            library type
#' @slot adapt5             adaptor sequence at 5 prime end
#' @slot adapt3             adaptor sequence at 3 prime end
#' @slot per_r1_adaptor     percentage of r1 adatpor in the sample sheet
#' @slot per_r2_adaptor     percentage of r2 adatpor in the sample sheet
#' @slot refseq             reference sequence
#' @slot pamseq             sequence with pam variants
#' @slot libcounts          QUANTS library-dependent counts, per sequence per count
#' @slot allcounts          QUANTS library-independent counts, per sequence per count
#' @slot valiant_meta       VaLiAnT meta file
#' @slot vep_anno           vep consequence annotation file
#' @slot meta_mseqs         non-redundant mseq in VaLiAnT meta file
#' @slot missing_meta_seqs  missing sequenced in library compared to VaLiAnT meta file
#' @slot libstats           summaries of library dependent counts
#' @slot allstats           summaries of library independent counts
#' @slot libstats_qc        qc stats of library dependent counts
#' @slot allstats_qc        qc stats of library independent counts
setClass("SGE",
    slots = list(
        sample = "character",
        libname = "character",
        libtype = "character",
        adapt5 = "character",
        adapt3 = "character",
        per_r1_adaptor = "numeric",
        per_r2_adaptor = "numeric",
        refseq = "character",
        pamseq = "character",
        libcounts = "data.frame",
        allcounts = "data.frame",
        valiant_meta = "data.frame",
        vep_anno = "data.frame",
        meta_mseqs = "character",
        missing_meta_seqs = "character",
        libstats = "data.frame",
        allstats = "data.frame",
        libstats_qc = "data.frame",
        allstats_qc = "data.frame"
    ),
    prototype = list(
        sample = character(),
        libname = character(),
        libtype = character(),
        adapt5 = character(),
        adapt3 = character(),
        per_r1_adaptor = numeric(),
        per_r2_adaptor = numeric(),
        refseq = character(),
        pamseq = character(),
        libcounts = data.frame(),
        allcounts = data.frame(),
        valiant_meta = data.frame(),
        vep_anno = data.frame(),
        meta_mseqs = character(),
        missing_meta_seqs = character(),
        libstats = data.frame(),
        allstats = data.frame(),
        libstats_qc = data.frame(),
        allstats_qc = data.frame()
    )
)

#' Create a new SGE object
#'
#' @export
#' @name create_sge_object
#' @param file_libcount            QUANTS library-dependent count file, per sequence per count
#' @param file_allcount            QUANTS library-independent count file, per sequence per count
#' @param file_valiant_meta        VaLiAnT meta file
#' @param file_vep_anno            vep annotation file
#' @param file_libcount_hline      line number of header in library-dependent count file
#' @param file_allcount_hline      line number of header in library-independent count file
#' @param file_valiant_meta_hline  line number of header in VaLiAnT meta file
#' @param file_vep_anno_hline      line number of header in vep annotation file
#' @param file_libcount_cols       a vector of numbers of selected columns in library-dependent count file, default is none
#' @param file_allcount_cols       a vector of numbers of selected columns in library-independent count file, default is none
#' @param file_valiant_meta_cols   a vector of numbers of selected columns in VaLiAnT meta file, default is none
#' @param file_vep_anno_cols       a vector of numbers of selected columns in vep annotation file, default is none
#' @return An object of class SGE
create_sge_object <- function(file_libcount,
                              file_allcount,
                              file_valiant_meta,
                              file_vep_anno = NULL,
                              file_libcount_hline = 3,
                              file_allcount_hline = 3,
                              file_valiant_meta_hline = 1,
                              file_vep_anno_hline = 1,
                              file_libcount_cols = vector(),
                              file_allcount_cols = vector(),
                              file_valiant_meta_cols = vector(),
                              file_vep_anno_cols = vector()) {
    # Read files
    libcounts <- read_sge_file(file_libcount, file_libcount_hline, file_libcount_cols)
    allcounts <- read_sge_file(file_allcount, file_allcount_hline, file_allcount_cols)
    valiant_meta <- read_sge_file(file_valiant_meta, file_valiant_meta_hline, file_valiant_meta_cols)

    # vep is only required for screen qc
    if (is.null(file_vep_anno)) {
        vep_anno <- data.frame()
    } else {
        vep_anno <- read_sge_file(file_vep_anno, file_vep_anno_hline, file_vep_anno_cols)

        if (vep_anno[1, ]$unique_oligo_name %nin% valiant_meta$oligo_name) {
            vep_anno$unique_oligo_name <- sapply(vep_anno$unique_oligo_name, function (x) paste(head(unlist(strsplit(x, "_")), -2), collapse = "_"))
        }

        if (vep_anno[1, ]$seq %nin% libcounts$sequence) {
            if (revcomp(vep_anno[1, ]$seq) %in% libcounts$sequence) {
                vep_anno$seq <- sapply(vep_anno$seq, function (x) revcomp(x))
            }
        }
    }

    # initializing
    cols <- c("total_num_oligos",
              "total_num_unique_oligos",
              "total_counts",
              "max_counts",
              "min_counts",
              "median_counts",
              "mean_counts",
              "num_oligos_nocount",
              "num_oligos_lowcount",
              "max_len_oligos",
              "min_len_oligos")
    df_libstats <- data.frame(matrix(NA, 1, length(cols)))
    colnames(df_libstats) <- cols

    cols <- c("total_num_oligos",
              "total_num_unique_oligos",
              "total_counts",
              "max_counts",
              "min_counts",
              "median_counts",
              "mean_counts",
              "num_oligos_nocount",
              "num_oligos_lowcount",
              "max_len_oligos",
              "min_len_oligos")
    df_allstats <- data.frame(matrix(NA, 1, length(cols)))
    colnames(df_allstats) <- cols

    cols <- c("num_ref_reads",
              "per_ref_reads",
              "num_pam_reads",
              "per_pam_reads",
              "num_eff_reads",
              "per_eff_reads",
              "num_unmapped_reads",
              "per_unmapped_reads",
              "num_missing_var",
              "per_missing_var",
              "gini_coeff")
    df_libstats_qc <- data.frame(matrix(NA, 1, length(cols)))
    colnames(df_libstats_qc) <- cols

    cols <- c("num_ref_reads",
              "per_ref_reads",
              "num_pam_reads",
              "per_pam_reads",
              "num_eff_reads",
              "per_eff_reads",
              "num_unmapped_reads",
              "per_unmapped_reads",
              "num_missing_var",
              "per_missing_var",
              "gini_coeff")
    df_allstats_qc <- data.frame(matrix(NA, 1, length(cols)))
    colnames(df_allstats_qc) <- cols

    # Create the object
    sge_object <- new("SGE",
        libcounts = libcounts,
        allcounts = allcounts,
        valiant_meta = valiant_meta,
        vep_anno = vep_anno,
        libstats = df_libstats,
        allstats = df_allstats,
        libstats_qc = df_libstats_qc,
        allstats_qc = df_allstats_qc)

    sge_object@libcounts <- as.data.table(sge_object@libcounts)
    sge_object@allcounts <- as.data.table(sge_object@allcounts)
    sge_object@valiant_meta <- as.data.table(sge_object@valiant_meta)
    sge_object@vep_anno <- as.data.table(sge_object@vep_anno)

    # Return the object
    return(sge_object)
}

#' A class representing a sample QC object
#'
#' @export
#' @name sampleQC
#' @slot cutoffs                  a data frame of cutoffs using in sample QC
#' @slot samples                  a list of SGE objects
#' @slot samples_ref              a list of SGE objects which are the references for screen QC
#' @slot counts                   a list of sample libraray-independent counts
#' @slot lengths                  a list of sequence lengths
#' @slot seq_clusters             a list of dataframes of sequences and cluster IDs
#' @slot accepted_counts          a list of filtered counts of all the samples
#' @slot library_counts           a list of library counts of all the samples
#' @slot unmapped_counts          a list of unmapped counts against meta library sequences of all the samples
#' @slot library_counts_chr       a list of chromosomes for library counts
#' @slot library_counts_pos       a list of library counts of all the samples sorted by position in meta
#' @slot library_counts_anno      a data frame of library counts of all the samples, annotated with consequences
#' @slot library_counts_pos_anno  a data frame of library counts of all the samples, annotated with consequences, sorted by position
#' @slot stats                    a data frame of samples and stats, eg. total no, filtered no.
#' @slot bad_seqs_bycluster       a list of filter-out sequences by cluster
#' @slot bad_seqs_bydepth         a list of filter-out sequences by depth
#' @slot bad_seqs_bylib           a list of filter-out sequences by library mapping
#' @slot filtered_samples         a vector of filtered sample names
setClass("sampleQC",
    slots = list(
        cutoffs = "data.frame",
        samples = "list",
        samples_ref = "list",
        counts = "list",
        lengths = "list",
        seq_clusters = "list",
        accepted_counts = "list",
        library_counts = "list",
        unmapped_counts = "list",
        library_counts_chr = "list",
        library_counts_pos = "list",
        library_counts_anno = "data.frame",
        library_counts_pos_anno = "data.frame",
        stats = "data.frame",
        bad_seqs_bycluster = "list",
        bad_seqs_bydepth = "list",
        bad_seqs_bylib = "list",
        filtered_samples = "character"
    ),
    prototype = list(
        cutoffs = data.frame(),
        samples = list(),
        samples_ref = list(),
        counts = list(),
        lengths = list(),
        seq_clusters = list(),
        accepted_counts = list(),
        library_counts = list(),
        unmapped_counts = list(),
        library_counts_chr = list(),
        library_counts_pos = list(),
        library_counts_anno = data.frame(),
        library_counts_pos_anno = data.frame(),
        stats = data.frame(),
        bad_seqs_bycluster = list(),
        bad_seqs_bydepth = list(),
        bad_seqs_bylib = list(),
        filtered_samples = character()
    )
)

#' Create a new sample QC object
#'
#' @export
#' @name create_sampleqc_object
#' @param samples a list of SGE objects
#' @return An object of class sampleQC
create_sampleqc_object <- function(samples) {
    # checking
    if (length(samples) == 0) {
         stop(paste0("====> Error: no sample found in the input!"))
    }

    # initializing
    num_samples <- length(samples)
    sample_names <- character()
    for (s in samples) {
        sample_names <- append(sample_names, s@sample)
    }
    if (anyDuplicated(sample_names) != 0) {
        dup_names <- paste(unique(sample_names[duplicated(sample_names)]), collapse = ",")
        stop(paste0("====> Error: duplicated sample names:", " ", dup_names))
    }

    list_counts <- list()
    list_lengths <- list()
    for (s in samples) {
        counts <- s@allcounts[, c("sequence", "count")]

        lengths <- s@allcounts[, "sequence", drop = FALSE]
        lengths$length <- nchar(lengths$sequence) - nchar(s@adapt5) - nchar(s@adapt3)

        # in case negative length after deduction
        lengths$length[lengths$length < 0] <- 0

        list_counts[[s@sample]] <- counts
        list_lengths[[s@sample]] <- lengths
    }

    cols <- c("per_r1_adaptor",
              "per_r2_adaptor",
              "total_reads",
              "excluded_reads",
              "accepted_reads",
              "library_seqs",
              "missing_meta_seqs",
              "per_missing_meta_seqs",
              "library_reads",
              "per_library_reads",
              "unmapped_reads",
              "per_unmapped_reads",
              "ref_reads",
              "per_ref_reads",
              "pam_reads",
              "per_pam_reads",
              "median_cov",
              "library_cov",
              "gini_coeff_before_qc",
              "gini_coeff_after_qc",
              "qcpass_total_reads",
              "qcpass_missing_per",
              "qcpass_accepted_reads",
              "qcpass_mapping_per",
              "qcpass_ref_per",
              "qcpass_library_per",
              "qcpass_library_cov",
              "qcpass")
    df_stats <- data.frame(matrix(NA, num_samples, length(cols)))
    rownames(df_stats) <- sample_names
    colnames(df_stats) <- cols

    # Create the object
    sampleqc_object <- new("sampleQC",
        samples = samples,
        counts = list_counts,
        lengths = list_lengths,
        stats = df_stats)

    # Return the object
    return(sampleqc_object)
}

setClass("hclust")
setClass("prcomp")
#' A class representing a experiment QC object
#'
#' @export
#' @name experimentQC
#' @slot samples                    a list of SGE objects
#' @slot coldata                    a data frame of coldata for DESeq2
#' @slot ref_condition              the reference condition, like D4, others VS D4 in DESeq2
#' @slot vep_anno                   a data frame of consequence annotations (should be the same in all the samples for screen qc)
#' @slot accepted_counts            a data frame of accepted counts of all the samples
#' @slot library_counts_anno        a data frame of library counts of all the samples, annotated with consequences
#' @slot library_counts_pos_anno    a data frame of library counts of all the samples, annotated with consequences, sorted by position
#' @slot comparisons                a list of comparisons for degComps
#' @slot lib_deseq_rlog             a data frame of deseq rlog counts of all the samples using library counts
#' @slot lib_hclust_res             a hclust object for all the samples using library counts
#' @slot lib_corr_res               the correlation results for all the samples using library counts
#' @slot lib_pca_res                the pca results for all the samples using library counts
#' @slot lib_deseq_res              a list of deseq results of all the comparison against reference using library counts
#' @slot lib_deseq_res_anno         a list of deseq results with consequence annotations using library counts
#' @slot all_deseq_rlog             a data frame of deseq rlog counts of all the samples using all counts
#' @slot all_deseq_res              a list of deseq results of all the comparison against reference using all counts
#' @slot all_deseq_res_anno         a list of deseq results with consequence annotations using all counts
#' @slot all_deseq_res_anno_adj     a list of deseq results with consequence annotations using all counts and adjusted lfc and p value
setClass("experimentQC",
    slots = list(
        samples = "list",
        coldata = "data.frame",
        ref_condition = "character",
        vep_anno = "data.frame",
        accepted_counts = "data.frame",
        library_counts_anno = "data.frame",
        library_counts_pos_anno = "data.frame",
        comparisons = "list",
        lib_deseq_rlog = "data.frame",
        lib_hclust_res = "hclust",
        lib_corr_res = "matrix",
        lib_pca_res = "prcomp",
        lib_deseq_res = "list",
        lib_deseq_res_anno = "list",
        all_deseq_rlog = "data.frame",
        all_deseq_res = "list",
        all_deseq_res_anno = "list",
        all_deseq_res_anno_adj = "list"
    ),
    prototype = list(
        samples = list(),
        coldata = data.frame(),
        ref_condition = character(),
        vep_anno = data.frame(),
        accepted_counts = data.frame(),
        library_counts_anno = data.frame(),
        library_counts_pos_anno = data.frame(),
        comparisons = list(),
        lib_deseq_rlog = data.frame(),
        lib_hclust_res = hclust(dist(matrix(seq(1:9), nrow = 3))),
        lib_corr_res = matrix(),
        lib_pca_res = prcomp(as.data.frame(matrix(round(runif(n = 25, min = 1, max = 20), 0), nrow = 5))),
        lib_deseq_res = list(),
        lib_deseq_res_anno = list(),
        all_deseq_rlog = data.frame(),
        all_deseq_res =  list(),
        all_deseq_res_anno =  list(),
        all_deseq_res_anno_adj = list()
    )
)

#' Create a new experiment QC object
#'
#' @export
#' @name create_experimentqc_object
#' @param samqc_obj a sampleQC object
#' @param coldata   a data frame of coldata for DESeq2
#' @param refcond   the reference condition, eg. D4
#' @return An object of class sampleQC
create_experimentqc_object <- function(samqc_obj,
                                       coldata = NULL,
                                       refcond) {
    # checking
    if (is.null(coldata)) {
         stop(paste0("====> Error: no coldata found in the input!"))
    }

    if (refcond %nin% coldata$condition) {
        stop(paste0("====> Error: reference condition is not in the coldata!"))
    }

    # initializing
    if ("condition" %nin% colnames(coldata)) {
        stop(paste0("====> Error: coldata must have condition values!"))
    } else {
        coldata <- as.data.frame(coldata)

        coldata$condition <- factor(coldata$condition)
        coldata$condition <- factor(coldata$condition, levels = mixedsort(levels(coldata$condition)))

        coldata$replicate <- factor(coldata$replicate)
        coldata$replicate <- factor(coldata$replicate, levels = mixedsort(levels(coldata$replicate)))
    }

    conds <- levels(coldata$condition)
    ds_contrast <- list()
    for (i in 1:length(conds)) {
        if (conds[i] != refcond) {
            ds_contrast <- append(ds_contrast, paste0("condition_", conds[i], "_vs_", refcond))
        }
    }

    # Create the object
    experimentqc_object <- new("experimentQC",
        samples = samqc_obj@samples,
        coldata = coldata,
        ref_condition = refcond,
        vep_anno = samqc_obj@samples[[1]]@vep_anno,
        accepted_counts = merge_list_to_dt(samqc_obj@accepted_counts, "sequence", "count"),
        library_counts_anno = samqc_obj@library_counts_anno,
        library_counts_pos_anno = samqc_obj@library_counts_pos_anno,
        comparisons = ds_contrast)

    experimentqc_object@vep_anno <- as.data.table(experimentqc_object@vep_anno)
    experimentqc_object@accepted_counts <- as.data.table(experimentqc_object@accepted_counts)
    experimentqc_object@library_counts_anno <- as.data.table(experimentqc_object@library_counts_anno)
    experimentqc_object@library_counts_pos_anno <- as.data.table(experimentqc_object@library_counts_pos_anno)

    # Return the object
    return(experimentqc_object)
}
