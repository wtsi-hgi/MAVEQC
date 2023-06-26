#' A class representing a SGE object
#'
#' @export
#' @name SGE
#' @slot sample            the sample name
#' @slot libname           library name
#' @slot libtype           library type
#' @slot adapt5            adaptor sequence at 5 prime end
#' @slot adapt3            adaptor sequence at 3 prime end
#' @slot refseq            reference sequence
#' @slot pamseq            sequence with pam variants
#' @slot libcounts         QUANTS library-dependent counts, per sequence per count
#' @slot allcounts         QUANTS library-independent counts, per sequence per count
#' @slot valiant_meta      VaLiAnT meta file
#' @slot vep_anno          vep consequence annotation file
#' @slot meta_mseqs        non-redundant mseq in VaLiAnT meta file
#' @slot missing_meta_seqs missing sequenced in library compared to VaLiAnT meta file
#' @slot libstats          summaries of library dependent counts
#' @slot allstats          summaries of library independent counts
#' @slot libstats_qc       qc stats of library dependent counts
#' @slot allstats_qc       qc stats of library independent counts
setClass("SGE",
    slots = list(
        sample = "character",
        libname = "character",
        libtype = "character",
        adapt5 = "character",
        adapt3 = "character",
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
#' @param file_libcount           QUANTS library-dependent count file, per sequence per count
#' @param file_allcount           QUANTS library-independent count file, per sequence per count
#' @param file_valiant_meta       VaLiAnT meta file
#' @param file_vep_anno           vep annotation file
#' @param file_libcount_hline     line number of header in library-dependent count file
#' @param file_allcount_hline     line number of header in library-independent count file
#' @param file_valiant_meta_hline line number of header in VaLiAnT meta file
#' @param file_vep_anno_hline     line number of header in vep annotation file
#' @param file_libcount_cols      a vector of numbers of selected columns in library-dependent count file, default is none
#' @param file_allcount_cols      a vector of numbers of selected columns in library-independent count file, default is none
#' @param file_valiant_meta_cols  a vector of numbers of selected columns in VaLiAnT meta file, default is none
#' @param file_vep_anno_cols      a vector of numbers of selected columns in vep annotation file, default is none
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

    # Return the object
    return(sge_object)
}

#' A class representing a primary qc object
#'
#' @export
#' @name sampleQC
#' @slot samples                   a list of SGE objects
#' @slot samples_ref               a list of SGE objects which are the references for screen QC
#' @slot counts                    a list of sample counts
#' @slot lengths                   a list of sequence lengths
#' @slot seq_clusters              a list of dataframes of sequences and cluster IDs
#' @slot filtered_seqs             a vector of sequences with counts > clustering cutoff for screen QC in the references
#' @slot filtered_counts           a data frame of filtered counts of all the samples
#' @slot effective_counts          a data frame of effective counts of all the samples
#' @slot effective_counts_pos      a data frame of effective counts of all the samples sorted by position
#' @slot effective_counts_anno     a data frame of effective counts of all the samples, annotated with consequences
#' @slot effective_counts_pos_anno a data frame of effective counts of all the samples, annotated with consequences, sorted by position
#' @slot stats                     a data frame of samples and stats, eg. total no, filtered no.
#' @slot filtered_samples          a vector of filtered sample names
#' @slot bad_seqs                  a list of filter-out sequences
#' @slot deseq_rlog                a data frame of deseq rlog counts of all the samples using effective counts
#' @slot deseq_res                 a list of deseq results of all the comparison against reference
setClass("sampleQC",
    slots = list(
        samples = "list",
        samples_ref = "list",
        counts = "list",
        lengths = "list",
        seq_clusters = "list",
        filtered_seqs = "character",
        filtered_counts = "data.frame",
        effective_counts = "data.frame",
        effective_counts_pos = "data.frame",
        effective_counts_anno = "data.frame",
        effective_counts_pos_anno = "data.frame",
        stats = "data.frame",
        filtered_samples = "character",
        bad_seqs = "list",
        deseq_rlog = "data.frame",
        deseq_res = "list"
    ),
    prototype = list(
        samples = list(),
        samples_ref = list(),
        counts = list(),
        lengths = list(),
        seq_clusters = list(),
        filtered_seqs = character(),
        filtered_counts = data.frame(),
        effective_counts = data.frame(),
        effective_counts_pos = data.frame(),
        effective_counts_anno = data.frame(),
        effective_counts_pos_anno = data.frame(),
        stats = data.frame(),
        filtered_samples = character(),
        bad_seqs = list(),
        deseq_rlog = data.frame(),
        deseq_res = list()
    )
)

#' Create a new sampleQC object
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
        rownames(counts) <- counts$sequence
        counts <- subset(counts, select = -sequence)

        lengths <- s@allcounts[, "sequence", drop = FALSE]
        lengths$length <- nchar(lengths$sequence)

        list_counts[[s@sample]] <- counts
        list_lengths[[s@sample]] <- lengths
    }

    cols <- c("total_reads",
              "failed_reads",
              "filtered_reads",
              "missing_meta_seqs",
              "effective_reads",
              "per_effective_reads",
              "unmapped_reads",
              "per_unmapped_reads",
              "ref_reads",
              "per_ref_reads",
              "pam_reads",
              "per_pam_reads",
              "effective_cov",
              "gini_coeff_before_qc",
              "gini_coeff_after_qc",
              "qcpass_filtered_reads",
              "qcpass_mapping_per",
              "qcpass_effective_per",
              "qcpass_effective_cov",
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