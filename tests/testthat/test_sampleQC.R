mock_sampleqc_object <- function() {
    counts <- data.table::data.table(
        sequence = c("LIB1", "LIB2", "REF", "PAM", "NOISE"),
        # library reads % is calculated as (LIB1 + LIB2)/(LIB1 + LIB2 + REF + PAM), therefore is 0.45
        count = c(20, 25, 45, 10, 1) 
    )

    sample <- methods::new(
        "SGE",
        sample = "sample1",
        sample_meta = c(
            sample_name = "sample1",
            sample_info = "sample1_info",
            transcript_id = "ENST01",
            exon_num = "1"
        ),
        libname = "test_library",
        libtype = "screen",
        per_r1_adaptor = 0,
        per_r2_adaptor = 0,
        refseq = "REF",
        pamseq = "PAM",
        libcounts = data.frame(
            name = c("oligo1", "oligo2"),
            sequence = c("LIB1", "LIB2")
        ),
        allcounts = counts,
        valiant_meta = data.frame(
            oligo_name = c("oligo1", "oligo2"),
            mut_position = c(1, 2),
            ref_chr = "chr1",
            ref_strand = "+",
            ref_start = 1,
            ref_end = 2
        ),
        vep_anno = data.frame(
            unique_oligo_name = c("oligo1", "oligo2"),
            seq = c("LIB1", "LIB2"),
            summary_plot = c("Missense", "Synonymous")
        ),
        meta_mseqs = c("LIB1", "LIB2"),
        missing_meta_seqs = character(),
        allstats = data.frame(total_counts = sum(counts$count)),
        allstats_qc = data.frame(num_ref_reads = 45, num_pam_reads = 10),
        libstats_qc = data.frame(gini_coeff = 0.1)
    )
    sample@libcounts <- data.table::as.data.table(sample@libcounts)
    sample@allcounts <- data.table::as.data.table(sample@allcounts)
    sample@valiant_meta <- data.table::as.data.table(sample@valiant_meta)
    sample@vep_anno <- data.table::as.data.table(sample@vep_anno)

    stats_columns <- c(
        "per_r1_adaptor", "per_r2_adaptor", "total_reads", "excluded_reads",
        "accepted_reads", "library_seqs", "missing_meta_seqs",
        "per_missing_meta_seqs", "library_reads", "per_library_reads",
        "unmapped_reads", "per_unmapped_reads", "ref_reads", "per_ref_reads",
        "pam_reads", "per_pam_reads", "median_cov", "library_cov",
        "gini_coeff_before_qc", "gini_coeff_after_qc", "qcpass_total_reads",
        "qcpass_missing_per", "qcpass_accepted_reads", "qcpass_mapping_per",
        "qcpass_ref_per", "qcpass_library_per", "qcpass_library_cov", "qcpass"
    )
    stats <- data.frame(matrix(NA, nrow = 1, ncol = length(stats_columns)))
    names(stats) <- stats_columns
    rownames(stats) <- sample@sample

    methods::new(
        "sampleQC",
        samples = list(sample),
        samples_ref = list(sample),
        samples_meta = data.frame(
            sample_info = "sample1_info",
            transcript_id = "ENST01",
            exon_num = "1"
        ),
        stats = stats
    )
}

test_that("run_sample_qc uses the plasmid library read threshold", {
    result <- run_sample_qc(
        mock_sampleqc_object(),
        qc_type = "plasmid",
        cutoff_low_count = 0,
        screen_cutoff_library_per = 0.4,
        plasmid_cutoff_library_per = 0.5
    )

    expect_false(result@stats$qcpass_library_per)
    expect_equal(result@cutoffs$per_library_reads_plasmid, 0.5)
})

test_that("run_sample_qc uses the screen library read threshold", {
    result <- run_sample_qc(
        mock_sampleqc_object(),
        qc_type = "screen",
        cutoff_low_count = 0,
        screen_cutoff_library_per = 0.4,
        plasmid_cutoff_library_per = 0.5
    )

    expect_true(result@stats$qcpass_library_per)
    expect_equal(result@cutoffs$per_library_reads_screen, 0.4)
})
