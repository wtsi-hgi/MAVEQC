test_that("qcout_samqc_accepted sets threshold using qc_type", {

    # Shared set-up for plasmid and screen test
    out_dir <- withr::local_tempdir()

    sample <- methods::new(
        "SGE",
        sample = "sample1",
        libname = "test_library"
    )

    object <- methods::new(
        "sampleQC",
        cutoffs = data.frame(
            per_library_reads_screen = 0.4,
            per_library_reads_plasmid = 0.5
        ),
        samples = list(sample),
        samples_meta = data.frame(
            sample_info = "sample1_info",
            transcript_id = "ENST01",
            exon_num = "1"
        ),
        stats = data.frame( 
            per_library_reads = 0.45,
            per_ref_reads = 0.1,
            per_pam_reads = 0.1,
            per_unmapped_reads = 0.35,
            qcpass_library_per = TRUE,
            row.names = "sample1"
        )
    )

    # Screen test
    qcout_samqc_accepted(
        object,
        qc_type = "screen",
        out_dir = out_dir
    )
    screen_output <- read.delim(
        file.path(out_dir, "sample_qc_stats_accepted.tsv"),
        check.names = FALSE
    )
    expect_equal(screen_output[["Pass Threshold (%)"]], 40)

    # Plasmid test
    qcout_samqc_accepted(
        object,
        qc_type = "plasmid",
        out_dir = out_dir
    )
    plasmid_output <- read.delim(
        file.path(out_dir, "sample_qc_stats_accepted.tsv"),
        check.names = FALSE
    )
    expect_equal(plasmid_output[["Pass Threshold (%)"]], 50)
})