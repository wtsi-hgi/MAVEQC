library(testthat)

test_that("import_sge_files correctly reads '#' in manifest fields", {
    tmp_dir <- tempdir()

    # Generate files that the samplesheet points to
    writeLines(
        c(
            "oligo_name,ref_seq,pam_seq",
            "oligo1,ATGC,ATGC"
        ),
        file.path(tmp_dir, "#meta.csv")
    )

    writeLines(
        c(
            "##Command: pyquest",
            "##Version: 1.1.0",
            "#ID\tNAME\tSEQUENCE\tLENGTH\tCOUNT\tUNIQUE\tSAMPLE",
            "1\ta\tATCG\t4\t5\t1\ts1",
            "2\tb\tATGC\t4\t10\t1\ts1"
        ),
        file.path(tmp_dir, "lib#_#counts.tsv")
      )

    writeLines(
        c(
            "##Command: pyquest",
            "##Version: 1.1.0",
            "#SEQUENCE\tLENGTH\tCOUNT",
            "ATCG\t4\t5",
            "ATGC\t4\t10"
        ),
        file.path(tmp_dir, "query_counts.tsv")
      )

    writeLines(
        c(
            "unique_oligo_name\tGene\tSYMBOL\tTargeton_ID\tsgRNA_id\tFeature\tEXON\tSeq",
            "ENST01OLG1\tgene1\tg1\tAAAA\t001\tENST01\t-\tATCG"
        ),
        file.path(tmp_dir, "meta_consequences###.tsv")
    )

    # Generate samplesheet
    writeLines(
        c(
            "sample_name\tfastq_1\tlibrary_name\tlibrary_type\tvaliant_meta\tlibrary_dependent_count\tlibrary_independent_count\tper_r1_adaptor\tper_r2_adaptor\tvep_anno\treplicate\tref_time_point\tcondition\tadapt3\tadapt5",
            "AAAA\tsample1#.fastq.gz\tAAAA\tplasmid\t#meta.csv\tlib#_#counts.tsv\tquery_counts.tsv\t95\t0\tmeta_consequences###.tsv\tR1\tDay0\tDay0\tATCG\tATCG"
        ),
        file.path(tmp_dir, "samplesheet.tsv")
    )

    # Run function and get samplesheet df object
    sge_objects <- import_sge_files(dir_path = tmp_dir,
                                    sample_sheet = "samplesheet.tsv")

    samplesheet_df <- get("qc_samplesheet", envir = .GlobalEnv)

    # Check filenames are preserved
    expect_equal(samplesheet_df$fastq_1[1], "sample1#.fastq.gz")
    expect_equal(samplesheet_df$valiant_meta[1], "#meta.csv")
    expect_equal(samplesheet_df$library_dependent_count[1], "lib#_#counts.tsv")
    expect_equal(samplesheet_df$library_independent_count[1], "query_counts.tsv")
    expect_equal(samplesheet_df$vep_anno[1], "meta_consequences###.tsv")
})
