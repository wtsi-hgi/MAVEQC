library(testthat)
library(withr)

test_that("HTML Plasmid QC report is generated", {
    withr::local_tempdir()

    # Create files tsv
    samplesheet <- tempfile(pattern = "samplesheet", tmpdir = tempdir(), fileext = ".tsv")
    writeLines(
      c(
        "sample_name\tfastq_1\tlibrary_name\tlibrary_type\tvaliant_meta\tlibrary_dependent_count\tlibrary_independent_count\tper_r1_adaptor\tper_r2_adaptor\tvep_anno\treplicate\tref_time_point\tcondition\tadapt3\tadapt5",
        "AAAA\tsample1#.fastq.gz\tAAAA\tplasmid\t#meta.csv\tlib#_#counts.tsv\tquery_counts.tsv\t95\t0\tmeta_consequences###.tsv\tR1\tDay0\tDay0\tATCG\tATCG"
      ),
      samplesheet
    )

    writeLines(
        c(
          "num_total_reads\tper_missing_variants\tseq_low_count\tseq_low_sample_per\tnum_accepted_reads\tper_mapping_reads\tper_ref_reads\tper_library_reads\tlibrary_cov\tlow_abundance_per\tlow_abundance_lib_per",
          "100\t1\t1\t1\t100\t1\t1\t1\t1\t1\t1"
        ),
    file.path(tempdir(), "sample_qc_cutoffs.tsv"))

    writeLines(
        c(
          "Group\tSample\tSample Info\tSample Exon\t% Library Reads\t% Reference Reads\t% PAM Reads\t% Unmapped Reads\tPass Threshold (%)\tPass",
          "TEST\tAAAA\tAAAA_Day0_R1\tENST01\t50\t10\t10\t30\t40\tTRUE"
        ),
    file.path(tempdir(), "sample_qc_stats_accepted.tsv"))

    writeLines(
        c(
          "sample_name\tsample_info\tgene_id\tgene_name\ttranscript_id\texon_num\ttargeton_id\tsgrna_id",
          "AAAA\tAAAA_Day0_R1\tENST01\tgene1\tENST01\t-\tAAAA\t1"
        ),
    file.path(tempdir(), "sample_qc_meta.tsv"))

    writeLines(
        c(
          "Group\tSample\tSample Info\tSample Exon\tLibrary Sequences\tMissing Library Sequences\t% Missing Library Sequences\tPass Threshold (%)\tPass",
          "TEST\tAAAA\tAAAA_Day0_R1\tENST01\t10\t1\t0.1\t1\tTRUE"
        ),
    file.path(tempdir(), "sample_qc_stats_missing.tsv"))

    writeLines(
        c(
          "id\tname\tsequence\tlength\tcount\tunique\tsample\tis_ref\tis_pam",
          "1\tENST01_chr1:1_100_clinvar\tATCG\t100\t0\t1\tsample\t0\t0"
        ),
    file.path(tempdir(), "missing_variants_in_library.tsv"))

    writeLines(
        c(
          "Group\tSample\tSample Info\tSample Exon\tAccepted Reads\t% Accepted Reads\tExcluded Reads\t% Excluded Reads\tTotal Reads\tPass Threshold\tPass",
          "TEST\tAAAA\tAAAA_Day0_R1\tENST01\t1\t75\t1\t25\t1\t1000000\tTRUE"
        ),
    file.path(tempdir(), "sample_qc_stats_total.tsv"))

    writeLines(
        c(
          "Group\tSample\tSample Info\tSample Exon\tTotal Library Reads\tTotal Library Sequences\tLibrary Coverage\tMedian Coverage\tPass Threshold\tPass",
          "TEST\tAAAA\tAAAA_Day0_R1\tENST01\t1\t500\t1000\t5000\t100\tTRUE"
        ),
    file.path(tempdir(), "sample_qc_stats_coverage.tsv"))

    writeLines(
        c(
          "Group\tSample\tSample Info\tSample Exon\tChromosome\tStrand\tGenomic Start\tGenomic End\t% Low Abundance\tLow Abundance cutoff\tPass Threshold (%)\tPass",
          "TEST\tAAAA\tAAAA_Day0_R1\tENST01\tchr1\t+\t100\t200\t0.1\t5\t30\tTRUE"
        ),
    file.path(tempdir(), "sample_qc_stats_pos_coverage.tsv"))

    writeLines(
        c(
          "AAAA",
          "1000",
          "2000",
          "3000"
        ),
    file.path(tempdir(), "sample_qc_stats_pos_counts.tsv"))

    writeLines(
        c(
          "Group\tSample\tSample Info\tSample Exon\tGini Coefficient\tTotal Reads\t% Missing Variants\tAccepted Reads\t% Mapping Reads\t% Reference Reads\t% Library Reads\tLibrary Coverage\t% R1 Adaptor\t% R2 Adaptor\tQCPass_Library_Per",
          "TEST\tAAAA\tAAAA_Day0_R1\tENST01\t0.2\tTRUE\tTRUE\tTRUE\t50\tTRUE\tTRUE\tTRUE\t95\t0\t55"
        ),
    file.path(tempdir(), "sample_qc_results.tsv"))

    writeLines(
        c(
          "Group\tSample\tSample Info\tSample Exon\tTotal Reads\t% 0 ~ 50\t% 50 ~ 100\t% 100 ~ 150\t% 150 ~ 200\t% 200 ~ 250\t% 250 ~ 300\tPass Threshold (%)\tPass",
          "TEST\tAAAA\tAAAA_Day0_R1\tENST01\t1000\t0\t0.1\t0.5\t0.5\t2.5\t95\t90\tTRUE"
        ),
    file.path(tempdir(), "sample_qc_read_length.tsv"))

    # Create files png
    make_dummy_png <- function(name) {
        png(file.path(tempdir(), name))
        plot(1, 1)
        dev.off()
    }
    
    png_files <- c(
        "sample_qc_read_length.png",
        "sample_qc_stats_total.png",
        "sample_qc_stats_accepted.png",
        "sample_qc_position_cov.dots.png",
        "sample_qc_position_cov.dots.png"
    )
    
    lapply(png_files, make_dummy_png)

    # Run function
    create_qc_reports(
        samplesheet = samplesheet,
        qc_type = "plasmid",
        qc_dir = tempdir()
    )
  
    output_file <- file.path(tempdir(), "MAVEQC_report.html")
  
    # Check file exists and not empty
    expect_true(file.exists(output_file))
    expect_gt(file.info(output_file)$size, 0)
})


test_that("HTML Screen QC report is generated", {
    withr::local_tempdir()

    # Create files tsv
    samplesheet <- tempfile(pattern = "samplesheet", tmpdir = tempdir(), fileext = ".tsv")
    writeLines(
      c(
        "sample_name\tfastq_1\tlibrary_name\tlibrary_type\tvaliant_meta\tlibrary_dependent_count\tlibrary_independent_count\tper_r1_adaptor\tper_r2_adaptor\tvep_anno\treplicate\tref_time_point\tcondition\tadapt3\tadapt5",
        "AAAA\tsample1#.fastq.gz\tAAAA\tplasmid\t#meta.csv\tlib#_#counts.tsv\tquery_counts.tsv\t95\t0\tmeta_consequences###.tsv\tR1\tHPLE_Day4_R1\tDay4\tATCG\tATCG"
      ),
      samplesheet
    )

    writeLines(
        c(
          "num_total_reads\tper_missing_variants\tseq_low_count\tseq_low_sample_per\tnum_accepted_reads\tper_mapping_reads\tper_ref_reads\tper_library_reads\tlibrary_cov\tlow_abundance_per\tlow_abundance_lib_per",
          "100\t1\t1\t1\t100\t1\t1\t1\t1\t1\t1"
        ),
    file.path(tempdir(), "sample_qc_cutoffs.tsv"))

    writeLines(
        c(
          "Group\tSample\tSample Info\tSample Exon\t% Library Reads\t% Reference Reads\t% PAM Reads\t% Unmapped Reads\tPass Threshold (%)\tPass",
          "TEST\tAAAA\tAAAA_day4_R1\tENST01\t50\t10\t10\t30\t40\tTRUE"
        ),
    file.path(tempdir(), "sample_qc_stats_accepted.tsv"))

    writeLines(
        c(
          "sample_name\tsample_info\tgene_id\tgene_name\ttranscript_id\texon_num\ttargeton_id\tsgrna_id",
          "AAAA\tAAAA_day4_R1\tENST01\tgene1\tENST01\t-\tAAAA\t1"
        ),
    file.path(tempdir(), "sample_qc_meta.tsv"))

    writeLines(
        c(
          "Group\tSample\tSample Info\tSample Exon\tLibrary Sequences\tMissing Library Sequences\t% Missing Library Sequences\tPass Threshold (%)\tPass",
          "TEST\tAAAA\tAAAA_day4_R1\tENST01\t10\t1\t0.1\t1\tTRUE"
        ),
    file.path(tempdir(), "sample_qc_stats_missing.tsv"))

    writeLines(
        c(
          "id\tname\tsequence\tlength\tcount\tunique\tsample\tis_ref\tis_pam",
          "1\tENST01_chr1:1_100_clinvar\tATCG\t100\t0\t1\tsample\t0\t0"
        ),
    file.path(tempdir(), "missing_variants_in_library.tsv"))

    writeLines(
        c(
          "Group\tSample\tSample Info\tSample Exon\tAccepted Reads\t% Accepted Reads\tExcluded Reads\t% Excluded Reads\tTotal Reads\tPass Threshold\tPass",
          "TEST\tAAAA\tAAAA_day4_R1\tENST01\t1\t75\t1\t25\t1\t1000000\tTRUE"
        ),
    file.path(tempdir(), "sample_qc_stats_total.tsv"))

    writeLines(
        c(
          "Group\tSample\tSample Info\tSample Exon\tTotal Library Reads\tTotal Library Sequences\tLibrary Coverage\tMedian Coverage\tPass Threshold\tPass",
          "TEST\tAAAA\tAAAA_day4_R1\tENST01\t1\t500\t1000\t5000\t100\tTRUE"
        ),
    file.path(tempdir(), "sample_qc_stats_coverage.tsv"))

    writeLines(
        c(
          "Group\tSample\tSample Info\tSample Exon\tChromosome\tStrand\tGenomic Start\tGenomic End\t% Low Abundance\tLow Abundance cutoff\tPass Threshold (%)\tPass",
          "TEST\tAAAA\tAAAA_day4_R1\tENST01\tchr1\t+\t100\t200\t0.1\t5\t30\tTRUE"
        ),
    file.path(tempdir(), "sample_qc_stats_pos_coverage.tsv"))

    writeLines(
        c(
          "AAAA",
          "1000",
          "2000",
          "3000"
        ),
    file.path(tempdir(), "sample_qc_stats_pos_counts.tsv"))

    writeLines(
        c(
          "Group\tSample\tSample Info\tSample Exon\tGini Coefficient\tTotal Reads\t% Missing Variants\tAccepted Reads\t% Mapping Reads\t% Reference Reads\t% Library Reads\tLibrary Coverage\t% R1 Adaptor\t% R2 Adaptor\tQCPass_Library_Per",
          "TEST\tAAAA\tAAAA_day4_R1\tENST01\t0.2\tTRUE\tTRUE\tTRUE\t50\tTRUE\tTRUE\tTRUE\t95\t0\t55"
        ),
    file.path(tempdir(), "sample_qc_results.tsv"))

    writeLines(
        c(
          "Group\tSample\tSample Info\tSample Exon\tTotal Reads\t% 0 ~ 50\t% 50 ~ 100\t% 100 ~ 150\t% 150 ~ 200\t% 200 ~ 250\t% 250 ~ 300\tPass Threshold (%)\tPass",
          "TEST\tAAAA\tAAAA_day4_R1\tENST01\t1000\t0\t0.1\t0.5\t0.5\t2.5\t95\t90\tTRUE"
        ),
    file.path(tempdir(), "sample_qc_read_length.tsv"))

    writeLines(
        c(
          "Group\tSample\tSample Info\tSample Exon\tChromosome\tStrand\tGenomic Start\tGenomic End\t% Low Abundance (LOF)\t% Low Abundance (Others)\t% Low Abundance (ALL)\t% Low Abundance cutoff\tPass Threshold\tPass",
          "TEST\tAAAA_A_4_pel001\tAAAA_Day4_R1\tENST01\tchr1\t+\t1\t1000\t1\t2\t3\t0.005\t30\tTRUE"
        ),
    file.path(tempdir(), "sample_qc_stats_pos_percentage.tsv"))
  
    writeLines(
        c(
          "Sample\tCluster\tReplicate\tCondition\tAAAA_A_4_pel001\tPass",
          "AAAA_A_4_pel001\t1\tR1\tDay4\t0.95\tTRUE"
        ),
    file.path(tempdir(), "experiment_qc_corr.tsv"))

    writeLines(
        c(
          "oligo_name\tconsequence\tposition\tlog2FoldChange\tlfcSE\tpadj\tstat\tsequence\tadj_log2FoldChange\tadj_score\tadj_pval\tadj_fdr",
          "ENST01_chr1:1_100__CGC>TAC_aa\tMissense_Variant\t100\t-2.4601\t1\t1\tno impact\tATCG\t-2.4601\t-2.4601\t1\t1"
        ),
    file.path(tempdir(), "experiment_qc_deseq_fc.test.all.tsv"))

    writeLines(
        c(
          "consequence\tnumber of depleted\tnumber of no impact\tnumber of enriched\tlfc range\ttotal number\ttotal number of outside range",
           "Synonymous_Variant\t0\t155\t0\t-6 ~ 2\t155\t0"
          ),
    file.path(tempdir(), "experiment_qc_deseq_fc.test.all_sum.tsv"))

    # Create files png
    make_dummy_png <- function(name) {
        png(file.path(tempdir(), name))
        plot(1, 1)
        dev.off()
    }
    
    png_files <- c(
        "sample_qc_read_length.png",
        "sample_qc_stats_total.png",
        "sample_qc_stats_accepted.png",
        "sample_qc_position_cov.dots.png",
        "sample_qc_position_cov.dots.png",
        "sample_qc_position_anno.lof_dots.png",
        "experiment_qc_samples_tree.png",
        "experiment_qc_samples_corr.png",
        "experiment_qc_pca_samples.png",
        "experiment_qc_deseq_fc.test.all_beeswarm.png",
        "experiment_qc_deseq_fc.test.all_position.png"
    )
    
    lapply(png_files, make_dummy_png)

    # Run function
    create_qc_reports(
        samplesheet = samplesheet,
        qc_type = "screen",
        qc_dir = tempdir()
    )
  
    output_file <- file.path(tempdir(), "MAVEQC_report.html")
  
    # Check file exists and not empty
    expect_true(file.exists(output_file))
    expect_gt(file.info(output_file)$size, 0)
})
