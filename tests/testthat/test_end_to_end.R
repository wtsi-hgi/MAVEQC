# Import unit test library testthat
library(testthat)
library(withr)

testthat::test_that(
  "E2E Plasmid Run",
  {
    maveqc_dir <- normalizePath(file.path(getwd(), "../.."))
    withr::with_tempdir({
      input_dir <- file.path(maveqc_dir, "test", "screen")
      temp_dir <- file.path(tempdir(), "output")
      dir.create(temp_dir)
      library_type <- "plasmid"

      sge_objs <- import_sge_files(input_dir, "sample_sheet.tsv")

      # QC1
      samqc_obj <- create_sampleqc_object(sge_objs)
      samqc <- run_sample_qc(samqc_obj, library_type)

      qcplot_samqc_all(samqc, qc_type = library_type, plot_dir = temp_dir)
      qcout_samqc_all(samqc, qc_type = library_type, out_dir = temp_dir)

      # Report
      create_qc_reports(file.path(input_dir, "sample_sheet.tsv"), library_type, temp_dir)

      # Expected files
      expected_files <- c(
                          "MAVEQC_report.html",
                          "missing_variants_in_library.tsv",
                          "sample_qc_cutoffs.tsv",
                          "sample_qc_meta.tsv",
                          "sample_qc_position_cov.dots.png",
                          "sample_qc_read_length.png",
                          "sample_qc_read_length.tsv",
                          "sample_qc_results.tsv",
                          "sample_qc_stats_accepted.png",
                          "sample_qc_stats_accepted.tsv",
                          "sample_qc_stats_coverage.tsv",
                          "sample_qc_stats_missing.tsv",
                          "sample_qc_stats_pos_counts.tsv",
                          "sample_qc_stats_pos_coverage.tsv",
                          "sample_qc_stats_total.png",
                          "sample_qc_stats_total.tsv")
      expected_files <- paste(temp_dir, "/", expected_files, sep = "")

      expected_lines_count <- c(5561, 125, 2, 10, 1238, 778, 10, 10, 606, 10, 10, 10, 2771, 10, 185, 10)
      expected_file_sizes <- c(2386741, 46647, 245, 926, 194945, 116152, 1015, 1245, 100596, 897, 960, 827, 119405, 1042, 47551, 1095)

      # Expected files
      output_files <- list.files(temp_dir, full.names = TRUE)
      output_lines_count <- unname(sapply(output_files, function(f) length(readLines(f, warn = FALSE))))
      output_file_sizes <- file.info(output_files)$size
    })

    # Checks
    expect_true(all(file.exists(expected_files)))                              # files exist (& name & format)
    expect_equal(length(output_files), 16)                                     # number of files
    expect_equal(output_lines_count, expected_lines_count)                     # number of lines
    expect_equal(output_file_sizes, expected_file_sizes)                       # file size
  }
)

testthat::test_that(
  "E2E Screen Run",
  {
    maveqc_dir <- normalizePath(file.path(getwd(), "../.."))
    withr::with_tempdir({
      input_dir <- file.path(maveqc_dir, "test", "screen")
      temp_dir <- file.path(tempdir(), "output")
      dir.create(temp_dir)
      library_type <- "screen"

      sge_objs <- import_sge_files(input_dir, "sample_sheet.tsv")

      # QC1
      samqc_obj <- create_sampleqc_object(sge_objs)
      samqc <- run_sample_qc(samqc_obj, library_type)

      qcplot_samqc_all(samqc, qc_type = library_type, plot_dir = temp_dir)
      qcout_samqc_all(samqc, qc_type = library_type, out_dir = temp_dir)

      # QC2
      expqc_obj <- create_experimentqc_object(samqc)
      expqc <- run_experiment_qc(expqc_obj)

      qcplot_expqc_all(expqc, plot_dir = temp_dir)
      qcout_expqc_all(expqc, out_dir = temp_dir)

      # Report
      create_qc_reports(file.path(input_dir, "sample_sheet.tsv"), library_type, temp_dir)

      # Expected files
      expected_files <- c(
        "MAVEQC_report.html",
        "experiment_qc_corr.tsv",
        "experiment_qc_deseq_fc.condition_Day15_vs_Day4.all.tsv",
        "experiment_qc_deseq_fc.condition_Day15_vs_Day4.all_beeswarm.png",
        "experiment_qc_deseq_fc.condition_Day15_vs_Day4.all_position.png",
        "experiment_qc_deseq_fc.condition_Day15_vs_Day4.all_sum.tsv",
        "experiment_qc_deseq_fc.condition_Day7_vs_Day4.all.tsv",
        "experiment_qc_deseq_fc.condition_Day7_vs_Day4.all_beeswarm.png",
        "experiment_qc_deseq_fc.condition_Day7_vs_Day4.all_position.png",
        "experiment_qc_deseq_fc.condition_Day7_vs_Day4.all_sum.tsv",
        "experiment_qc_pca_samples.png",
        "experiment_qc_samples_corr.png",
        "experiment_qc_samples_tree.png",
        "failed_variants_by_cluster.tsv",
        "failed_variants_by_depth.tsv",
        "failed_variants_by_mapping.tsv",
        "missing_variants_in_library.tsv",
        "sample_qc_cutoffs.tsv",
        "sample_qc_meta.tsv",
        "sample_qc_position_anno.lof_dots.png",
        "sample_qc_position_cov.dots.png",
        "sample_qc_read_length.png",
        "sample_qc_read_length.tsv",
        "sample_qc_results.tsv",
        "sample_qc_stats_accepted.png",
        "sample_qc_stats_accepted.tsv",
        "sample_qc_stats_coverage.tsv",
        "sample_qc_stats_missing.tsv",
        "sample_qc_stats_pos_counts.tsv",
        "sample_qc_stats_pos_coverage.tsv",
        "sample_qc_stats_pos_percentage.tsv",
        "sample_qc_stats_total.png",
        "sample_qc_stats_total.tsv"
      )
      expected_files <- paste(temp_dir, "/", expected_files, sep = "")
      expected_lines_count <- c(5688, 10, 2078, 1029, 1253, 9, 2078, 566, 624, 9, 452, 1284, 170, 3, 1, 1, 125, 2, 10,
                                1887, 1271, 778, 10, 10, 543, 10, 10, 10, 2076, 10, 10, 185, 10)
      expected_file_sizes <- c(4062410, 869, 1039261, 154325, 181261, 420, 1019880, 93035, 98535, 415, 72745, 158956,
                               25593, 897, 39, 39, 46647, 245, 926, 271697, 193232, 116152, 1015, 1245, 98627, 895, 960,
                               827, 89213, 1033, 1196, 47554, 1098)

      # Output files
      output_files <- list.files(temp_dir, full.names = TRUE)
      output_lines_count <- unname(sapply(output_files, function(f) length(readLines(f, warn = FALSE))))
      output_file_sizes <- file.info(output_files)$size
    })

    # Checks
    expect_true(all(file.exists(expected_files)))                              # files exist (& name & format)
    expect_equal(length(output_files), 33)                                     # number of files
    expect_equal(output_lines_count, expected_lines_count)                     # number of lines (not empty)
    expect_equal(output_file_sizes, expected_file_sizes)                       # file size (not empty)
  }
)