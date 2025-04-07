# Define function that runs the analyses
generate_outputs <- function(library_type, output_dir) {
  # Prepare variables
  maveqc_dir <- normalizePath(file.path(getwd(), "../.."))
  input_dir <- file.path(maveqc_dir, "test", "screen")      # Currently only screen example data provided
  # Prepare SGE objects
  sge_objs <- suppressWarnings(import_sge_files(input_dir, "sample_sheet.tsv"))

  # Run QC1
  samqc_obj <- create_sampleqc_object(sge_objs)
  samqc <- run_sample_qc(samqc_obj, library_type)

  qcplot_samqc_all(samqc, qc_type = library_type, plot_dir = output_dir)
  qcout_samqc_all(samqc, qc_type = library_type, out_dir = output_dir)

  # Run QC2 (only for Screen QC analyses)
  if (library_type == "screen") {
    expqc_obj <- create_experimentqc_object(samqc)
    expqc <- suppressWarnings(run_experiment_qc(expqc_obj))

    suppressWarnings(qcplot_expqc_all(expqc, plot_dir = output_dir))
    qcout_expqc_all(expqc, out_dir = output_dir)
  }

  # Generate Report
  create_qc_reports(file.path(input_dir, "sample_sheet.tsv"), library_type, output_dir)
}

withr::with_tempdir("output_dir", {
  # E2E test Screen QC outputs
  testthat::test_that("End to end testing (Screen QC): correct output files are generated", {
    # Run Screen QC analysis
    library_type <- "screen"
    generate_outputs(library_type, output_dir)

    # Expected output files
    expected_file_names <- c(
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
    expected_file_names <- paste(output_dir, "/", expected_file_names, sep = "")
    expected_extensions <- tools::file_ext(expected_file_names)
    expected_file_names <- expected_file_names[order(expected_extensions, expected_file_names)]

    expected_files <- data.frame(
      file_name = expected_file_names,
      line_count = c(5688, rep(NA, count(expected_extensions == "png")), 10, 9, 2078, 9, 2078, 3, 1, 1, 125, 2, 10, 10,
                     10, 10, 10, 10, 2076, 10, 10, 10),
      file_size = c(rep(NA, count(expected_extensions == "html")), 154325, 181261, 93035, 98535, 72745, 158956, 25593,
                    271697, 193232, 116152, 98627, 47554, rep(NA, count(expected_extensions == "tsv")))
    )

    # Generated files
    generated_file_names <- list.files(output_dir, full.names = TRUE)
    generated_extensions <- tools::file_ext(generated_file_names)

    generated_files <- data.frame(
      file_name = generated_file_names,
      line_count = unname(sapply(generated_file_names, function(f) length(readLines(f, warn = FALSE)))),
      file_size = file.info(generated_file_names)$size
    )
    generated_files <- generated_files[order(generated_extensions, generated_files$file_name), ]

    # Checks
    # All expected files exit
    for (i in seq_along(expected_file_names)) {
      show_failure(expect_true(!!(expected_file_names[i]) %in% generated_file_names,
                               sprintf("Missing file: %s", expected_file_names[i])))
    }

    # No unexpected files generated
    for (i in seq_along(generated_file_names)) {
      show_failure(expect_true(!!(expected_file_names[i]) %in% generated_file_names,
                               sprintf("Unexpected file generated: %s", generated_file_names[i])))
    }

    # Check number of lines (tsv, html)
    expected_line_counts <- expected_files[grepl("\\.(html|tsv)$", expected_files$file_name),
                                           c("file_name", "line_count")]
    generated_line_counts <- generated_files[grepl("\\.(html|tsv)$", generated_files$file_name),
                                             c("file_name", "line_count")]

    for (i in seq_along(expected_line_counts$file_name)) {
      show_failure(expect_true(!!(expected_line_counts$line_count[i]) %in% generated_line_counts$line_count,
                               sprintf("%s (expected) != %s (generated) for %s",
                                       expected_line_counts$line_count[i],
                                       generated_line_counts$line_count[i],
                                       expected_line_counts$file_name[i])))
    }

    # Check file size (png)
    expected_file_sizes <- expected_files[grepl("\\.png$", expected_files$file_name),
                                          c("file_name", "file_size")]
    generated_file_sizes <- generated_files[grepl("\\.png$", generated_files$file_name),
                                            c("file_name", "file_size")]

    for (i in seq_along(expected_file_sizes$file_name)) {
      show_failure(expect_true(!!(expected_file_sizes$file_size[i]) %in% generated_file_sizes$file_size,
                               sprintf("%s (expected) != %s (generated) for %s",
                                       expected_file_sizes$file_size[i],
                                       generated_file_sizes$file_size[i],
                                       expected_file_sizes$file_name[i])))
    }
  })
})
