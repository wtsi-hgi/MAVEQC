# Import unit test library testthat
library(testthat)

# Define function that runs the analyses
generate_outputs <- function(library_type, output_dir) {
  # Prepare variables
  maveqc_dir <- normalizePath(file.path(getwd(), "../.."))
  input_dir <- file.path(maveqc_dir, "test", "screen")     # Currently only screen example data provided
  dir.create(temp_dir)

  # Prepare SGE objects
  sge_objs <- import_sge_files(input_dir, "sample_sheet.tsv")

  # Run QC1
  samqc_obj <- create_sampleqc_object(sge_objs)
  samqc <- run_sample_qc(samqc_obj, library_type)

  qcplot_samqc_all(samqc, qc_type = library_type, plot_dir = output_dir)
  qcout_samqc_all(samqc, qc_type = library_type, out_dir = output_dir)

  # Run QC2 (only for Screen QC analyses)
  if (library_type == "screen") {
    expqc_obj <- create_experimentqc_object(samqc)
    expqc <- run_experiment_qc(expqc_obj)

    qcplot_expqc_all(expqc, plot_dir = output_dir)
    qcout_expqc_all(expqc, out_dir = output_dir)
  }

  # Generate Report
  create_qc_reports(file.path(input_dir, "sample_sheet.tsv"), library_type, output_dir)
}

# Run Plasmid QC analysis
library_type <- "plasmid"
temp_dir <- file.path(tempdir(), library_type)
generate_outputs(library_type, temp_dir)

# Expected Plasmid QC output files
expected_files_plasmid <- c(
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
  "sample_qc_stats_total.tsv"
)
expected_files_plasmid <- paste(temp_dir, "/", expected_files_plasmid, sep = "")

# Generated Plasmid QC output files
output_files_plasmid <- list.files(temp_dir, full.names = TRUE)

# Test Plasmid QC outputs
test_that("Plasmid: all expected files exist, with the right name and extension", {
  # Check
  expect_true(all(file.exists(expected_files_plasmid)),
              info = "Plasmid: expected files are missing")
})

# Not really needed given the two tests below ()
test_that("Plasmid: correct number of output files are generated (no extra files)", {
  # Check
  expect_equal(length(output_files_plasmid), length(expected_files_plasmid),
               info = "Plasmid: number of output files not as expected")
})

test_that("Plasmid: output files contain the expected number of lines", {
  # Output (same file order in expected and generated)
  expected_lines_count <- c(5561, 125, 2, 10, 1238, 778, 10, 10, 606, 10, 10, 10, 2771, 10, 185, 10)
  expected_lines_count <- setNames(expected_lines_count, expected_files_plasmid)
  expected_lines_count <- expected_lines_count[output_files_plasmid]
  output_lines_count <- sapply(output_files_plasmid, function(f) length(readLines(f, warn = FALSE)))

  # Check
  expect_equal(output_lines_count, expected_lines_count,
               info = "Plasmid: line counts not as expected")
})

test_that("Plasmid: output file sizes match expected values", {
  # Output
  expected_file_sizes <- c(2386700, 46647, 245, 926, 194945, 116152, 1015, 1245, 100596, 897, 960, 827, 119405, 1042,
                           47551, 1095)
  expected_file_sizes <- setNames(expected_file_sizes, expected_files_plasmid)
  expected_file_sizes <- expected_file_sizes[output_files_plasmid]
  output_file_sizes <- file.info(output_files_plasmid)$size
  output_file_sizes <- setNames(output_file_sizes, output_files_plasmid)

  # Check

  expect_true(all(output_file_sizes >= expected_file_sizes),
               info = "Plasmid: file sizes not as expected")
})

# Run Screen QC analysis
library_type <- "screen"
temp_dir <- file.path(tempdir(), library_type)
generate_outputs(library_type, temp_dir)

# Expected Screen QC output files
expected_files_screen <- c(
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
expected_files_screen <- paste(temp_dir, "/", expected_files_screen, sep = "")

# Generated Screen QC output files
output_files_screen <- list.files(temp_dir, full.names = TRUE)

# Test Screen QC outputs
test_that("Screen: all expected files exist, with the right name and extension", {
  # Check
  expect_true(all(file.exists(expected_files_screen)),
              info = "Screen: expected files are missing")
})

test_that("Screen: correct number of output files are generated (no extra files)", {
  # Check
  expect_equal(length(output_files_screen), length(expected_files_screen),
               info = "Screen: number of output files not as expected")
})

test_that("Screen: output files contain the expected number of lines", {
  # Output
  expected_lines_count <- c(10, 1029, 1253, 9, 2078, 566, 624, 9, 2078, 452, 1284,170, 3, 1, 1, 5688, 125, 2, 10, 1887,
                            1271, 778, 10, 10, 543, 10, 10, 10, 2076, 10, 10, 185, 10)
  expected_lines_count <- setNames(expected_lines_count, output_files_screen)
  expected_lines_count <- expected_lines_count[output_files_screen]
  output_lines_count <- sapply(output_files_screen, function(f) length(readLines(f, warn = FALSE)))

  # Check
  expect_equal(output_lines_count, expected_lines_count,
               info = "Screen: line counts not as expected")
})

test_that("Screen: output file sizes match expected values", {
  # Output
  expected_file_sizes <- c(4062400, 869, 1039261, 154325, 181261, 420, 1019880, 93035, 98535, 415, 72745, 158956,
                           25593, 897, 39, 39, 46647, 245, 926, 271697, 193232, 116152, 1015, 1245, 98627, 895, 960,
                           827, 89213, 1033, 1196, 47554, 1098)
  expected_file_sizes <- setNames(expected_file_sizes, expected_files_screen)
  expected_file_sizes <- expected_file_sizes[expected_files_screen]
  output_file_sizes <- file.info(expected_files_screen)$size
  output_file_sizes <- setNames(output_file_sizes, expected_files_screen)

  # Check
  expect_true(all(output_file_sizes >= expected_file_sizes),
              info = "Screen: file sizes not as expected")
})