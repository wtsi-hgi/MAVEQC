# Define function that runs the analyses
generate_outputs <- function(library_type, input_dir, output_dir) {
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

# E2E test Screen QC outputs
test_that("End to end testing (Screen QC): correct output files are generated", {

    # Prepare variables
    library_type <- "screen"                                   # Currently only screen example data provided
    maveqc_dir <- normalizePath(file.path(getwd(), "../.."))
    input_dir <- file.path(maveqc_dir, "test", "screen")

    withr::with_tempdir({
        output_dir <- getwd()

        # Generate Screen QC outputs
        generate_outputs(library_type, input_dir, output_dir)  # Currently only screen example data provided

        # Expected output files
        expected_files <- as.data.frame(
            rbind(
                  c("MAVEQC_report.html", 5688, NA),
                  c("experiment_qc_corr.tsv", 10, 869),
                  c("experiment_qc_deseq_fc.condition_Day15_vs_Day4.all.tsv", 2078, 1039261),
                  c("experiment_qc_deseq_fc.condition_Day15_vs_Day4.all_beeswarm.png", NA, 154325),
                  c("experiment_qc_deseq_fc.condition_Day15_vs_Day4.all_position.png", NA, 181261),
                  c("experiment_qc_deseq_fc.condition_Day15_vs_Day4.all_sum.tsv", 9, 420),
                  c("experiment_qc_deseq_fc.condition_Day7_vs_Day4.all.tsv", 2078, 1019880),
                  c("experiment_qc_deseq_fc.condition_Day7_vs_Day4.all_beeswarm.png", NA, 93035),
                  c("experiment_qc_deseq_fc.condition_Day7_vs_Day4.all_position.png", NA, 98535),
                  c("experiment_qc_deseq_fc.condition_Day7_vs_Day4.all_sum.tsv", 9, 415),
                  c("experiment_qc_pca_samples.png", NA, 72745),
                  c("experiment_qc_samples_corr.png", NA, 158956),
                  c("experiment_qc_samples_tree.png", NA, 25593),
                  c("failed_variants_by_cluster.tsv", 3, 897),
                  c("failed_variants_by_depth.tsv", 1, 39),
                  c("failed_variants_by_mapping.tsv", 1, 39),
                  c("missing_variants_in_library.tsv", 125, 46647),
                  c("sample_qc_cutoffs.tsv", 2, 245),
                  c("sample_qc_meta.tsv", 10, 926),
                  c("sample_qc_position_anno.lof_dots.png", NA, 271697),
                  c("sample_qc_position_cov.dots.png", NA, 193232),
                  c("sample_qc_read_length.png", NA, 116152),
                  c("sample_qc_read_length.tsv", 10, 1015),
                  c("sample_qc_results.tsv", 10, 1245),
                  c("sample_qc_stats_accepted.png", NA, 98627),
                  c("sample_qc_stats_accepted.tsv", 10, 895),
                  c("sample_qc_stats_coverage.tsv", 10, 960),
                  c("sample_qc_stats_missing.tsv", 10, 827),
                  c("sample_qc_stats_pos_counts.tsv", 2076, 89213),
                  c("sample_qc_stats_pos_coverage.tsv", 10, 1033),
                  c("sample_qc_stats_pos_percentage.tsv", 10, 1196),
                  c("sample_qc_stats_total.png", NA, 47554),
                  c("sample_qc_stats_total.tsv", 10, 1098))
        )
        names(expected_files) <- c("file_name", "line_count", "file_size_bytes")
        expected_files$file_name <- paste(output_dir, "/", expected_files$file_name, sep = "")
        expected_files[c("line_count", "file_size_bytes")] <- lapply(expected_files[c("line_count", "file_size_bytes")],
                                                                     as.integer)

        # Generated files
        generated_file_names <- list.files(output_dir, full.names = TRUE)
        generated_files <- data.frame(
            file_name = generated_file_names,
            line_count = unname(sapply(generated_file_names, function(f) length(readLines(f, warn = FALSE)))),
            file_size_bytes = file.info(generated_file_names)$size
        )

        # Checks
        # All expected files exist
        for (i in seq_along(expected_files$file_name)) {
            expect_true(!!(expected_files$file_name[i]) %in% generated_files$file_name)
        }

        # No unexpected files generated
        for (i in seq_along(generated_files$file_name)) {
            expect_true(!!(generated_files$file_name[i]) %in% expected_files$file_name)
        }

        # Check number of lines (tsv, html)
        expected_line_counts <- expected_files[grepl("\\.(html|tsv)$", expected_files$file_name),
                                               c("file_name", "line_count")]
        generated_line_counts <- generated_files[grepl("\\.(html|tsv)$", generated_files$file_name),
                                                 c("file_name", "line_count")]

        for (i in seq_along(expected_line_counts$file_name)) {
            expect_true(!!(expected_line_counts$line_count[i]) %in% generated_line_counts$line_count,
                        label = sprintf("%s (expected) != %s (generated) for %s",
                                        expected_line_counts$line_count[i],
                                        generated_line_counts$line_count[i],
                                        expected_line_counts$file_name[i]))
        }

        # Check file size (png, tsv; html will change depending on date & time)
        expected_file_size_bytess <- expected_files[grepl("\\.(png|tsv)$", expected_files$file_name),
                                                    c("file_name", "file_size_bytes")]
        generated_file_size_bytess <- generated_files[grepl("\\.(png|tsv)$", generated_files$file_name),
                                                      c("file_name", "file_size_bytes")]

        for (i in seq_along(expected_file_size_bytess$file_name)) {
            expect_true(!!(expected_file_size_bytess$file_size_bytes[i]) %in% generated_file_size_bytess$file_size_bytes,
                        label = sprintf("%s (expected) != %s (generated) for %s",
                                        expected_file_size_bytess$file_size_bytes[i],
                                        generated_file_size_bytess$file_size_bytes[i],
                                        expected_file_size_bytess$file_name[i]))
        }
    })
})
