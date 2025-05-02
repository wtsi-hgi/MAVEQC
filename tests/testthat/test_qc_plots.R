test_that("qcplot_samqc_total barplot is created with facets", {

    # Mock object
    dummy_stats <- data.frame(
        excluded_reads = seq(4, 96, 4),
        accepted_reads = rev(seq(4, 96, 4)),
        row.names = c(paste0("sample", 24:13), paste0("sample", 1:12))
    )
    mock_sample_qc <- new("sampleQC", stats = dummy_stats)

    # Generate and build the plot, and extract data
    generated_plot <- qcplot_samqc_total(mock_sample_qc, return_plot = TRUE)
    pbuilt <- ggplot2::ggplot_build(generated_plot)
    plot_data <- pbuilt$data[[1]]

    # Check that output is a ggplot object
    expect_s3_class(generated_plot, "ggplot")

    # Check the correct faceting variable is used
    facet_vars <- names(generated_plot$facet$params$facets)
    expect_true("facet_group" %in% facet_vars)

    # Check sample order
    x_labels <- levels(generated_plot$data$samples)
    expect_equal(x_labels, c(paste0("sample", 1:24), paste0("dummy", 1:16)))

    # Check number of panels
    layout <- pbuilt$layout
    expect_equal(length(unique(layout$layout$PANEL)), 2)  # 2 facet groups

    # Check 20 samples per panel
    bars_per_panel <- aggregate(x ~ PANEL, data = plot_data, FUN = function(x) length(unique(x)))
    expect_true(all(bars_per_panel$x == c(20, 4)))
  })


test_that("qcplot_samqc_accepted barplot is created with facets", {

    # Mock object
    dummy_stats <- data.frame(
        per_unmapped_reads = rep(c(.50, .25, .15, .10), 6),
        per_ref_reads = rep(c(.25, .50, .10, .15), 6),
        per_pam_reads = rep(c(.15, .10, .50, .25), 6),
        per_library_reads = rep(c(.10, .15, .25, .50), 6),
        total_reads = rep(100, 24),
        library_reads = 1:24,
        library_cov = 1:24,
        row.names = c(paste0("sample", 24:13), paste0("sample", 1:12))
    )
    mock_sample_qc <- new("sampleQC", stats = dummy_stats)

    # Generate and build the plot, and extract data
    generated_plot <- qcplot_samqc_accepted(mock_sample_qc, return_plot = TRUE)
    pbuilt <- ggplot2::ggplot_build(generated_plot)
    plot_data <- pbuilt$data[[1]]

    # Check that output is a ggplot object
    expect_s3_class(generated_plot, "ggplot")

    # Check the correct faceting variable is used
    facet_vars <- names(generated_plot$facet$params$facets)
    expect_true("facet_group" %in% facet_vars)

    # Check sample order
    x_labels <- levels(generated_plot$data$samples)
    expect_equal(x_labels, c(paste0("sample", 1:24), paste0("dummy", 1:16)))

    # Check number of panels
    layout <- pbuilt$layout
    expect_equal(length(unique(layout$layout$PANEL)), 2)

    # Check 20 samples per panel
    bars_per_panel <- aggregate(x ~ PANEL, data = plot_data, FUN = function(x) length(unique(x)))
    expect_true(all(bars_per_panel$x == c(20, 4)))
})
