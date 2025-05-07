test_that("qcplot_samqc_total barplot is created with facets (Plasmid QC)", {

    # Mock object
    test_stats <- data.frame(
        excluded_reads = seq(4, 96, 4),
        accepted_reads = rev(seq(4, 96, 4)),
        row.names = c(paste0("sample", 24:13), paste0("sample", 1:12))
    )
    mock_sample_qc <- new("sampleQC", stats = test_stats)

    # Generate and build the plot, and extract data
    generated_plot <- qcplot_samqc_total(mock_sample_qc, qc_type = "plasmid")
    pbuilt <- ggplot2::ggplot_build(generated_plot)
    plot_data <- pbuilt$data[[1]]

    # Check that output is a ggplot object
    expect_s3_class(generated_plot, "ggplot")

    # Check the correct faceting variable is used
    facet_var <- names(generated_plot$facet$params$facets)
    expect_true("facet_group" %in% facet_var)

    # Check sample order
    plot_sample_labels <- levels(generated_plot$data$plot_samples)
    expect_equal(plot_sample_labels, c(paste0("sample", 1:24), paste0("filler", 1:16)))

    # Check number of facets
    facet_n <- length(unique(pbuilt$layout$layout$PANEL))
    expect_equal(facet_n, 2)

    # Check 20 samples per facet
    bars_per_panel <- aggregate(x ~ PANEL, data = plot_data, FUN = function(x) length(unique(x)))$x
    expect_equal(bars_per_panel, c(20, 4))
  })


test_that("qcplot_samqc_total barplot is created with facets (Screen QC)", {

    # Mock object
    test_stats <- data.frame(
        excluded_reads = seq(4, 36, 4),
        accepted_reads = 100 - seq(4, 36, 4),
        row.names = c(paste0("sample", 9:5), paste0("sample", 1:4))
    )
    mock_sample_qc <- new("sampleQC", stats = test_stats)

    # Generate and build the plot, and extract data
    generated_plot <- qcplot_samqc_total(mock_sample_qc, qc_type = "screen")
    pbuilt <- ggplot2::ggplot_build(generated_plot)
    plot_data <- pbuilt$data[[1]]

    # Check that output is a ggplot object
    expect_s3_class(generated_plot, "ggplot")

    # Check the correct faceting variable is used
    facet_var <- names(generated_plot$facet$params$facets)
    expect_true("facet_group" %in% facet_var)

    # Check sample order
    plot_sample_labels <- levels(generated_plot$data$plot_samples)
    expect_equal(plot_sample_labels, paste0("sample", 1:9))

    # Check number of facets
    facet_n <- length(unique(pbuilt$layout$layout$PANEL))
    expect_equal(facet_n, 1)

    # Check 9 samples per facet
    bars_per_panel <- aggregate(x ~ PANEL, data = plot_data, FUN = function(x) length(unique(x)))$x
    expect_equal(bars_per_panel, 9)
  })


test_that("qcplot_samqc_accepted barplot is created with facets (Plasmid QC)", {

    # Mock object
    test_stats <- data.frame(
        per_unmapped_reads = rep(c(.50, .25, .15, .10), 6),
        per_ref_reads = rep(c(.25, .50, .10, .15), 6),
        per_pam_reads = rep(c(.15, .10, .50, .25), 6),
        per_library_reads = rep(c(.10, .15, .25, .50), 6),
        total_reads = rep(100, 24),
        library_reads = 1:24,
        library_cov = 1:24,
        row.names = c(paste0("sample", 24:13), paste0("sample", 1:12))
    )
    mock_sample_qc <- new("sampleQC", stats = test_stats)

    # Generate and build the plot, and extract data
    generated_plot <- qcplot_samqc_accepted(mock_sample_qc, qc_type = "plasmid")
    pbuilt <- ggplot2::ggplot_build(generated_plot)
    plot_data <- pbuilt$data[[1]]

    # Check that output is a ggplot object
    expect_s3_class(generated_plot, "ggplot")

    # Check the correct faceting variable is used
    facet_var <- names(generated_plot$facet$params$facets)
    expect_true("facet_group" %in% facet_var)

    # Check sample order
    plot_sample_labels <- levels(generated_plot$data$plot_samples)
    expect_equal(plot_sample_labels, c(paste0("sample", 1:24), paste0("filler", 1:16)))

    # Check number of facets
    facet_n <- length(unique(pbuilt$layout$layout$PANEL))
    expect_equal(facet_n, 2)

    # Check 20 samples per facet
    bars_per_panel <- aggregate(x ~ PANEL, data = plot_data, FUN = function(x) length(unique(x)))$x
    expect_equal(bars_per_panel, c(20, 4))
})


test_that("qcplot_samqc_accepted barplot is created with facets (Screen QC)", {

    # Mock object
    test_stats <- data.frame(
        per_unmapped_reads = rep(c(.50, .25, .15), 3),
        per_ref_reads = rep(c(.49, .50, .10), 3),
        per_pam_reads = rep(c(.01, .10, .50), 3),
        per_library_reads = rep(c(.00, .15, .25), 3),
        total_reads = rep(100, 9),
        library_reads = 1:9,
        library_cov = 1:9,
        row.names = c(paste0("sample", 9:5), paste0("sample", 1:4))
    )
    mock_sample_qc <- new("sampleQC", stats = test_stats)

    # Generate and build the plot, and extract data
    generated_plot <- qcplot_samqc_accepted(mock_sample_qc, qc_type = "screen", plot_dir = ".")
    pbuilt <- ggplot2::ggplot_build(generated_plot)
    plot_data <- pbuilt$data[[1]]

    # Check that output is a ggplot object
    expect_s3_class(generated_plot, "ggplot")

    # Check the correct faceting variable is used
    facet_var <- names(generated_plot$facet$params$facets)
    expect_true("facet_group" %in% facet_var)

    # Check sample order
    plot_sample_labels <- levels(generated_plot$data$plot_samples)
    expect_equal(plot_sample_labels, paste0("sample", 1:9))

    # Check number of facet
    facet_n <- length(unique(pbuilt$layout$layout$PANEL))
    expect_equal(facet_n, 1)

    # Check 20 samples per facet
    bars_per_panel <- aggregate(x ~ PANEL, data = plot_data, FUN = function(x) length(unique(x)))$x
    expect_equal(bars_per_panel, 9)
})
