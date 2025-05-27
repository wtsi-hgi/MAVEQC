test_that("add_filler_samples generates filler samples when needed", {
    # Test data
    test_df <- data.frame(samples = paste0("sample", 1:5))
    output <- add_filler_samples(test_df, bars_per_facet = 3)
    output_df <- output[[1]]
    output_names <- output[[2]]

    # Tests
    expect_equal(nrow(output_df), 6)
    expect_equal(colnames(output_df), c("samples", "facet_group"))
    expect_equal(output_df$facet_group, c(rep(1, 3), rep(2, 3)))
    expect_equal(output_df$samples, c(test_df$samples, "filler1"))
    expect_equal(names(output_names), "filler1")
    expect_equal(unname(output_names), "")
})

test_that("add_filler_samples does not generate filler samples when divisible", {
    # Test data
    test_df <- data.frame(samples = paste0("s", 1:6))
    output <- add_filler_samples(test_df, bars_per_facet = 2)
    output_df <- output[[1]]
    output_names <- output[[2]]

    # Tests
    expect_equal(nrow(output_df), 6)
    expect_equal(colnames(output_df), c("samples", "facet_group"))
    expect_equal(output_df$facet_group, c(rep(1, 2), rep(2, 2), rep(3, 2)))
    expect_equal(length(output_names), 0)
})
