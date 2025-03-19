# Import unit test library testthat
library(testthat)

test_that("Test create_sampleqc_object works correctly", {
   
    # create a relative path to the plasmid data
    plasmid_test_data_path <- test_path("testdata", "plasmid")

    # import_sge_files function with path to plasmid data sample sheet
    sge_objs <- import_sge_files(plasmid_test_data_path, "sample_sheet.tsv")

    sampleqc_object <- create_sampleqc_object(sge_objs)

    # test expected sample vaules in sampleqc object
    expect_s4_class(sampleqc_object, "sampleQC")
    expect_equal(length(sampleqc_object@samples), 1)
    expect_equal(sampleqc_object@samples_meta$sample_name, c("testsample_hdr741"))
    expect_equal(ncol(sampleqc_object@stats), 28)
    expect_equal(rownames(sampleqc_object@stats), c("testsample_hdr741"))
    expect_equal(sampleqc_object@counts[[1]]$sequence, "ATGTTACAGGTTTAGTTTCTGTCGGTCGCCTCAAGGGTGTAGTCTT")
    expect_equal(sampleqc_object@counts[[1]]$count, 1)
})
