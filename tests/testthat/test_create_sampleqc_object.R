# Import unit test library testthat
library(testthat)

test_that("Test create_sampleqc_object works correctly", {
   
    # create a relative path to the plasmid data
    plasmid_test_data_path <- test_path("testdata", "plasmid")

    # mock the import_sge_files function with mocked plasmid data
    mocked_sge_objs <- import_sge_files(plasmid_test_data_path, "sample_sheet.tsv")

    # mock the create_sampleqc_object function
    mocked_sampleqc_object <- create_sampleqc_object(mocked_sge_objs)

    # test the mocked sampleqc object
    expect_s4_class(mocked_sampleqc_object, "sampleQC", info = "Expected object to be of class sampleQC")
    expect_equal(length(mocked_sampleqc_object@samples), 1)
    expect_equal(mocked_sampleqc_object@samples_meta$sample_name, c("testsample_hdr741"), info = "Expected sample name to be testsample_hdr741")
    expect_equal(ncol(mocked_sampleqc_object@stats), 28)
    expect_equal(rownames(mocked_sampleqc_object@stats), c("testsample_hdr741"))
    expect_equal(mocked_sampleqc_object@counts[[1]]$sequence, "ATGTTACAGGTTTAGTTTCTGTCGGTCGCCTCAAGGGTGTAGTCTT")
    expect_equal(mocked_sampleqc_object@counts[[1]]$count, 1)
})
