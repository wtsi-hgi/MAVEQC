# Import unit test library testthat
library(testthat)
library(withr)

# Helper function to import SGE files
import_test_sge_files <- function() {
    plasmid_test_data_path <- test_path("testdata", "plasmid")

    # To create a test fixture for for this test need creating sample test data
    # Such as the following files with required columns and extensions.
    # - sample_sheet.tsv
    # - _meta_consequences.tsv
    # - _meta.csv
    # - _lib_counts.tsv.gz
    # - _query_counts.tsv.gz
    # Mocking all the values and objects from the above files is intricate. hence,
    # import_sge_files is used to import the sample sheet and create the sge objects.
    # sge_objs is not mocked it is the actual object created from the test data.

    import_sge_files(plasmid_test_data_path, "sample_sheet.tsv")
}

test_that("create_sampleqc_object returns expected object of class 'sampleQC'", {
    # Use the helper function to import SGE files
    sge_objs <- import_test_sge_files()
    
    # Ensure the test data is cleaned up after the test
    withr::defer(rm(sge_objs))

    # Test expected sample values in sampleqc object
    sampleqc_object <- create_sampleqc_object(sge_objs)

    expect_s4_class(sampleqc_object, "sampleQC")

    expect_equal(class(sampleqc_object@samples), "list")
    expect_equal(length(sampleqc_object@samples), 2)
    
    expect_equal(ncol(sampleqc_object@stats), 28)

    expect_equal(sampleqc_object@samples_meta$sample_name, c("test_hdr741", "test_hdr742"))

    expect_equal(sampleqc_object@counts[[1]]$sequence[1], "GAAGAAGAACTGTGACTCATCCTGAAAACCTCTTTGAGGATTGATAGCATTTCATCAGAATCGCCATTTTGATCTACGGTTTTTGATTGCTTAGATTTTGGCAATTTTTTAGGATTAGGATTATCTATAGCACTGTCAGAAGATTTACGCTTATCCTTTTTTCTCACTGGAACTTATAGTTTTTGTTTCTCCTTAACTGTTTCATTACATTCTTCATCTGAATTAGATGTTACAGGTTTAGTTTCTGTCGGTCGCCTCAAGGGTGTAGTCTT")

    
    expect_equal(class(sampleqc_object@samples_ref), "list")
    expect_equal(length(sampleqc_object@samples_ref), 2)

    expect_equal(length(sampleqc_object@lengths), 2)
    expect_equal(sampleqc_object@lengths[[2]]$length, c(236, 150))

    expect_equal(sampleqc_object@samples_ref[[1]]@libcounts$length, c(272, 272))
    expect_equal(sampleqc_object@samples_ref[[1]]@libcounts$count, c(1546465, 4145))
  
})
