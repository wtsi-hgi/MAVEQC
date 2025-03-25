# Import unit test library testthat
library(testthat)
library(withr)

# Helper function to import SGE files
import_test_sge_files <- function() {

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
    
    plasmid_test_data_path <- test_path("testdata", "plasmid")

    import_sge_files(plasmid_test_data_path, "sample_sheet.tsv")
}

test_that("create_sampleqc_object returns expected object of class 'sampleQC'", {

    # Define the expected column names in the stats data frame
    stats_col <- c("per_r1_adaptor", "per_r2_adaptor", "total_reads", "excluded_reads",
                   "accepted_reads", "library_seqs", "missing_meta_seqs", "per_missing_meta_seqs",
                   "library_reads", "per_library_reads", "unmapped_reads", "per_unmapped_reads",
                   "ref_reads", "per_ref_reads", "pam_reads", "per_pam_reads", "median_cov",
                   "library_cov", "gini_coeff_before_qc", "gini_coeff_after_qc", "qcpass_total_reads",
                   "qcpass_missing_per", "qcpass_accepted_reads", "qcpass_mapping_per", "qcpass_ref_per",
                   "qcpass_library_per", "qcpass_library_cov", "qcpass")

    # Use the helper function to import SGE files
    sge_objs <- import_test_sge_files()

    # Ensure the test data is cleaned up after the test
    withr::defer(rm(sge_objs))

    # Test expected sample values in sampleqc object
    sampleqc_object <- create_sampleqc_object(sge_objs)

    # Check if create_sampleqc_object returns an error when no samples are found
    expect_error(
        create_sampleqc_object(list()),
        "====> Error: no sample found in the input!"
    )

    # Check if the returned object is of class 'sampleQC'
    expect_s4_class(sampleqc_object, "sampleQC")

    expect_equal(class(sampleqc_object@samples), "list")
    expect_equal(length(sampleqc_object@samples), 2)
    
    # Check if the stats data frame has the expected number of rows and columns
    expect_equal(ncol(sampleqc_object@stats), 28)

    # Check if the column names in the stats data frame are as expected
    expect_equal(colnames(sampleqc_object@stats), stats_col)

    expect_equal(sampleqc_object@samples_meta$sample_name, c("test_hdr741", "test_hdr742"))

    expect_equal(sampleqc_object@counts[[1]]$sequence[1], "GAAGAAGAACTGTGACTCATCCTGAAAACCTCTTTGAGGATTGATAGCATTTCATCAGAATCGCCATTTTGATCTACGGTTTTTGATTGCTTAGATTTTGGCAATTTTTTAGGATTAGGATTATCTATAGCACTGTCAGAAGATTTACGCTTATCCTTTTTTCTCACTGGAACTTATAGTTTTTGTTTCTCCTTAACTGTTTCATTACATTCTTCATCTGAATTAGATGTTACAGGTTTAGTTTCTGTCGGTCGCCTCAAGGGTGTAGTCTT")

    
    expect_equal(class(sampleqc_object@samples_ref), "list")
    expect_equal(length(sampleqc_object@samples_ref), 2)

    expect_equal(length(sampleqc_object@lengths), 2)
    expect_equal(sampleqc_object@lengths[[2]]$length, c(236, 150))

    expect_equal(sampleqc_object@samples_ref[[1]]@libcounts$length, c(272, 272))
    expect_equal(sampleqc_object@samples_ref[[1]]@libcounts$count, c(1546465, 4145))
  
})
