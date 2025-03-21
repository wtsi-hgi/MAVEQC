# Import unit test library testthat
library(testthat)

test_that("create_sampleqc_object returns expected object of class 'sampleQC'", {
   
    # create a relative path to the plasmid data
    plasmid_test_data_path <- test_path("testdata", "plasmid")

    # import_sge_files function with path to plasmid data sample sheet
    sge_objs <- import_sge_files(plasmid_test_data_path, "sample_sheet.tsv")

    sampleqc_object <- create_sampleqc_object(sge_objs)

    # test expected sample vaules in sampleqc object
    expect_s4_class(sampleqc_object, "sampleQC")

    expect_equal(class(sampleqc_object@samples), "list")
    expect_equal(sampleqc_object@lengths[[2]]$length, c(236, 150))
    expect_equal(length(sampleqc_object@samples), 2)
    expect_equal(sampleqc_object@samples_meta$sample_name, c("test_hdr741", "test_hdr742"))
    expect_equal(ncol(sampleqc_object@stats), 28)
    expect_equal(sampleqc_object@counts[[1]]$sequence[1], "GAAGAAGAACTGTGACTCATCCTGAAAACCTCTTTGAGGATTGATAGCATTTCATCAGAATCGCCATTTTGATCTACGGTTTTTGATTGCTTAGATTTTGGCAATTTTTTAGGATTAGGATTATCTATAGCACTGTCAGAAGATTTACGCTTATCCTTTTTTCTCACTGGAACTTATAGTTTTTGTTTCTCCTTAACTGTTTCATTACATTCTTCATCTGAATTAGATGTTACAGGTTTAGTTTCTGTCGGTCGCCTCAAGGGTGTAGTCTT")
    expect_equal(length(sampleqc_object@lengths), 2)
    expect_equal(length(sampleqc_object@samples_ref), 2)

    expect_equal(class(sampleqc_object@samples_ref), "list")
    expect_equal(sampleqc_object@samples_ref[[1]]@libcounts$length, c(272, 272))
    expect_equal(sampleqc_object@samples_ref[[1]]@libcounts$count, c(1546465, 4145))
    
})
