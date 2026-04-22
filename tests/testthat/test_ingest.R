library(testthat)

test_that("import_sge_files correctly reads '#' in manifest fields", {

  # Generate files
  tmp_dir <- tempdir()

  writeLines(
    c(
      "oligo_name,species,assembly,gene_id,transcript_id,src_type,ref_chr,ref_strand,ref_start,ref_end,revc,ref_seq,pam_seq,vcf_alias,vcf_var_id,mut_position,ref,new,ref_aa,alt_aa,mut_type,mutator,oligo_length,mseq,mseq_no_adapt,pam_mut_annot,pam_mut_sgrna_id,mave_nt,mave_nt_ref,vcf_var_in_const",
      "oligo1,homo sapiens,GRCh38,GENE1,AA1,ref,chr1,+,100,200,0,ATGC,ATGC,alias1,var1,120,A,T,K,N,1del,mut1,50,ATGC,ATGC,annot1,guide1,A,A,del1"
      ),
    file.path(tmp_dir, "#meta.csv")
  )

  writeLines(
  c(
    "##Command: pyquest",
    "##Version: 1.1.0",
    "#ID\tNAME\tSEQUENCE\tLENGTH\tCOUNT\tUNIQUE\tSAMPLE",
    "1\ta\tATCG\t4\t5\t1\ts1",
    "2\tb\tATGC\t4\t10\t1\ts1"
    ),
   file.path(tmp_dir, "lib#_counts.tsv")
  )

  writeLines(
    c(
      "##Command: pyquest",
      "##Version: 1.1.0",
      "#SEQUENCE\tLENGTH\tCOUNT",
      "ATCG\t4\t5",
      "ATGC\t4\t10"
      ),
     file.path(tmp_dir, "query_counts.tsv")
    )

  writeLines(
    c(
        "sample_name\tfastq_1\tlibrary_name\tlibrary_type\tvaliant_meta\tlibrary_dependent_count\tlibrary_independent_count\tper_r1_adaptor\tper_r2_adaptor\tvep_anno\treplicate\tref_time_point\tcondition\tadapt3\tadapt5",
        "AAAA\tsample1#.fastq.gz\tAAAA\tplasmid\t#meta.csv\tlib#_counts.tsv\tquery_counts.tsv\t95\t0\tmeta_consequences.tsv\tR1\tDay0\tDay0\tATCG\tATCG"
    ),
    file.path(tmp_dir, "samplesheet.tsv")
  )

  # Run function
  sge_objects <- import_sge_files(dir_path = tmp_dir, sample_sheet = "samplesheet.tsv")

  # Check filenames are preserved
  df <- get("qc_samplesheet", envir = .GlobalEnv)

  expect_equal(df$fastq_1[1], "sample1#.fastq.gz")
  expect_equal(df$valiant_meta[1], "#meta.csv")
  expect_equal(df$library_dependent_count[1], "lib#_counts.tsv")
  expect_equal(df$library_independent_count[1], "query_counts.tsv")
  expect_equal(df$vep_anno[1], "meta_consequences.tsv")
})
