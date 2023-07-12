#' create QC reports
#'
#' @export
#' @name create_qc_reports
#' @param samplesheet  the path of sample sheet file
#' @param qctype       screen or plasmid
#' @param qcdir        the directory of QC plots and outs
create_qc_reports <- function(samplesheet = NULL,
                              qctype = "screen",
                              qcdir = NULL) {
        #----------#
        # checking #
        #----------#
        if (is.null(samplesheet)) {
            stop(paste0("====> Error: please provide the path of sample sheet file!"))
        }

        if (is.null(qcdir)) {
            stop(paste0("====> Error: qcdir is not provided, no output directory."))
        }

        #------------------#
        # creating reports #
        #------------------#
        report_path <- paste0(qcdir, "/", "MAVEQC_report.Rmd")

        sink(report_path)
        cat("---", "\n", sep = "")
        cat("title: \"MAVE QC Report\"", "\n", sep = "")
        cat("date: \"`r Sys.Date()`\"", "\n", sep = "")
        cat("output: html_document", "\n", sep = "")
        cat("---", "\n", sep = "")
        cat("\n", sep = "")

        cat("```{r setup, include = FALSE}", "\n", sep = "")
        cat("knitr::opts_chunk$set(echo = TRUE)", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("\n", sep = "")

        cat("---", "\n", sep = "")
        cat("\n", sep = "")

        cat("## Introduction", "\n", sep = "")
        cat("A R package of MAVE QC", "\n", sep = "")
        cat("\n", sep = "")

        cat("---", "\n", sep = "")
        cat("\n", sep = "")

        cat("## Analysis", "\n", sep = "")
        cat("\n", sep = "")

        cat("---", "\n", sep = "")
        cat("\n", sep = "")

        if (qctype == "screen") {
            cat("### Screen QC", "\n", sep = "")
        } else {
            cat("### Plasmid QC", "\n", sep = "")
        }
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        cat("#### Sample Sheet", "\n", sep = "")
        cat("\n", sep = "")
        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("library(reactable)", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", samplesheet, "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.75rem\")))", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        cat("#### Run Sample QC", "\n", sep = "")
        cat("\n", sep = "")
        cat("* ##### Read Length Distrubtion", "\n", sep = "")
        cat("\n", sep = "")
        cat("![](", qcdir, "/sample_qc_read_length.png)", "\n", sep = "")
        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", qcdir, "/sampleqc_read_length.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.75rem\")))", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        cat("* ##### Total Reads", "\n", sep = "")
        cat("\n", sep = "")
        cat("![](", qcdir, "/sample_qc_stats_total.png)", "\n", sep = "")
        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", qcdir, "/sampleqc_stats_total.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.75rem\")))", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        cat("* ##### Accepted Reads", "\n", sep = "")
        cat("\n", sep = "")
        cat("![](", qcdir, "/sample_qc_stats_accepted.png)", "\n", sep = "")
        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", qcdir, "/sampleqc_stats_accepted.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.75rem\")))", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", qcdir, "/sampleqc_stats_coverage.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.75rem\")))", "\n", sep = "")       
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        cat("* ##### Genomic Coverage", "\n", sep = "")
        cat("\n", sep = "")
        cat("![](", qcdir, "/sample_qc_position_cov.dots.png)", "\n", sep = "")
        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", qcdir, "/sampleqc_stats_pos_coverage.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.75rem\")))", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        if (qctype == "screen") {
            cat("* ##### Genomic Position Percentage", "\n", sep = "")
            cat("\n", sep = "")
            cat("![](", qcdir, "/sample_qc_position_anno.lof_dots.png)", "\n", sep = "")
            cat("```{r, echo = FALSE}", "\n", sep = "")
            cat("df <- as.data.frame(read.table(\"", qcdir, "/sampleqc_stats_pos_percentage.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
            cat("reactable(df, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = FALSE,", "\n", sep = "")
            cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.75rem\")))", "\n", sep = "")
            cat("```", "\n", sep = "")
            cat("<br>", "\n", sep = "")
            cat("\n", sep = "")

            cat("#### Run Experiment QC", "\n", sep = "")
            cat("\n", sep = "")
            cat("* ##### Sample Correlations", "\n", sep = "")
            cat("\n", sep = "")
            cat("![](", qcdir, "/sample_qc_samples_corr.png)", "\n", sep = "")
            cat("<br>", "\n", sep = "")
            cat("\n", sep = "")

            cat("* ##### Sample PCA", "\n", sep = "")
            cat("\n", sep = "")
            cat("![](", qcdir, "/sample_qc_pca_samples.png)", "\n", sep = "")
            cat("<br>", "\n", sep = "")
            cat("\n", sep = "")
        }

        sink()

        rmarkdown::render(paste0(qcdir, "/MAVEQC_report.Rmd"))
}