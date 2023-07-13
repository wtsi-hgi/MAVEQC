#' create QC reports
#'
#' @export
#' @name create_qc_reports
#' @param samplesheet    the path of sample sheet file
#' @param qctype         screen or plasmid
#' @param qcdir          the directory of QC plots and outs
create_qc_reports <- function(samplesheet = NULL,
                              qctype = "screen",
                              qcdir = NULL) {
        #----------#
        # checking #
        #----------#
        if (is.null(samplesheet)) {
            stop(paste0("====> Error: please provide the path of sample sheet file!"))
        }

        if (qctype %nin% c("screen", "plasmid")) {
            stop(paste0("====> Error: wrong qctype, please use screen or plasmid!"))
        }

        if (is.null(qcdir)) {
            stop(paste0("====> Error: qcdir is not provided, no output directory."))
        }

        if (!file.exists(paste0(qcdir, "/sample_qc_cutoffs.tsv"))) {
            stop(paste0("====> Error: sample_qc_cutoffs.tsv is not in ", qcdir, ". Please use qcout_samqc_cutoffs to create it."))
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
        cat("knitr::opts_chunk$set(echo = TRUE, fig.align = \"center\")", "\n", sep = "")
        cat("library(reactable)", "\n", sep = "")
        cat("outdir <- \"", qcdir, "\"", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("\n", sep = "")

        cat("---", "\n", sep = "")
        cat("\n", sep = "")

        cat("## 1. Introduction", "\n", sep = "")
        cat("MAVEQC is a flexible R-package that provides QC analysis of Saturation Genome Editing (SGE) experimental data. ",
            "Available under GPL 3.0 from https://github.com/wtsi-hgi/MAVEQC", "\n", sep = "")
        cat("\n", sep = "")

        cat("---", "\n", sep = "")
        cat("\n", sep = "")

        if (qctype == "screen") {
            cat("### 2. Screen QC", "\n", sep = "")
        } else {
            cat("### 2. Plasmid QC", "\n", sep = "")
        }
        cat("Displays QC plots and statistics for all samples for QC.", "\n", sep = "")
        cat("\n", sep = "")
        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("samqc_cutoffs <- as.data.frame(read.table(\"", qcdir, "/sample_qc_cutoffs.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        cat("#### 2.1. Sample Sheet", "\n", sep = "")
        cat("\n", sep = "")
        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", samplesheet, "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")))", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        cat("#### 2.2. Run Sample QC", "\n", sep = "")
        cat("\n", sep = "")
        cat("##### 2.2.1. Read Length Distrubtion", "\n", sep = "")
        cat("Displays the percentage of reads for each sample, based on 50 nucleotide increments.", "\n", sep = "")
        cat("\n", sep = "")

        cat("```{r, echo = FALSE, out.height = \"50%\", out.width = \"50%\"}", "\n", sep = "")
        cat("knitr::include_graphics(paste0(outdir, \"/sample_qc_read_length.png\"), rel_path = FALSE)", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", qcdir, "/sample_qc_read_length.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
        cat("          columns = list(\"Group\" = colDef(minWidth = 150),", "\n", sep = "")
        cat("                         \"Sample\" = colDef(minWidth = 150),", "\n", sep = "")
        cat("                         \"Total Reads\" = colDef(format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"Pass\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" })))", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        cat("##### 2.2.2. Total Reads", "\n", sep = "")
        cat("Displays the total number of reads per sample. ",
            "Filtering based on 1-dimensional Kmean clustering that excludes unique sequences with low read counts.", "\n", sep = "")
        cat("\n", sep = "")
        cat("* **Accepted reads**: Total read count for all unique sequences with sufficient reads.", "\n", sep = "")
        cat("* **Excluded reads**: Total read count for all unique sequences with insufficient reads.", "\n", sep = "")
        cat("\n", sep = "")

        cat("```{r, echo = FALSE, out.height = \"50%\", out.width = \"50%\"}", "\n", sep = "")
        cat("knitr::include_graphics(paste0(outdir, \"/sample_qc_stats_total.png\"), rel_path = FALSE)", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", qcdir, "/sample_qc_stats_total.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
        cat("          columns = list(\"Group\" = colDef(minWidth = 150),", "\n", sep = "")
        cat("                         \"Sample\" = colDef(minWidth = 150),", "\n", sep = "")
        cat("                         \"Accepted Reads\" = colDef(format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"Excluded Reads\" = colDef(format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"Total Reads\" = colDef(format = colFormat(separators = TRUE),", "\n", sep = "")
        cat("                                                  style = function(value) { ", "\n", sep = "")
        cat("                                                              if (value < samqc_cutoffs$total_reads) {", "\n", sep = "")
        cat("                                                                  color <- \"red\"", "\n", sep = "")
        cat("                                                                  fweight <- \"bold\"", "\n", sep = "")
        cat("                                                              } else {", "\n", sep = "")
        cat("                                                                  color <- \"forestgreen\"", "\n", sep = "")
        cat("                                                                  fweight <- \"plain\" }", "\n", sep = "")
        cat("                                                              list(color = color, fontWeight = fweight)}),", "\n", sep = "")
        cat("                         \"Pass Threshold\" = colDef(format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"Pass\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" })))", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        cat("##### 2.2.3. Accepted Reads", "\n", sep = "")
        cat("Displays the percentage of library reads vs non-library reads (ie. Reference, PAM and Unmapped) for Accepted Reads.", "\n", sep = "")
        cat("\n", sep = "")
        cat("* **Library Reads**: Percentage reads mapping to template oligo sequences, including intended variants.", "\n", sep = "")
        cat("* **Reference Reads**: Percentage reads mapping to Reference.", "\n", sep = "")
        cat("* **PAM_Reads**: Percentage reads mapping to PAM/Protospacer Protection Edits (PPEs), without intended variant.", "\n", sep = "")
        cat("* **Unmapped Reads**: Percentage of Unmapped Reads.", "\n", sep = "")
        cat("\n", sep = "")

        cat("```{r, echo = FALSE, out.height = \"50%\", out.width = \"50%\"}", "\n", sep = "")
        cat("knitr::include_graphics(paste0(outdir, \"/sample_qc_stats_accepted.png\"), rel_path = FALSE)", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", qcdir, "/sample_qc_stats_accepted.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
        cat("          columns = list(\"Group\" = colDef(minWidth = 150),", "\n", sep = "")
        cat("                         \"Sample\" = colDef(minWidth = 150),", "\n", sep = "")
        cat("                         \"% Library Reads\" = colDef(style = function(value) {", "\n", sep = "")
        cat("                                                                  if (value < samqc_cutoffs$library_percent * 100) {", "\n", sep = "")
        cat("                                                                      color <- \"red\"", "\n", sep = "")
        cat("                                                                      fweight <- \"bold\"", "\n", sep = "")
        cat("                                                                  } else {", "\n", sep = "")
        cat("                                                                      color <- \"forestgreen\"", "\n", sep = "")
        cat("                                                                      fweight <- \"plain\" }", "\n", sep = "")
        cat("                                                                  list(color = color, fontWeight = fweight)}),", "\n", sep = "")
        cat("                         \"Pass\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" })))", "\n", sep = "")
        cat("```", "\n", sep = "")

        cat("\n", sep = "")
        cat("Defines the mean read count per template oligo sequence.", "\n", sep = "")
        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", qcdir, "/sample_qc_stats_coverage.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
        cat("          columns = list(\"Group\" = colDef(minWidth = 150),", "\n", sep = "")
        cat("                         \"Sample\" = colDef(minWidth = 150),", "\n", sep = "")
        cat("                         \"Total Library Reads\" = colDef(format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"Total Template Oligo Sequences\" = colDef(format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"Library Coverage\" = colDef(format = colFormat(separators = TRUE),", "\n", sep = "")
        cat("                                                       style = function(value) {", "\n", sep = "")
        cat("                                                                   if (value < samqc_cutoffs$library_cov) {", "\n", sep = "")
        cat("                                                                       color <- \"red\"", "\n", sep = "")
        cat("                                                                       fweight <- \"bold\"", "\n", sep = "")
        cat("                                                                   } else {", "\n", sep = "")
        cat("                                                                       color <- \"forestgreen\"", "\n", sep = "")
        cat("                                                                       fweight <- \"plain\" }", "\n", sep = "")
        cat("                                                                   list(color = color, fontWeight = fweight)}),", "\n", sep = "")
        cat("                         \"Pass\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" })))", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        cat("##### 2.2.4. Genomic Coverage", "\n", sep = "")
        cat("Distribution of variants across targeton region based on log2(count+1) values.", "\n", sep = "")
        cat("\n", sep = "")

        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("knitr::include_graphics(paste0(outdir, \"/sample_qc_position_cov.dots.png\"), rel_path = FALSE)", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", qcdir, "/sample_qc_stats_pos_coverage.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
        cat("          columns = list(\"Group\" = colDef(minWidth = 150),", "\n", sep = "")
        cat("                         \"Sample\" = colDef(minWidth = 150),", "\n", sep = "")
        cat("                         \"Genomic Start\" = colDef(format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"Genomic End\" = colDef(format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"% Low Abundance\" = colDef(minWidth = 150,", "\n", sep = "")
        cat("                                                      style = function(value) {", "\n", sep = "")
        cat("                                                                  if (value > (1 - samqc_cutoffs$low_abundance_lib_per) * 100) {", "\n", sep = "")
        cat("                                                                      color <- \"red\"", "\n", sep = "")
        cat("                                                                      fweight <- \"bold\"", "\n", sep = "")
        cat("                                                                  } else {", "\n", sep = "")
        cat("                                                                      color <- \"forestgreen\"", "\n", sep = "")
        cat("                                                                      fweight <- \"plain\" }", "\n", sep = "")
        cat("                                                                  list(color = color, fontWeight = fweight)}),", "\n", sep = "")
        cat("                         \"Low Abundance cutoff\" = colDef(minWidth = 150),", "\n", sep = "")
        cat("                         \"Pass\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" })))", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        if (qctype == "screen") {
            cat("##### 2.2.5. Genomic Position Percentage", "\n", sep = "")
            cat("Displays distribution of \"LOF\" (loss-of-function) vs all \"Other\" variants across the targeton region, ",
                "based on read percentages for reference timepoint. ",
                "Requires concordant distribution of LOF and Other variants.", "\n", sep = "")
            cat("\n", sep = "")

            cat("```{r, echo = FALSE, out.height = \"65%\", out.width = \"65%\"}", "\n", sep = "")
            cat("knitr::include_graphics(paste0(outdir, \"/sample_qc_position_anno.lof_dots.png\"), rel_path = FALSE)", "\n", sep = "")
            cat("```", "\n", sep = "")
            cat("```{r, echo = FALSE}", "\n", sep = "")
            cat("df <- as.data.frame(read.table(\"", qcdir, "/sample_qc_stats_pos_percentage.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
            cat("reactable(df, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
            cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
            cat("          columns = list(\"Group\" = colDef(minWidth = 150),", "\n", sep = "")
            cat("                         \"Sample\" = colDef(minWidth = 150),", "\n", sep = "")
            cat("                         \"Genomic Start\" = colDef(format = colFormat(separators = TRUE)),", "\n", sep = "")
            cat("                         \"Genomic End\" = colDef(format = colFormat(separators = TRUE)),", "\n", sep = "")
            cat("                         \"% Low Abundance (LOF)\" = colDef(minWidth = 150),", "\n", sep = "")
            cat("                         \"% Low Abundance (Others)\" = colDef(minWidth = 150),", "\n", sep = "")
            cat("                         \"% Low Abundance (ALL)\" = colDef(minWidth = 150,", "\n", sep = "")
            cat("                                                            style = function(value) {", "\n", sep = "")
            cat("                                                                        if (value > (1 - samqc_cutoffs$low_abundance_lib_per) * 100) {", "\n", sep = "")
            cat("                                                                            color <- \"red\"", "\n", sep = "")
            cat("                                                                            fweight <- \"bold\"", "\n", sep = "")
            cat("                                                                        } else {", "\n", sep = "")
            cat("                                                                            color <- \"forestgreen\"", "\n", sep = "")
            cat("                                                                            fweight <- \"plain\" }", "\n", sep = "")
            cat("                                                                        list(color = color, fontWeight = fweight)}),", "\n", sep = "")
            cat("                         \"% Low Abundance cutoff\" = colDef(minWidth = 150),", "\n", sep = "")
            cat("                         \"Pass\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" })))", "\n", sep = "")
            cat("```", "\n", sep = "")
            cat("<br>", "\n", sep = "")
            cat("\n", sep = "")

            cat("#### 2.3. Run Experiment QC", "\n", sep = "")
            cat("\n", sep = "")
            cat("##### 2.3.1. Sample Correlations", "\n", sep = "")
            cat("\n", sep = "")
            cat("```{r, echo = FALSE, out.height = \"50%\", out.width = \"50%\"}", "\n", sep = "")
            cat("knitr::include_graphics(paste0(outdir, \"/sample_qc_samples_corr.png\"), rel_path = FALSE)", "\n", sep = "")
            cat("```", "\n", sep = "")
            cat("<br>", "\n", sep = "")
            cat("\n", sep = "")

            cat("##### 2.3.2. Sample PCA", "\n", sep = "")
            cat("\n", sep = "")
            cat("```{r, echo = FALSE, out.height = \"75%\", out.width = \"75%\"}", "\n", sep = "")
            cat("knitr::include_graphics(paste0(outdir, \"/sample_qc_pca_samples.png\"), rel_path = FALSE)", "\n", sep = "")
            cat("```", "\n", sep = "")
            cat("<br>", "\n", sep = "")
            cat("\n", sep = "")

            cat("##### 2.3.3. Folder Change", "\n", sep = "")
            cat("\n", sep = "")
            cat("```{r, echo = FALSE, out.height = \"75%\", out.width = \"75%\"}", "\n", sep = "")
            cat("figs <- list.files(path = outdir, pattern = \"sample_qc_deseq_fc.*.violin.png\", full.names = TRUE)", "\n", sep = "")
            cat("knitr::include_graphics(figs, rel_path = FALSE)", "\n", sep = "")
            cat("```", "\n", sep = "")
            cat("<br>", "\n", sep = "")
            cat("\n", sep = "")
        }

        sink()

        rmarkdown::render(paste0(qcdir, "/MAVEQC_report.Rmd"))
        invisible(file.remove(paste0(qcdir, "/MAVEQC_report.Rmd")))
}