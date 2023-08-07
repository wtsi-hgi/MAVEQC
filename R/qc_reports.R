#' barplot in reactable
#'
#' @name bar_style
#' @param length   the length of bar
#' @param height   the height of bar
#' @param fill     the bar color
#' @param align    the alignment of bar
#' @param color    the font color
#' @param fweight  the font weight
bar_style <- function(length = 1,
                      height = "80%",
                      fill = "#00FFFF7F",
                      align = c("right", "left"),
                      color = "black",
                      fweight = "plain") {
    align <- match.arg(align)
    if (align == "left") {
        position <- paste0(length * 100, "%")
        image <- sprintf("linear-gradient(90deg, %1$s %2$s, transparent %2$s)", fill, position)
    } else {
        position <- paste0(100 - length * 100, "%")
        image <- sprintf("linear-gradient(90deg, transparent %1$s, %2$s %1$s)", position, fill)
    }
    list(
        backgroundImage = image,
        backgroundSize = paste("100%", height),
        backgroundRepeat = "no-repeat",
        backgroundPosition = "center",
        color = color,
        fontWeight = fweight
    )
}

#' create QC reports
#'
#' @export
#' @name create_qc_reports
#' @param samplesheet    the path of sample sheet file
#' @param qc_type         screen or plasmid
#' @param qc_dir          the directory of QC plots and outs
create_qc_reports <- function(samplesheet = NULL,
                              qc_type = c("plasmid", "screen"),
                              qc_dir = NULL) {
        #----------#
        # checking #
        #----------#
        if (is.null(samplesheet)) {
            stop(paste0("====> Error: please provide the path of sample sheet file!"))
        }

        if (is.null(qc_dir)) {
            stop(paste0("====> Error: qc_dir is not provided, no output directory."))
        }

        if (!file.exists(paste0(qc_dir, "/sample_qc_cutoffs.tsv"))) {
            stop(paste0("====> Error: sample_qc_cutoffs.tsv is not in ", qc_dir, ". Please use qcout_samqc_cutoffs to create it."))
        }

        qc_type <- match.arg(qc_type)

        #------------------#
        # creating reports #
        #------------------#
        package_version <- paste0("MAVEQC", "-v", packageVersion("MAVEQC"))
        report_path <- paste0(qc_dir, "/", "MAVEQC_report.Rmd")
        sink(report_path)

        cat("---", "\n", sep = "")
        cat("title: \"MAVE QC Report\"", "\n", sep = "")
        cat("author: \"", package_version, "\"", "\n", sep = "")
        cat("date: \"`r Sys.time()`\"", "\n", sep = "")
        cat("output:", "\n", sep = "")
        cat("    html_document:", "\n", sep = "")
        cat("        toc: true", "\n", sep = "")
        cat("        toc_depth: 4", "\n", sep = "")
        cat("        theme: united", "\n", sep = "")
        cat("        highlight: tango", "\n", sep = "")
        cat("---", "\n", sep = "")
        cat("\n", sep = "")

        cat("```{r setup, include = FALSE}", "\n", sep = "")
        cat("knitr::opts_chunk$set(echo = TRUE, fig.align = \"center\")", "\n", sep = "")
        cat("library(reactable)", "\n", sep = "")
        cat("outdir <- \"", qc_dir, "\"", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("\n", sep = "")

        cat("---", "\n", sep = "")
        cat("\n", sep = "")

        cat("```{js, echo = FALSE}", "\n", sep = "")
        cat("function formatNumber(num, precision = 1) {", "\n", sep = "")
        cat("    const map = [", "\n", sep = "")
        cat("        { suffix: 'T', threshold: 1e12 },", "\n", sep = "")
        cat("        { suffix: 'B', threshold: 1e9 },", "\n", sep = "")
        cat("        { suffix: 'M', threshold: 1e6 },", "\n", sep = "")
        cat("        { suffix: 'K', threshold: 1e3 },", "\n", sep = "")
        cat("        { suffix: '', threshold: 1 },", "\n", sep = "")
        cat("    ];", "\n", sep = "")
        cat("    const found = map.find((x) => Math.abs(num) >= x.threshold);", "\n", sep = "")
        cat("    if (found) {", "\n", sep = "")
        cat("        const formatted = (num / found.threshold).toFixed(precision) + found.suffix;", "\n", sep = "")
        cat("        return formatted;", "\n", sep = "")
        cat("    }", "\n", sep = "")
        cat("    return num;", "\n", sep = "")
        cat("}", "\n", sep = "")
        cat("", "\n", sep = "")
        cat("function rangeFilter(column, state) {", "\n", sep = "")
        cat("    let min = Infinity", "\n", sep = "")
        cat("    let max = 0", "\n", sep = "")
        cat("    state.data.forEach(function(row) {", "\n", sep = "")
        cat("        const value = row[column.id]", "\n", sep = "")
        cat("        if (value < min) {", "\n", sep = "")
        cat("            min = Math.floor(value)", "\n", sep = "")
        cat("        }", "\n", sep = "")
        cat("        if (value > max) {", "\n", sep = "")
        cat("            max = Math.ceil(value)", "\n", sep = "")
        cat("        }", "\n", sep = "")
        cat("    })", "\n", sep = "")
        cat("", "\n", sep = "")
        cat("    const filterValue = column.filterValue || min", "\n", sep = "")
        cat("    const input = React.createElement('input', {", "\n", sep = "")
        cat("        type: 'range',", "\n", sep = "")
        cat("        value: filterValue,", "\n", sep = "")
        cat("        min: min,", "\n", sep = "")
        cat("        max: max,", "\n", sep = "")
        cat("        onChange: function(event) {", "\n", sep = "")
        cat("            column.setFilter(event.target.value || undefined)", "\n", sep = "")
        cat("        },", "\n", sep = "")
        cat("        style: { width: '100%', marginRight: '8px' },", "\n", sep = "")
        cat("        'aria-label': 'Filter ' + column.name", "\n", sep = "")
        cat("    })", "\n", sep = "")
        cat("", "\n", sep = "")
        cat("    return React.createElement(", "\n", sep = "")
        cat("        'div',", "\n", sep = "")
        cat("        { style: { display: 'flex', alignItems: 'center', height: '100%' } },", "\n", sep = "")
        cat("        [input, formatNumber(filterValue)]", "\n", sep = "")
        cat("    )", "\n", sep = "")
        cat("}", "\n", sep = "")
        cat("", "\n", sep = "")
        cat("function filterMinValue(rows, columnId, filterValue) {", "\n", sep = "")
        cat("    return rows.filter(function(row) {", "\n", sep = "")
        cat("        return row.values[columnId] >= filterValue", "\n", sep = "")
        cat("    })", "\n", sep = "")
        cat("}", "\n", sep = "")
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

        if (qc_type == "screen") {
            cat("## 2. Screen QC", "\n", sep = "")
        } else {
            cat("## 2. Plasmid QC", "\n", sep = "")
        }
        cat("Displays QC plots and statistics for all samples for QC.", "\n", sep = "")
        cat("\n", sep = "")
        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("samqc_cutoffs <- as.data.frame(read.table(\"", qc_dir, "/sample_qc_cutoffs.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        cat("### 2.1. Sample Sheet", "\n", sep = "")
        cat("\n", sep = "")
        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", samplesheet, "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
        cat("          filterable = TRUE, minRows = 10, defaultColDef = colDef(minWidth = 150),", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")))", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        cat("### 2.2. Run Sample QC", "\n", sep = "")
        cat("\n", sep = "")
        cat("#### 2.2.1. Read Length Distrubtion", "\n", sep = "")
        cat("Displays the percentage of reads for each sample, based on 50 nucleotide increments.", "\n", sep = "")
        cat("\n", sep = "")

        cat("```{r, echo = FALSE, out.height = \"50%\", out.width = \"50%\"}", "\n", sep = "")
        cat("knitr::include_graphics(paste0(outdir, \"/sample_qc_read_length.png\"), rel_path = FALSE)", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", qc_dir, "/sample_qc_read_length.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
        cat("          filterable = TRUE, minRows = 10, defaultColDef = colDef(minWidth = 150),", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
        cat("          columns = list(\"% 0 ~ 50\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeFilter\"),", "\n", sep = "")
        cat("                                               style = function(value) {bar_style(length = value/100)}),", "\n", sep = "")
        cat("                         \"% 50 ~ 100\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeFilter\"),", "\n", sep = "")
        cat("                                                 style = function(value) {bar_style(length = value/100)}),", "\n", sep = "")
        cat("                         \"% 100 ~ 150\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeFilter\"),", "\n", sep = "")
        cat("                                                  style = function(value) {bar_style(length = value/100)}),", "\n", sep = "")
        cat("                         \"% 150 ~ 200\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeFilter\"),", "\n", sep = "")
        cat("                                                  style = function(value) {bar_style(length = value/100)}),", "\n", sep = "")
        cat("                         \"% 200 ~ 250\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeFilter\"),", "\n", sep = "")
        cat("                                                  style = function(value) {bar_style(length = value/100)}),", "\n", sep = "")
        cat("                         \"% 250 ~ 300\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeFilter\"),", "\n", sep = "")
        cat("                                                  style = function(value) {bar_style(length = value/100)}),", "\n", sep = "")
        cat("                         \"Total Reads\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeFilter\"),", "\n", sep = "")
        cat("                                                  format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"Pass\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" })))", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        cat("#### 2.2.2. Total Reads", "\n", sep = "")
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
        cat("df <- as.data.frame(read.table(\"", qc_dir, "/sample_qc_stats_total.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
        cat("          filterable = TRUE, minRows = 10, defaultColDef = colDef(minWidth = 150),", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
        cat("          columns = list(\"Accepted Reads\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeFilter\"),", "\n", sep = "")
        cat("                                                     format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"% Accepted Reads\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeFilter\"),", "\n", sep = "")
        cat("                                                       style = function(value) {bar_style(length = value/100)}),", "\n", sep = "")
        cat("                         \"Excluded Reads\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeFilter\"),", "\n", sep = "")
        cat("                                                     format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"% Excluded Reads\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeFilter\"),", "\n", sep = "")
        cat("                                                       style = function(value) {bar_style(length = value/100)}),", "\n", sep = "")
        cat("                         \"Total Reads\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeFilter\"),", "\n", sep = "")
        cat("                                                  format = colFormat(separators = TRUE),", "\n", sep = "")
        cat("                                                  style = function(value) { ", "\n", sep = "")
        cat("                                                              if (value < samqc_cutoffs$total_reads) {", "\n", sep = "")
        cat("                                                                  color <- \"red\"", "\n", sep = "")
        cat("                                                                  fweight <- \"bold\"", "\n", sep = "")
        cat("                                                              } else {", "\n", sep = "")
        cat("                                                                  color <- \"forestgreen\"", "\n", sep = "")
        cat("                                                                  fweight <- \"bold\" }", "\n", sep = "")
        cat("                                                              list(color = color, fontWeight = fweight)}),", "\n", sep = "")
        cat("                         \"Pass Threshold\" = colDef(format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"Pass\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" })))", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        cat("#### 2.2.3. Accepted Reads", "\n", sep = "")
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
        cat("df <- as.data.frame(read.table(\"", qc_dir, "/sample_qc_stats_accepted.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
        cat("          filterable = TRUE, minRows = 10, defaultColDef = colDef(minWidth = 150),", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
        cat("          columns = list(\"% Library Reads\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeFilter\"),", "\n", sep = "")
        cat("                                                      style = function(value) {", "\n", sep = "")
        cat("                                                                  if (value < samqc_cutoffs$library_percent * 100) {", "\n", sep = "")
        cat("                                                                      bar_style(length = value/100, color = \"red\", fweight = \"bold\")", "\n", sep = "")
        cat("                                                                  } else {", "\n", sep = "")
        cat("                                                                      bar_style(length = value/100, color = \"forestgreen\", fweight = \"bold\") }}),", "\n", sep = "")
        cat("                         \"% Reference Reads\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeFilter\"),", "\n", sep = "")
        cat("                                                        style = function(value) {bar_style(length = value/100)}),", "\n", sep = "")
        cat("                         \"% PAM Reads\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeFilter\"),", "\n", sep = "")
        cat("                                                  style = function(value) {bar_style(length = value/100)}),", "\n", sep = "")
        cat("                         \"% Unmapped Reads\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeFilter\"),", "\n", sep = "")
        cat("                                                       style = function(value) {bar_style(length = value/100)}),", "\n", sep = "")
        cat("                         \"Pass\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" })))", "\n", sep = "")
        cat("```", "\n", sep = "")

        cat("\n", sep = "")
        cat("Defines the mean read count per template oligo sequence.", "\n", sep = "")
        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", qc_dir, "/sample_qc_stats_coverage.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
        cat("          filterable = TRUE, minRows = 10, defaultColDef = colDef(minWidth = 150),", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
        cat("          columns = list(\"Total Library Reads\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeFilter\"),", "\n", sep = "")
        cat("                                                          format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"Total Template Oligo Sequences\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeFilter\"),", "\n", sep = "")
        cat("                                                                     format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"Library Coverage\" = colDef(filterMethod = JS(\"filterMinValue\"), filterInput = JS(\"rangeFilter\"),", "\n", sep = "")
        cat("                                                       format = colFormat(separators = TRUE),", "\n", sep = "")
        cat("                                                       style = function(value) {", "\n", sep = "")
        cat("                                                                   if (value < samqc_cutoffs$library_cov) {", "\n", sep = "")
        cat("                                                                       color <- \"red\"", "\n", sep = "")
        cat("                                                                       fweight <- \"bold\"", "\n", sep = "")
        cat("                                                                   } else {", "\n", sep = "")
        cat("                                                                       color <- \"forestgreen\"", "\n", sep = "")
        cat("                                                                       fweight <- \"bold\" }", "\n", sep = "")
        cat("                                                                   list(color = color, fontWeight = fweight)}),", "\n", sep = "")
        cat("                         \"Pass\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" })))", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        cat("#### 2.2.4. Genomic Coverage", "\n", sep = "")
        cat("Distribution of variants across targeton region based on log2(count+1) values.", "\n", sep = "")
        cat("\n", sep = "")

        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("knitr::include_graphics(paste0(outdir, \"/sample_qc_position_cov.dots.png\"), rel_path = FALSE)", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("```{r, echo = FALSE}", "\n", sep = "")
        cat("df <- as.data.frame(read.table(\"", qc_dir, "/sample_qc_stats_pos_coverage.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
        cat("reactable(df, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
        cat("          filterable = TRUE, minRows = 10, defaultColDef = colDef(minWidth = 150),", "\n", sep = "")
        cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
        cat("          columns = list(\"Genomic Start\" = colDef(format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"Genomic End\" = colDef(format = colFormat(separators = TRUE)),", "\n", sep = "")
        cat("                         \"% Low Abundance\" = colDef(minWidth = 150,", "\n", sep = "")
        cat("                                                      style = function(value) {", "\n", sep = "")
        cat("                                                                  if (value > (1 - samqc_cutoffs$low_abundance_lib_per) * 100) {", "\n", sep = "")
        cat("                                                                      color <- \"red\"", "\n", sep = "")
        cat("                                                                      fweight <- \"bold\"", "\n", sep = "")
        cat("                                                                  } else {", "\n", sep = "")
        cat("                                                                      color <- \"forestgreen\"", "\n", sep = "")
        cat("                                                                      fweight <- \"bold\" }", "\n", sep = "")
        cat("                                                                  list(color = color, fontWeight = fweight)}),", "\n", sep = "")
        cat("                         \"Low Abundance cutoff\" = colDef(minWidth = 150),", "\n", sep = "")
        cat("                         \"Pass\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" })))", "\n", sep = "")
        cat("```", "\n", sep = "")
        cat("<br>", "\n", sep = "")
        cat("\n", sep = "")

        if (qc_type == "screen") {
            cat("#### 2.2.5. Genomic Position Percentage", "\n", sep = "")
            cat("Displays distribution of \"LOF\" (loss-of-function) vs all \"Other\" variants across the targeton region, ",
                "based on read percentages for reference timepoint. ",
                "Requires concordant distribution of LOF and Other variants.", "\n", sep = "")
            cat("\n", sep = "")

            cat("```{r, echo = FALSE, out.height = \"65%\", out.width = \"65%\"}", "\n", sep = "")
            cat("knitr::include_graphics(paste0(outdir, \"/sample_qc_position_anno.lof_dots.png\"), rel_path = FALSE)", "\n", sep = "")
            cat("```", "\n", sep = "")
            cat("```{r, echo = FALSE}", "\n", sep = "")
            cat("df <- as.data.frame(read.table(\"", qc_dir, "/sample_qc_stats_pos_percentage.tsv", "\", header = TRUE, sep = \"\\t\", check.names = FALSE))", "\n", sep = "")
            cat("reactable(df, highlight = TRUE, bordered = TRUE,  striped = TRUE, compact = TRUE, wrap = TRUE,", "\n", sep = "")
            cat("          filterable = TRUE, minRows = 10, defaultColDef = colDef(minWidth = 150),", "\n", sep = "")
            cat("          theme = reactableTheme(style = list(fontFamily = \"-apple-system\", fontSize = \"0.85em\")),", "\n", sep = "")
            cat("          columns = list(\"Genomic Start\" = colDef(format = colFormat(separators = TRUE)),", "\n", sep = "")
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
            cat("                                                                            fweight <- \"bold\" }", "\n", sep = "")
            cat("                                                                        list(color = color, fontWeight = fweight)}),", "\n", sep = "")
            cat("                         \"% Low Abundance cutoff\" = colDef(minWidth = 150),", "\n", sep = "")
            cat("                         \"Pass\" = colDef(cell = function(value) { if (value) \"\\u2705\" else \"\\u274c\" })))", "\n", sep = "")
            cat("```", "\n", sep = "")
            cat("<br>", "\n", sep = "")
            cat("\n", sep = "")

            cat("### 2.3. Run Experiment QC", "\n", sep = "")
            cat("\n", sep = "")
            cat("#### 2.3.1. Sample Correlations", "\n", sep = "")
            cat("\n", sep = "")
            cat("```{r, echo = FALSE, out.height = \"50%\", out.width = \"50%\"}", "\n", sep = "")
            cat("knitr::include_graphics(paste0(outdir, \"/sample_qc_samples_tree.png\"), rel_path = FALSE)", "\n", sep = "")
            cat("knitr::include_graphics(paste0(outdir, \"/sample_qc_samples_corr.png\"), rel_path = FALSE)", "\n", sep = "")
            cat("```", "\n", sep = "")
            cat("<br>", "\n", sep = "")
            cat("\n", sep = "")

            cat("#### 2.3.2. Sample PCA", "\n", sep = "")
            cat("\n", sep = "")
            cat("```{r, echo = FALSE, out.height = \"75%\", out.width = \"75%\"}", "\n", sep = "")
            cat("knitr::include_graphics(paste0(outdir, \"/sample_qc_pca_samples.png\"), rel_path = FALSE)", "\n", sep = "")
            cat("```", "\n", sep = "")
            cat("<br>", "\n", sep = "")
            cat("\n", sep = "")

            cat("#### 2.3.3. Folder Change (by category)", "\n", sep = "")
            cat("\n", sep = "")
            cat("```{r, echo = FALSE, out.height = \"75%\", out.width = \"75%\"}", "\n", sep = "")
            cat("figs <- list.files(path = outdir, pattern = \"sample_qc_deseq_fc.*.violin.png\", full.names = TRUE)", "\n", sep = "")
            cat("figs <- mixedsort(figs)", "\n", sep = "")
            cat("knitr::include_graphics(figs, rel_path = FALSE)", "\n", sep = "")
            cat("```", "\n", sep = "")
            cat("<br>", "\n", sep = "")
            cat("\n", sep = "")

            cat("#### 2.3.4. Folder Change (by position)", "\n", sep = "")
            cat("\n", sep = "")
            cat("```{r, echo = FALSE, out.height = \"75%\", out.width = \"75%\"}", "\n", sep = "")
            cat("figs <- list.files(path = outdir, pattern = \"sample_qc_deseq_fc.*.position.png\", full.names = TRUE)", "\n", sep = "")
            cat("figs <- mixedsort(figs)", "\n", sep = "")
            cat("knitr::include_graphics(figs, rel_path = FALSE)", "\n", sep = "")
            cat("```", "\n", sep = "")
            cat("<br>", "\n", sep = "")
            cat("\n", sep = "")
        }

        sink()

        rmarkdown::render(paste0(qc_dir, "/MAVEQC_report.Rmd"), clean = TRUE, quiet = TRUE)
        invisible(file.remove(paste0(qc_dir, "/MAVEQC_report.Rmd")))
}