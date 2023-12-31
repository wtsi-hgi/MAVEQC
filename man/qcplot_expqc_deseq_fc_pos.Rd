% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc_plots.R
\name{qcplot_expqc_deseq_fc_pos}
\alias{qcplot_expqc_deseq_fc_pos}
\title{initialize function}
\usage{
qcplot_expqc_deseq_fc_pos(object, ...)

\S4method{qcplot_expqc_deseq_fc_pos}{experimentQC}(
  object,
  eqc_type = c("lib", "all"),
  cons = c("Synonymous_Variant", "LOF", "Missense_Variant", "Intronic_Variant",
    "Inframe_Deletion", "Splice_Variant", "Splice_Polypyrimidine_Tract_Variant",
    "Others"),
  plot_dir = NULL
)
}
\arguments{
\item{object}{experimentQC object}

\item{eqc_type}{library counts or all counts}

\item{cons}{a vector of all the consequences in the vep annotation file}

\item{plot_dir}{the output plot directory}
}
\description{
initialize function

create fold change and consequence plot
}
