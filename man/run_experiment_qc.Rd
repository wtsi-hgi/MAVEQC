% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/experimentQC.R
\name{run_experiment_qc}
\alias{run_experiment_qc}
\title{initialize function}
\usage{
run_experiment_qc(object, ...)

\S4method{run_experiment_qc}{experimentQC}(
  object,
  pcut = maveqc_config$expqc_padj,
  dcut = maveqc_config$expqc_lfc_depleted,
  ecut = maveqc_config$expqc_lfc_enriched,
  ntop = maveqc_config$expqc_top_variants
)
}
\arguments{
\item{object}{experimentQC object}

\item{pcut}{the padj cutoff}

\item{dcut}{the depleted log2 fold change cutoff}

\item{ecut}{the enriched log2 fold change cutoff}

\item{ntop}{the number of top variances}
}
\value{
object
}
\description{
initialize function

run DESeq2 for the list of samples
}
