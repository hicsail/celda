% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diffExp.R
\name{diffExp}
\alias{diffExp}
\title{Gene expression markers for cell clusters using MAST}
\usage{
diffExp(counts, celda.mod, c1, c2 = NULL, only.pos = FALSE,
  log2fc.threshold = NULL, fdr.threshold = 1)
}
\arguments{
\item{counts}{A numeric count matrix.}

\item{celda.mod}{An object of class celda_C or celda_CG.}

\item{c1}{Numeric. Cell cluster(s) to define markers for.}

\item{c2}{Numeric. Second cell cluster(s) for comparison. If NULL (default) - use all
other cells for comparison.}

\item{only.pos}{Logical. Only return markers with positive log2fc (FALSE by default).}

\item{log2fc.threshold}{Numeric. Only return markers whose the absolute values of log2fc
are greater than this threshold (NULL by default).}

\item{fdr.threshold}{Numeric. Only return markers whose false discovery rates (FDRs) are less
than this threshold (1 by default).}
}
\value{
Data frame containing a ranked list (based on the absolute value of log2fc) of putative markers,
   and associated statistics (p-value, log2fc and FDR).
}
\description{
Finds markers (differentially expressed genes) for cell clusters
   using MAST: a flexible statistical framework for assessing transcriptional
   changes and characterizing heterogeneity in single-cell RNA sequencing data
   (Finak et al, Genome Biology, 2015)
}
