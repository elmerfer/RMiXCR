% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RMiXCR_reports.R
\name{Report}
\alias{Report}
\title{Report}
\usage{
Report(clonsFile, type = c("full","min"), otherFields, fileType = c("Excel","tab"))
}
\arguments{
\item{clonsFile}{the full path to the clone file (subjsect.R1.fastq.clones.clns) as returned by \code{\link{RunMiXCRseqClones}}}

\item{type}{string "full" (default) or "min"}

\item{otherFields}{other fields (not implemented yet)}

\item{fileType}{Excel or tab (tab delimited)}
}
\description{
In order to export alignment results or clones from a binary file (.clns) to a human-readable tab or excel file
}
\details{
The resulting tab-delimited text file or excel file (type param) will contain columns with different types of information. 
If no options are specified (otherFields), the default set of columns - which is sufficient in most cases - will be exported. 
The possible columns include. See details in [MiXCR main page](https://mixcr.readthedocs.io/en/master/export.html): 
aligned sequences, qualities, all or just best hit for V, D, J and C genes, corresponding alignments, nucleotide and amino acid sequences of gene region present in sequence, etc. When exporting clones, the additional columns include: clone count, clone fraction etc.
One can customize the list of fields that will be exported by passing parameters to export commands (not implemented yet).
}
