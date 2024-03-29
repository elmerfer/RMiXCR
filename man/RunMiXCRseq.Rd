% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RMiXCR_runners.R
\name{RunMiXCRseq}
\alias{RunMiXCRseq}
\title{RunMiXCRseq}
\usage{
RunMiXCRseq(sbj, species = c("hsa","mmu"), extendedAlignments = TRUE )
}
\arguments{
\item{sbj}{sbj (character) The full path name of the sbj.R1.fastq file (the paired read file should be named as sbj.R2.fastq)}

\item{species}{character any of "hsa","mmu" (default hsa)}

\item{extendedAlignments}{default FALSE,}
}
\value{
if success, the it builds three files sbj.alignments.vdjca and sbj.clones.clns and clones.xlsx and return a vector with the full path to files
#' if fail, returns NA
}
\description{
It runs MiXCR for RNAseq data according to MiXCR main page advise. In this implementation it runs first \code{\link{RunMiXCRseqAlignment}} and the \code{\link{RunMiXCRseqClones}}
one after the other. See the other helps for more details
}
\examples{
/dontrun{
      subj.file <- "/home/.../subject.R1.fastq"
      out.vdja.and.clone.files <- RunMiXCRseq(sbj = subj.file)
      print(out.vdja.and.clone.files[1])#vdja file path
      print(out.vdja.and.clone.files[2])#clones file path
      print(out.vdja.and.clone.files[3])#clones excel report file path
}
}
