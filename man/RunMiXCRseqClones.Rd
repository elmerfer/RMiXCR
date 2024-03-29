% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RMiXCR_runners.R
\name{RunMiXCRseqClones}
\alias{RunMiXCRseqClones}
\title{RunMiXCRseqClones}
\usage{
RunMiXCRseqClones(alignmentFile, nThreads = 4L)
}
\arguments{
\item{alignmentFile}{the full path name of the vdja file returned by \code{\link{RunMiXCRseqAligment}}}

\item{nThreads}{(integer) number of cores (4 default as it is suggested as the optimal value)}
}
\value{
it will generate a file names xxx.R1.fastq.clones.clns and it will return the full path to such file
}
\description{
The assemble command builds clonotypes from alignments obtained with \code{\link{RunMiXCR}} or
\code{\link{RunMiXCRseqAligment}}. Clonotypes assembly is performed for a chosen assembling feature (e.g. CDR3 by default).
}
