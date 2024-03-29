% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RMiXCR_runners.R
\name{RunMiXCR}
\alias{RunMiXCR}
\title{RunMiXCR}
\usage{
\dontrun{
fasta.file <- "/home/mysubjects/subject1.R1.fastq"
RunMiXCR(sbj=fasta.file, 
         species="hsa", #human
         assemblePartial=FALSE, 
         extendedAlignments=FALSE, 
         nThreads = 6L #six cores
         )

}
}
\arguments{
\item{sbj}{(character) The full path name of the sbj.R1.fastq file (the paired read file should be named as sbj.R2.fastq)}

\item{species}{any of "hsa","mmu" (human - mouse)}

\item{assemblePartial}{(logical, default FALSE) if partial read should be assembled (Perform two rounds of contig assembly)}

\item{extendedAlignments}{(logical, default TRUE) if TCR extension alignment should be performed (in case of incomplete TCR alignments)}

\item{nThreads}{(integer) number of cores}
}
\value{
a list with clones, vmin, vdetails.clones, TRA.clones,TRB.clones,TRG.clones,TRD.clones data frames
it also save an Excel file as full.path.sbj.alignments.vdjca.xlsx
}
\description{
RunMiXCR
}
