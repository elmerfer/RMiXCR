
.mixcrTest <- function(cmd){
  software <- RMiXCR:::.OpenConfigFile()
  system2(command= "java",
          args = c(software$mixcr$base.args, cmd))
  
}

#'RunMiXCR
#'
#'@param sbj (character) The full path name of the sbj.R1.fastq file (the paired read file should be named as sbj.R2.fastq)
#'@param species any of "hsa","mmu" (human - mouse)
#'@param assemblePartial (logical, default FALSE) if partial read should be assembled (Perform two rounds of contig assembly)
#'@param extendedAlignments (logical, default TRUE) if TCR extension alignment should be performed (in case of incomplete TCR alignments)
#'@param nThreads (integer) number of cores
#'@export
#'@return a list with clones, vmin, vdetails.clones, TRA.clones,TRB.clones,TRG.clones,TRD.clones data frames
#'it also save an Excel file as full.path.sbj.alignments.vdjca.xlsx
#'@usage
#'\dontrun{
#' fasta.file <- "/home/mysubjects/subject1.R1.fastq"
#' RunMiXCR(sbj=fasta.file, 
#'          species="hsa", #human
#'          assemblePartial=FALSE, 
#'          extendedAlignments=FALSE, 
#'          nThreads = 6L #six cores
#'          )
#'
#'}
#'

RunMiXCR <- function(sbj, species = c("hsa","mmu"), assemblePartial=TRUE,extendedAlignments=FALSE, nThreads){
  software <- RMiXCR:::.OpenConfigFile()  
  cat(paste0("\nAlignment"))
  if(missing(nThreads)){
    nt <- parallel::detectCores()
    nThreads <- ifelse(nt-2 < 1, 1, nt)
    message(paste0("\nRunning with ", nThreads, " Cores"))
  }
  species <- match.arg(species, c("hsa","mmu") )
  file1 <- sbj
  file2 <- stringr::str_replace(file1,".R1.fasta",".R2.fasta")
  align.ofile <- file.path(dirname(file1), paste0(basename(file1),".alignments.vdjca"))
  cat("\nAlign sequencing reads against reference V, D, J and C genes.")
  system2(command= "java",
          args = c(software$mixcr$base.args,"align",
                   paste0("--species ",species),
                   "--parameters rna-seq",
                   "-OallowPartialAlignments=true",
                   paste0("--threads ",nThreads),
                   file1,
                   file2,
                   align.ofile))
  if(file.exists(align.ofile)==FALSE){
    stop(align.ofile)
  }      
  
  if(assemblePartial){
    # mixcr assemblePartial alignments.vdjca alignmentsRescued.vdjca
    cat("\nAssembling parial reads")
    apartial.input.file <- align.ofile
    apartial.out.file <- file.path(dirname(file1), paste0(basename(file1),".alignments.partial.vdjca")) 
    for( i in 1:2){## perform two rounds of contig assembly 
      system2(command= "java",
              args = c(software$mixcr$base.args,
                       "assemblePartial",
                       apartial.input.file,
                       apartial.out.file))
      apartial.input.file <- apartial.out.file
    }
    if(extendedAlignments){
      cat("\nExtending TCR alignments")
      system2(command= "java",
              args = c(software$mixcr$base.args,
                       "extendAlignments",
                       apartial.out.file,
                       apartial.out.file))
    }
    align.ofile <- apartial.out.file
  }
  
  assemble.ofile <- file.path(dirname(file1), paste0(basename(file1),".clones.clns"))
  cat(paste0("\nAssembling"))
  system2(command = "java",
          args = c(software$mixcr$base.args,
                   "assemble",
                   paste0("--threads ",nThreads),
                   "-OmaxBadPointsPercent=0",
                   align.ofile,
                   assemble.ofile
          ))
  if(file.exists(assemble.ofile)==FALSE){
    stop(assemble.ofile)
  }
  cat(paste0("\nExporting clones to Excel file"))
  # java -jar -Xmx6g $MIXCR_BIN_DIR/mixcr.jar exportClones --chains TCR ${OUTPUT_FILE}.clones.clns ${OUTPUT_FILE}.clones.txt
  system2(command = "java",
          args = c(software$mixcr$base.args,
                   "exportClones", "--chains TCR",
                   assemble.ofile, 
                   stringr::str_replace(assemble.ofile,".clns",".tab")
          ))
  # java -jar -Xmx6g $MIXCR_BIN_DIR/mixcr.jar exportClones --chains TCR -p min ${OUTPUT_FILE}.clones.clns ${OUTPUT_FILE}.vmin.clones.txt
  system2(command = "java",
          args = c(software$mixcr$base.args,
                   "exportClones", "--chains TCR",
                   "-p min",
                   assemble.ofile, 
                   stringr::str_replace(assemble.ofile,".clns","vmin.tab")
          ))
  # java -jar -Xmx6g $MIXCR_BIN_DIR/mixcr.jar exportClones -count -vHit -dHit -jHit -cHit -vHits -dHits -jHits -cHits -nFeature CDR3 -aaFeature CDR3 ${OUTPUT_FILE}.clones.clns ${OUTPUT_FILE}.vdetails.clones.txt
  system2(command = "java",
          args = c(software$mixcr$base.args,
                   "exportClones", "-count",
                   "-vhit", "-dhit", "-jHit", "-cHit", "-vHits", "-dHits", "-jHits", "-cHits",
                   "-nFeature CDR3", "-aaFeature CDR3",
                   assemble.ofile, 
                   stringr::str_replace(assemble.ofile,".clns",".vdetails.clones.tab")
          ))
  # java -jar -Xmx6g $MIXCR_BIN_DIR/mixcr.jar exportClones -count -vHit -dHit -jHit -cHit -vHits -dHits -jHits -cHits -nFeature CDR3 -c TRA -aaFeature CDR3 ${OUTPUT_FILE}.clones.clns ${OUTPUT_FILE}.TRA.clones.txt
  system2(command = "java",
          args = c(software$mixcr$base.args,
                   "exportClones", "-count",
                   "-vHit", "-dHit", "-jHit", "-cHit", "-vHits", "-dHits", "-jHits", "-cHits",
                   "-nFeature CDR3", "-c TRA", "-aaFeature CDR3",
                   assemble.ofile, 
                   stringr::str_replace(assemble.ofile,".clns",".TRA.clones.tab")
          ))
  
  # java -jar -Xmx6g $MIXCR_BIN_DIR/mixcr.jar exportClones -count -vHit -dHit -jHit -cHit -vHits -dHits -jHits -cHits -nFeature CDR3 -c TRB -aaFeature CDR3 ${OUTPUT_FILE}.clones.clns ${OUTPUT_FILE}.TRB.clones.txt
  system2(command = "java",
          args = c(software$mixcr$base.args,
                   "exportClones", "-count",
                   "-vHit", "-dHit", "-jHit", "-cHit", "-vHits", "-dHits", "-jHits", "-cHits",
                   "-nFeature CDR3", "-c TRB", "-aaFeature CDR3",
                   assemble.ofile, 
                   stringr::str_replace(assemble.ofile,".clns",".TRB.clones.tab")
          ))
  # java -jar -Xmx6g $MIXCR_BIN_DIR/mixcr.jar exportClones -count -vHit -dHit -jHit -cHit -vHits -dHits -jHits -cHits -nFeature CDR3 -c TRG -aaFeature CDR3 ${OUTPUT_FILE}.clones.clns ${OUTPUT_FILE}.TRG.clones.txt
  system2(command = "java",
          args = c(software$mixcr$base.args,
                   "exportClones", "-count",
                   "-vHit", "-dHit", "-jHit", "-cHit", "-vHits", "-dHits", "-jHits", "-cHits",
                   "-nFeature CDR3", "-c TRG", "-aaFeature CDR3",
                   assemble.ofile, 
                   stringr::str_replace(assemble.ofile,".clns",".TRG.clones.tab")
          ))
  # java -jar -Xmx6g $MIXCR_BIN_DIR/mixcr.jar exportClones -count -vHit -dHit -jHit -cHit -vHits -dHits -jHits -cHits -nFeature CDR3 -c TRD -aaFeature CDR3 ${OUTPUT_FILE}.clones.clns ${OUTPUT_FILE}.TRD.clones.txt
  system2(command = "java",
          args = c(software$mixcr$base.args,
                   "exportClones", "-count",
                   "-vHit", "-dHit", "-jHit", "-cHit", "-vHits", "-dHits", "-jHits", "-cHits",
                   "-nFeature CDR3", "-c TRD", "-aaFeature CDR3",
                   assemble.ofile, 
                   stringr::str_replace(assemble.ofile,".clns",".TRD.clones.tab")
          ))
  
  to.excel <- lapply(c(".tab","vmin.tab",".vdetails.clones.tab",".TRA.clones.tab",".TRB.clones.tab",".TRG.clones.tab",".TRD.clones.tab"), function(x){
    if(file.exists(stringr::str_replace(assemble.ofile,".clns",x))){
      return(read.table(stringr::str_replace(assemble.ofile,".clns",x),h=T, sep="\t"))  
    }else{
      return(NA)
    }
    
  })
  aux.files.subfix <- c(".tab","vmin.tab",".vdetails.clones.tab",".TRA.clones.tab",".TRB.clones.tab",".TRG.clones.tab",".TRD.clones.tab")
  names(to.excel) <- stringr::str_remove(c(".tab","vmin.tab",".vdetails.clones.tab",".TRA.clones.tab",".TRB.clones.tab",".TRG.clones.tab",".TRD.clones.tab"),".tab")
  openxlsx::write.xlsx(to.excel, stringr::str_replace(assemble.ofile,".clns",".xlsx"))
  file.remove(stringr::str_replace(assemble.ofile,".clns",aux.files.subfix))
  return(invisible(aux.files.subfix))      
}

#' RNAseqAlignment
#'@description  This function aligns raw sequencing reads to reference V, D, J and C genes of T- and B- cell receptors.
#'@usage RunMiXCRseqAlignment(sbj, species = c("hsa","mmu","rat"), extendedAlignments=FALSE, nThreads = 4L)
#'@param sbj full path to the subject.R1.fastq sample file. The paired read should be named subject.R2.fastq
#'@param species any of "hsa","mmu", "rat" (human - mouse, rat)
#'@param extendedAlignments (logical, default TRUE) if TCR extension alignment should be performed (in case of incomplete TCR alignments)
#'@param nThreads (integer) number of cores (default 4, since it is suggested as the optimum value)
#'@export
#'@return 
#'It will generate file named "subject.R1.fastq.alignments.vdjca"
#'the full path to this file is returned to be used as input to \code{\link{RunMiXCRseqClones}}
#'@examples 
#'/dontrun{
#'       subj.file <- "/home/.../subject.R1.fastq"
#'       out.vdja.file <- RunMiXCR(sbj = subj.file)
#'       print(out.vdja.file)
#'}
#'
#'
RunMiXCRseqAlignment <- function(sbj, species = c("hsa","mmu","rat"), extendedAlignments=FALSE, nThreads = 4L){
  software <- RMiXCR:::.OpenConfigFile()  
  cat(paste0("\nAlignment"))
  if(missing(nThreads)){
    nt <- parallel::detectCores()
    nThreads <- min(ifelse(nt-2 < 1, 1, nt),5)
    message(paste0("\nRunning with ", nThreads, " Cores"))
  }
  species <- match.arg(species, c("hsa","mmu","rat") )
  file1 <- sbj
  file2 <- stringr::str_replace(file1,".R1.fasta",".R2.fasta")
  align.ofile <- file.path(dirname(file1), paste0(basename(file1),".alignments.vdjca"))
  if(file.exists(align.ofile)){
    stop(paste0("\nFile laready exists", align.ofile,"\nPls rename or remove"))
  }
  cat("\nAlign sequencing reads against reference V, D, J and C genes.")
  system2(command= "java",
          args = c(software$mixcr$base.args,"align",
                   paste0("--species ",species),
                   "--parameters rna-seq",
                   "-OallowPartialAlignments=true","-f",
                   paste0("--threads ",nThreads),
                   file1,
                   file2,
                   align.ofile))
  if(file.exists(align.ofile)==FALSE){
    stop(paste0("\nSomething went wrong with ",align.ofile))
  }      
  
  
  cat(paste0("\nAlignment Partial Round 1"))
  
  apartial.input.file <- align.ofile
  rescued.1.file <- stringr::str_replace(apartial.input.file,".alignments.vdjca",".alignments.rescued1.vdjca") 
    
  system2(command= "java",
              args = c(software$mixcr$base.args,
                       "assemblePartial","-f",
                       align.ofile,
                       rescued.1.file))
  cat(paste0("\nAlignment Partial Round 2"))
  system2(command= "java",
          args = c(software$mixcr$base.args,
                   "assemblePartial","-f",
                   rescued.1.file,
                   align.ofile))
  
  file.remove(rescued.1.file)
   if(extendedAlignments){
      cat("\nExtending TCR alignments")
      system2(command= "java",
              args = c(software$mixcr$base.args,
                       "extendAlignments","-f",
                       align.ofile,
                       align.ofile))
    }
    
  attr(align.ofile,"extended") <- ifelse(extendedAlignments,TRUE,FALSE)
  return(align.ofile)
}

#' RunMiXCRseq
#' @description 
#' It runs MiXCR for RNAseq data according to MiXCR main page advise. In this implementation it runs first \code{\link{RunMiXCRseqAlignment}} and the \code{\link{RunMiXCRseqClones}}
#' one after the other. See the other helps for more details
#' @usage 
#' RunMiXCRseq(sbj, species = c("hsa","mmu"), extendedAlignments = TRUE )
#' @param sbj sbj (character) The full path name of the sbj.R1.fastq file (the paired read file should be named as sbj.R2.fastq)
#' @param species character any of "hsa","mmu" (default hsa)
#' @param extendedAlignments default TRUE, 
#' @export
#' @return if success, the it builds three files sbj.alignments.vdjca and sbj.clones.clns and clones.xlsx and return a vector with the full path to files
#' #' if fail, returns NA
#' @examples 
#'/dontrun{
#'       subj.file <- "/home/.../subject.R1.fastq"
#'       out.vdja.and.clone.files <- RunMiXCRseq(sbj = subj.file)
#'       print(out.vdja.and.clone.files[1])#vdja file path
#'       print(out.vdja.and.clone.files[2])#clones file path
#'       print(out.vdja.and.clone.files[3])#clones excel report file path
#'}
RunMiXCRseq <- function(sbj, species = c("hsa","mmu"), extendedAlignments = TRUE ){
  
  out.alignment <- RunMiXCRseqAlignment(sbj = sbj, species = species, extendedAlignments = extendedAlignments)
  out.clones <- RunMiXCRseqClones(alignmentFile = out.alignment)
  out.exel.report <- RunMiXCRreport(out.clones,fileType = "Excel")
  files <- c(Aligment= out.alignment, Clones=out.clones, ClonesExcelReport=out.excel.report)
  if(all(file.exists(files))==TRUE){
    message(paste0("\n files created",basename(files) ))
    return(files)
  }
  return(NA)
  
  
}

#' RNAseqClones
#' @description The assemble command builds clonotypes from alignments obtained with \code{\link{RunMiXCR}} or
#' \code{\link{RunMiXCRseqAligment}}. Clonotypes assembly is performed for a chosen assembling feature (e.g. CDR3 by default).
#'@param alignmentFile the full path name of the vdja file returned by \code{\link{RunMiXCRseqAligment}}
#'@param nThreads (integer) number of cores (4 default as it is suggested as the optimal value)
#'@export
#'@return it will generate a file names xxx.R1.fastq.clones.clns and it will return the full path to such file
#'@exaples
#'\dontrun{
#'      subj.file <- "/home/.../subject.R1.fastq"
#'       out.vdja.file <- RunMiXCR(sbj = subj.file)
#'       clns.file <- RNAseqClones(out.vdja.file)
#'       stropifnot(file.exists(clns.file)==TRUE)
#'}
#'
RunMiXCRseqClones <- function(alignmentFile, nThreads = 4L){
  software <- RMiXCR:::.OpenConfigFile()  
 
  if(missing(nThreads)){
    nt <- parallel::detectCores()
    nThreads <- min(ifelse(nt-2 < 1, 1, nt),5)
    message(paste0("\nRunning with ", nThreads, " Cores"))
  }  
  if(file.exists(alignmentFile)==FALSE){
    stop(paste0("\nERROR file not found\n",alignmentFile))
  }
  assemble.ofile <- stringr::str_replace(alignmentFile,".alignments.vdjca",".clones.clns")
  cat(paste0("\nAssembling"))
  system2(command = "java",
          args = c(software$mixcr$base.args,
                   "assemble",
                   paste0("--threads ",nThreads),
                   "-OmaxBadPointsPercent=0","-f",
                   alignmentFile,
                   assemble.ofile
          ))
  if(file.exists(assemble.ofile)==FALSE){
    stop(paste0("\nSomething went wrong with ",assemble.ofile))
  }
  return(assemble.ofile)
}

#' RunMiXCR.ExportExcel
#' @param clonsFile (string) the clone file. No IG clones
#' @export
RunMiXCRexportExcel <- function(clonsFile ){
  software <- RMiXCR:::.OpenConfigFile()
  cat(paste0("\nExporting clones to Excel file"))
  # java -jar -Xmx6g $MIXCR_BIN_DIR/mixcr.jar exportClones --chains TCR ${OUTPUT_FILE}.clones.clns ${OUTPUT_FILE}.clones.txt
  clonsFile <- clonsFile[stringr::str_detect(clonsFile,".clones.clns")]
  if(length(clonsFile)<1){
    stop("ERROR, clones.clns file not found")
  }
  if(file.exists(clonsFile)){
    assemble.ofile <- clonsFile
    st <- system2(command = "java",
                  args = c(software$mixcr$base.args,
                           "exportClones", "--chains TCR",
                           assemble.ofile, 
                           stringr::str_replace(assemble.ofile,".clns",".tab")
                  ), stderr = TRUE )
    if(any(stringr::str_detect(st, "ERROR"))){
      cat("\nSomething HAPPEND on clns")
    }
    # java -jar -Xmx6g $MIXCR_BIN_DIR/mixcr.jar exportClones --chains TCR -p min ${OUTPUT_FILE}.clones.clns ${OUTPUT_FILE}.vmin.clones.txt
    
    st <- system2(command = "java",
                  args = c(software$mixcr$base.args,
                           "exportClones", "--chains TCR",
                           "-p min",
                           assemble.ofile, 
                           stringr::str_replace(assemble.ofile,".clns",".vmin.tab")
                  ),  stderr = TRUE)
    if(any(stringr::str_detect(st, "ERROR"))){
      cat("\nSomething HAPPEND on vmin")
    }
    gen.files <- stringr::str_replace(assemble.ofile,".clns",".vmin.tab")
    # java -jar -Xmx6g $MIXCR_BIN_DIR/mixcr.jar exportClones -count -vHit -dHit -jHit -cHit -vHits -dHits -jHits -cHits -nFeature CDR3 -aaFeature CDR3 ${OUTPUT_FILE}.clones.clns ${OUTPUT_FILE}.vdetails.clones.txt
    system2(command = "java",
            args = c(software$mixcr$base.args,
                     "exportClones", "-count",
                     "-vhit", "-dhit", "-jHit", "-cHit", "-vHits", "-dHits", "-jHits", "-cHits",
                     "-nFeature CDR3", "-aaFeature CDR3",
                     assemble.ofile, 
                     stringr::str_replace(assemble.ofile,".clns",".vdetails.clones.tab")
            ), stderr = TRUE)
    gen.files <- c(gen.files,stringr::str_replace(assemble.ofile,".clns",".vdetails.clones.tab"))
    # java -jar -Xmx6g $MIXCR_BIN_DIR/mixcr.jar exportClones -count -vHit -dHit -jHit -cHit -vHits -dHits -jHits -cHits -nFeature CDR3 -c TRA -aaFeature CDR3 ${OUTPUT_FILE}.clones.clns ${OUTPUT_FILE}.TRA.clones.txt
    chains <- c("IGH","IGK","TRA","TRB","TRD","TRG")
    gen.files <- c(gen.files,unlist(lapply(chains, function(ch){
      sys.out <- system2(command = "java",
              args = c(software$mixcr$base.args,
                       "exportClones", "-count",
                       "-vHit", "-dHit", "-jHit", "-cHit", "-vHits", "-dHits", "-jHits", "-cHits",
                       "-nFeature CDR3", paste0("-c ",ch), "-aaFeature CDR3",
                       assemble.ofile, 
                       stringr::str_replace(assemble.ofile,".clns",paste0(".",ch,".clones.tab"))
              ),stderr = TRUE)
      ret.val <- stringr::str_replace(assemble.ofile,".clns",paste0(".",ch,".clones.tab"))
    
      attr(ret.val,"status") <- attr(sys.out,"status")
      attr(ret.val,"errmsg") <- attr(sys.out,"errmsg")
      return(ret.val)
    })))
    
    
    gen.files <- gen.files[file.exists(gen.files)]
    to.excel <- lapply(gen.files, function(x){
        return(read.table(x,h=T, sep="\t"))  
    })
    aux.files.subfix <- c(".tab","vmin.tab",".vdetails.clones.tab",".TRA.clones.tab",".TRB.clones.tab",".TRG.clones.tab",".TRD.clones.tab")
    names(to.excel) <- stringr::str_remove(c("clones.tab","vmin.tab","vdetails.clones.tab","TRA.clones.tab","TRB.clones.tab","TRG.clones.tab","TRD.clones.tab"),".tab")
    openxlsx::write.xlsx(to.excel, stringr::str_replace(assemble.ofile,".clns",".xlsx"))
    file.remove(c(gen.files,stringr::str_replace(assemble.ofile,".clns",".tab")))
    if(file.exists(stringr::str_replace(assemble.ofile,".clns",".xlsx"))){
      message(paste0("\nFile saved as :",stringr::str_replace(assemble.ofile,".clns",".xlsx")))
      return(invisible(aux.files.subfix))          
    }
    return(NA)
    
  }else{
    message(paste0("\nFile not found :",clonsFile))
    return(NA)
  }
  
}

#' RunMiXCRreport
#' @description 
#' In order to export alignment results or clones from a binary file (.clns) to a human-readable tab or excel file
#' RunMiXCRreport
#' @usage 
#' RunMiXCRreport(clonsFile, type = c("full","min"), otherFields, fileType = c("Excel","tab"))
#' @param clonsFile the full path to the clone file (subjsect.R1.fastq.clones.clns) as returned by \code{\link{RunMiXCRseqClones}}
#' @param type string "full" (default) or "min"
#' @param otherFields other fields (not implemented yet)
#' @details 
#' The resulting tab-delimited text file or excel file (type param) will contain columns with different types of information. 
#' If no options are specified (otherFields), the default set of columns - which is sufficient in most cases - will be exported. 
#' The possible columns include. See details in [MiXCR main page](https://mixcr.readthedocs.io/en/master/export.html): 
#' aligned sequences, qualities, all or just best hit for V, D, J and C genes, corresponding alignments, nucleotide and amino acid sequences of gene region present in sequence, etc. When exporting clones, the additional columns include: clone count, clone fraction etc.
#' One can customize the list of fields that will be exported by passing parameters to export commands (not implemented yet). 
#' 
#' @export
#'\dontrun{
#'      subj.file <- "/home/.../subject.R1.fastq"
#'       out.vdja.file <- RunMiXCR(sbj = subj.file)
#'       clns.file <- RNAseqClones(out.vdja.file)
#'       clone.excel.file <- RunMiXCRReport(clns.file, fileType = "Excel")
#'       require(openxlsx)
#'       CloneTable <- read.xlsx(clone.excel.file)
#'       View(CloneTable)
#'}
#' 
RunMiXCRreport <- function(clonsFile, type = c("full","min"), otherFields, fileType = c("Excel","tab")){
  software <- RMiXCR:::.OpenConfigFile()
  cat(paste0("\nExporting clones to Excel file"))
  # java -jar -Xmx6g $MIXCR_BIN_DIR/mixcr.jar exportClones --chains TCR ${OUTPUT_FILE}.clones.clns ${OUTPUT_FILE}.clones.txt
  clonsFile <- clonsFile[stringr::str_detect(clonsFile,".clones.clns")]
  if(length(clonsFile)<1){
    stop("ERROR, clones.clns file not found")
  }
  type <- match.arg(type[1], c("full","min"))
  fileType <- match.arg(toupper(fileType[1]) , toupper(c("Excel","tab")))
  assemble.ofile <- clonsFile
  out.tab.file <- stringr::str_replace(assemble.ofile,".clns",paste0("_",type,".tab"))
  st <- system2(command = "java",
                args = c(software$mixcr$base.args,
                         "exportClones", 
                         paste0("--preset ", type),
                         assemble.ofile, 
                         out.tab.file), 
                stderr = TRUE )
  if(file.exists(out.tab.file)==FALSE){
    stop("Fail in create clones file")
  }
  if(fileType == "TAB") return(out.tab.file)
  tab.file <- read.table(out.tab.file,h=T, sep="\t")
  openxlsx::write.xlsx(tab.file, file = stringr::str_replace(out.tab.file,".tab",".xlsx"))
  if(file.exists(stringr::str_replace(out.tab.file,".tab",".xlsx"))){
    file.remove(out.tab.file)
  }else{
    cat(paste0("\nExcel file failed"))
  }
  return(stringr::str_replace(out.tab.file,".tab","xlsx"))
}

