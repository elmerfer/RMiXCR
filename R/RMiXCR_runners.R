
.mixcrTest <- function(cmd){
  software <- RMiXCR:::.OpenConfigFile()
  system2(command= "java",
          args = c(software$mixcr$base.args, cmd))
  
}

#'RunMiXCR
#'
#'@param sbj full path to sample fastq R1 and R2. It should be named as xxxx_R1.fastq and xxxx_R2.fastq
#'@param species any of "hsa","mmu" (human - mouse)
#'@param assemblePartial (logical, dafault FALSE) if parial read should be assembled
#'@param extendedAlignments (logical, dafault FALSE) if TCR extension alignment should be performed
#'@param nThreads (integer) number of cores
#'@export
#'@return a list with clones, vmin, vdetails.clones, TRA.clones,TRB.clones,TRG.clones,TRD.clones data frames
#'it also save an Excel file as full.path.sbj.alignments.vdjca.xlsx
#'

RunMiXCR <- function(sbj, species = c("hsa","mmu"), assemblePartial=FALSE,extendedAlignments=FALSE, nThreads){
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
    for( i in 1:2){
      system2(command= "java",
              args = c(software$mixcr$base.args,
                       "assemblePartial",
                       apartial.input.file,
                       apartial.out.file))
      apartial.input.file <- apartial.out.file
    }
    if(extendedAlignments){
      cat("\nExtending TCR alñignments")
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

#' RunMiXCR.RNAseqAlignment
#'@param sbj full path to sample fastq R1 and R2. It should be named as xxxx_R1.fastq and xxxx_R2.fastq
#'@param species any of "hsa","mmu", "rat" (human - mouse, rat)
#'#'@param extendedAlignments (logical, default TRUE) if TCR extension alignment should be performed
#'@param nThreads (integer) number of cores (default 4, since it is suggested as the optimum value)
#'@export
#'@return the output file full path name.
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
#' It runs the MiXCR package for RNAseq experiments
#' @param sbj full path fastq(.gz) R1 sample
#' @param species character any of "hsa","mmu" (default hsa)
#' @export
#' @return if success, the it builds two files sbj.alignments.vdjca and sbj.clones.clns and return a vector with both full path file names
#' if fail, returns NA
RunMiXCRseq <- function(sbj, species = c("hsa","mmu"), extendedAlignments = TRUE ){
  
  out.alignment <- RunMiXCRseqAlignment(sbj = sbj, species = species, extendedAlignments = extendedAlignments)
  out.clones <- RunMiXCRseqClones(alignmentFile = out.alignment)
  files <- c(out.alignment, out.clones)
  if(all(file.exists(files))==TRUE){
    message(paste0("\n files created",basename(files) ))
    return(files)
  }
  return(NA)
  
  
}

#' RunMiXCR.RNAseqClones
#'@param alignmentFile full path to sample fastq R1 and R2. It should be named as xxxx_R1.fastq and xxxx_R2.fastq
#'@param nThreads (integer) number of cores (4 default as it is suggested as the optimal value)
#'@export
#'@return the full path name of the clones file
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
#' @param clonsFile the full path to the clone file
#' @param type string "full" (default) or "min"
#' @param otherFields other fileds (not implemented yet)
#' @export
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
