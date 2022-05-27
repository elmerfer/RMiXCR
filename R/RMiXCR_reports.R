
#' Report
#' @description 
#' In order to export alignment results or clones from a binary file (.clns) to a human-readable tab or excel file
#' @details 
#' The resulting tab-delimited text file or excel file (type param) will contain columns with different types of information. 
#' If no options are specified (otherFields), the default set of columns - which is sufficient in most cases - will be exported. 
#' The possible columns include. See details in [MiXCR main page](https://mixcr.readthedocs.io/en/master/export.html): 
#' aligned sequences, qualities, all or just best hit for V, D, J and C genes, corresponding alignments, nucleotide and amino acid sequences of gene region present in sequence, etc. When exporting clones, the additional columns include: clone count, clone fraction etc.
#' One can customize the list of fields that will be exported by passing parameters to export commands (not implemented yet). 
#' @usage 
#' Report(clonsFile, type = c("full","min"), otherFields, fileType = c("Excel","tab"))
#' @param clonsFile the full path to the clone file (subjsect.R1.fastq.clones.clns) as returned by \code{\link{RunMiXCRseqClones}}
#' @param type string "full" (default) or "min"
#' @param otherFields other fields (not implemented yet)
#' @param fileType Excel or tab (tab delimited)
#' @export

Report <- function(clonsFile, type = c("full","min"), otherFields, fileType = c("Excel","tab")){
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

