# ##check on consolo the java version
# java.version <- system2(command = "java", args = "--version", stdout = TRUE)
# java.version <- unlist(stringr::str_split(java.version," "))[2]
# 
# mixcr.path <- "/media/respaldo4t/mixcr3/mixcr.jar"
# system2(command= "java",
#         args = c(paste0("-Xmx4g"),
#                  paste0("-Xms3g"),
#                  paste0("-jar ",mixcr.path),"analyze","--help"))
# 
# out <- system2(command= "java",
#         args = c(paste0("-Xmx4g"),
#                  paste0("-Xms3g"),
#                  paste0("-jar ",mixcr.path),
#                  # "-v",
#                  "analyze",
#                  "shotgun",
#                  "--species mmu",
#                  "--starting-material rna",
#                  "--only-productive",
#                  "/media/respaldo4t/Consultorias/Montes/mouse/D120T33.R1.fastq.gz",
#                  "/media/respaldo4t/Consultorias/Montes/mouse/D120T33.R2.fastq.gz",
#                  "analysis",
#                  "--export "
#                  # 
#                  
#                  ),
#         stdout = TRUE
# )
# 
# destdir <- "/media/respaldo4t/Consultorias/Montes/mouse/"
# sample <- "D120T33"
# analysis.files <- list.files(getwd(),full.names = TRUE)
# analysis.files <- analysis.files[stringr::str_detect(analysis.files,"analysis.")]
# destination.files <- file.path(destdir, paste0(sample,"_",basename(analysis.files)))
# file.copy(analysis.files,file.path(destdir, paste0(sample,"_",basename(analysis.files))))
# file.remove(analysis.files)
# 
# 
# system2(command= "java",
#                args = c(paste0("-Xmx4g"),
#                         paste0("-Xms3g"),
#                         paste0("-jar ",mixcr.path), "analyze","help" ))
#                         
# targets <- openxlsx::read.xlsx("/media/respaldo4t/Consultorias/Montes/mouse/TablaURLS.xlsx")
# log.file <- "/media/respaldo4t/Consultorias/Montes/mouse/TablaURLS.xlsx"
# lapply(1:nrow(targets), function(id){
#   fileID1 <- rev(unlist(stringr::str_split(targets$url1[id],"/")))[2]
#   fileID2 <- rev(unlist(stringr::str_split(targets$url2[id],"/")))[2]
#   fileOUT <- targets$Name[id]
#   cat(paste0("\n-----------------" ))
#   cat(paste0("\n",fileID1, ": R1" ))
#   arg2 <- paste0("wget --load-cookies /tmp/cookies.txt \"https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=FILEID' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\\1\\n/p'",")","&id=FILEID\"")
#   arg2 <- stringr::str_replace_all(arg2, "FILEID",fileID1)
#   arg2 <- paste(arg2," -O ",paste0("/media/respaldo4t/Consultorias/Montes/mouse/","D120T",fileOUT,".R1.fastq.gz")," && rm -rf /tmp/cookies.txt")
#   system(command = arg2, ignore.stdout = TRUE)
#   cat(paste0("\n",fileID1, ": R2" ))
#   
#   arg2 <- paste0("wget --load-cookies /tmp/cookies.txt \"https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=FILEID' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\\1\\n/p'",")","&id=FILEID\"")
#   arg2 <- stringr::str_replace_all(arg2, "FILEID",fileID2)
#   arg2 <- paste(arg2," -O ",paste0("/media/respaldo4t/Consultorias/Montes/mouse/","D120T",fileOUT,".R2.fastq.gz")," && rm -rf /tmp/cookies.txt")
#   system(command = arg2, ignore.stdout = TRUE)
# 
# })


mixcr.path <- "/media/respaldo4t/mixcr3/mixcr.jar"
# java -jar -Xmx6g $MIXCR_BIN_DIR/mixcr.jar align --species hsa --chains TCR --parameters rna-seq -OallowPartialAlignments=true --threads 8 ${INPUT_FASTQ}.R1.fastq.gz ${INPUT_FASTQ}.R2.fastq.gz ${OUTPUT_FILE}.alignments.vdjca
species <- "mmu"#hsa for human
nThreads <- 6
file1 <- "/media/respaldo4t/Consultorias/Montes/mouse/D120T33.R1.fastq.gz"
file2 <- stringr::str_replace(file1,".R1.fasta",".R2.fasta")
align.ofile <- file.path(dirname(file1), paste0(basename(file1),".alignments.vdjca"))
base.args <- c(paste0("-Xmx4g"),paste0("-Xms3g"),paste0("-jar ",mixcr.path))

RunMiXCR <- function(sbj, species = c("hsa","mmu"), nTreads){
      software <- RMiXCR:::.OpenConfigFile()  
      cat(paste0("\nAlignment"))
      if(missing(nTreads)){
        nt <- parallel::detectCores()
        nTreads <- ifelse(nt-2 < 1, 1, nt)
        message(paste0("\nRunning with ", nt, " Cores"))
      }
      species <- match.arg(species, c("hsa","mmu") )
      file1 <- sbj
      file2 <- stringr::str_replace(file1,".R1.fasta",".R2.fasta")
      align.ofile <- file.path(dirname(file1), paste0(basename(file1),".alignments.vdjca"))
      system2(command= "java",
              args = c(software$mixcr$base.args,"align",
                       paste0("--species ",species),
                       "--chains TCR",
                       "--parameters rna-seq",
                       "-OallowPartialAlignments=true",
                       paste0("--threads ",nThreads),
                       file1,
                       file2,
                       align.ofile))
      if(file.exists(align.ofile)==FALSE){
        stop(align.ofile)
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
cat(paste0("\nAlignment"))
system2(command= "java",
        args = c(base.args,"align",
                 paste0("--species ",species),
                 "--chains TCR",
                 "--parameters rna-seq",
                 "-OallowPartialAlignments=true",
                 paste0("--threads ",nThreads),
                 file1,
                 file2,
                 align.ofile))
if(file.exists(align.ofile)==FALSE){
  stop(align.ofile)
}
# java -jar -Xmx6g $MIXCR_BIN_DIR/mixcr.jar assemble --threads 8 -OmaxBadPointsPercent=0 ${OUTPUT_FILE}.alignments.vdjca ${OUTPUT_FILE}.clones.clns
assemble.ofile <- file.path(dirname(file1), paste0(basename(file1),".clones.clns"))
cat(paste0("\nAssembling"))
system2(command = "java",
        args = c(base.args,
                 "assemble",
                 paste0("--threads ",nThreads),
                 "-OmaxBadPointsPercent=0",
                 align.ofile,
                 assemble.ofile
                 ))
if(file.exists(assemble.ofile)==FALSE){
  stop(assemble.ofile)
}
#Exporting clones to tab-delimited file
cat(paste0("\nExporting clones to tab-delimited file"))
# java -jar -Xmx6g $MIXCR_BIN_DIR/mixcr.jar exportClones --chains TCR ${OUTPUT_FILE}.clones.clns ${OUTPUT_FILE}.clones.txt
system2(command = "java",
        args = c(base.args,
                 "exportClones", "--chains TCR",
                 assemble.ofile, 
                 stringr::str_replace(assemble.ofile,".clns",".tab")
                 ))
# java -jar -Xmx6g $MIXCR_BIN_DIR/mixcr.jar exportClones --chains TCR -p min ${OUTPUT_FILE}.clones.clns ${OUTPUT_FILE}.vmin.clones.txt
system2(command = "java",
        args = c(base.args,
                 "exportClones", "--chains TCR",
                 "-p min",
                 assemble.ofile, 
                 stringr::str_replace(assemble.ofile,".clns","vmin.tab")
        ))
# java -jar -Xmx6g $MIXCR_BIN_DIR/mixcr.jar exportClones -count -vHit -dHit -jHit -cHit -vHits -dHits -jHits -cHits -nFeature CDR3 -aaFeature CDR3 ${OUTPUT_FILE}.clones.clns ${OUTPUT_FILE}.vdetails.clones.txt
system2(command = "java",
        args = c(base.args,
                 "exportClones", "-count",
                 "-vhit", "-dhit", "-jHit", "-cHit", "-vHits", "-dHits", "-jHits", "-cHits",
                 "-nFeature CDR3", "-aaFeature CDR3",
                 assemble.ofile, 
                 stringr::str_replace(assemble.ofile,".clns",".vdetails.clones.tab")
        ))
# java -jar -Xmx6g $MIXCR_BIN_DIR/mixcr.jar exportClones -count -vHit -dHit -jHit -cHit -vHits -dHits -jHits -cHits -nFeature CDR3 -c TRA -aaFeature CDR3 ${OUTPUT_FILE}.clones.clns ${OUTPUT_FILE}.TRA.clones.txt
system2(command = "java",
        args = c(base.args,
                 "exportClones", "-count",
                 "-vHit", "-dHit", "-jHit", "-cHit", "-vHits", "-dHits", "-jHits", "-cHits",
                 "-nFeature CDR3", "-c TRA", "-aaFeature CDR3",
                 assemble.ofile, 
                 stringr::str_replace(assemble.ofile,".clns",".TRA.clones.tab")
        ))
                 
# java -jar -Xmx6g $MIXCR_BIN_DIR/mixcr.jar exportClones -count -vHit -dHit -jHit -cHit -vHits -dHits -jHits -cHits -nFeature CDR3 -c TRB -aaFeature CDR3 ${OUTPUT_FILE}.clones.clns ${OUTPUT_FILE}.TRB.clones.txt
system2(command = "java",
        args = c(base.args,
                 "exportClones", "-count",
                 "-vHit", "-dHit", "-jHit", "-cHit", "-vHits", "-dHits", "-jHits", "-cHits",
                 "-nFeature CDR3", "-c TRB", "-aaFeature CDR3",
                 assemble.ofile, 
                 stringr::str_replace(assemble.ofile,".clns",".TRB.clones.tab")
        ))
# java -jar -Xmx6g $MIXCR_BIN_DIR/mixcr.jar exportClones -count -vHit -dHit -jHit -cHit -vHits -dHits -jHits -cHits -nFeature CDR3 -c TRG -aaFeature CDR3 ${OUTPUT_FILE}.clones.clns ${OUTPUT_FILE}.TRG.clones.txt
system2(command = "java",
        args = c(base.args,
                 "exportClones", "-count",
                 "-vHit", "-dHit", "-jHit", "-cHit", "-vHits", "-dHits", "-jHits", "-cHits",
                 "-nFeature CDR3", "-c TRG", "-aaFeature CDR3",
                 assemble.ofile, 
                 stringr::str_replace(assemble.ofile,".clns",".TRG.clones.tab")
        ))
# java -jar -Xmx6g $MIXCR_BIN_DIR/mixcr.jar exportClones -count -vHit -dHit -jHit -cHit -vHits -dHits -jHits -cHits -nFeature CDR3 -c TRD -aaFeature CDR3 ${OUTPUT_FILE}.clones.clns ${OUTPUT_FILE}.TRD.clones.txt
system2(command = "java",
        args = c(base.args,
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


                     