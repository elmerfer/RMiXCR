
.CheckJava <- function(){
  message("\nChecking Java Version")
  jv <- try(system2(command = "java", args = "--version", stdout = TRUE))
  if(class(jv)=="try-error"){
    stop("java not installed , the version of java should be > 1.8")
  }
  jv <- as.numeric(unlist(stringr::str_split(unlist(stringr::str_split(jv," "))[2],"\\.")))
  if(jv[1]<1){
    stop("the version of java should be > 1.8")
  }
  if(jv[1]==1 & jv[2]<8){
    stop("the version of java should be > 1.8")
  }
  return(paste0(jv, collapse = "."))
}

#' InstallMiXCR
#' @param where (character) the path to the directory where I want to install MiXCR
#' @description The software MiXCR will be instaled at ./home/where/Software/mixcr-3.0.13/
#' if the directory "where" do not exists, the it will be created as well as directory "Software/mixcr-3.0.13/" below "where"
#' so it will look like
#' 
#' ./home/
#' 
#' ......where/
#' 
#' ...........Software/mixcr-3.0.13/ 
#' @export
#' @usage 
#' \dontrun{
#' InstallMiXCR("myfavoriteplace")
#' }
#' 

InstallMiXCR <- function(where){
  jv <- .CheckJava()
  message(paste0("\nJava version ",jv, "is OK"))
  software <- RMiXCR:::.OpenConfigFile()
  if(dir.exists(where)==FALSE){
    dir.create(where)
  }
  if(length(software)==0){
    if(dir.exists(file.path(where,"Software"))==FALSE){
      dir.create(file.path(where,"Software"))
    }
  }
  software$main <- file.path(where,"Software")
  
  cat("\nDownloading the MiXCR version 3.0.13")
  tmp.destfile <- tempfile()
  download.file(url = "https://github.com/milaboratory/mixcr/releases/download/v3.0.13/mixcr-3.0.13.zip",
                method = "wget",
                destfile = tmp.destfile)
  files <- zip::zip_list( tmp.destfile  )
  
  zip::unzip(tmp.destfile, exdir = software$main)
  software$mixcr$path <- file.path(software$main, files$filename[1])
  software$mixcr$command <- file.path(software$main, files$filename[stringr::str_detect(files$filename,"mixcr.jar")])
  software$mixcr$base.args <- c(paste0("-Xmx4g"),paste0("-Xms3g"),paste0("-jar ",software$mixcr$command))
  if(file.exists(software$mixcr$command)){
      system2(command = "java", args = software$mixcr$base.args)
  }
  RMiXCR:::.OpenConfigFile(software)
  file.remove(tmp.destfile)
}
