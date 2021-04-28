#' .OpenConfigFile
#' Internal. save and restore config info
#' @param config config list structure with software$softName...
#'
.OpenConfigFile <- function(config){
  config.file <- file.path( .libPaths()[1],"RMiXCR.RData")
  if(file.exists(config.file)){
    if(missing(config)){##si lo llama vacio, devuelve lo que hay en el file
      config <- readRDS(config.file)
      return(invisible(config))
    }else{
      saveRDS(config,config.file)
      return(invisible(config))
    }
  }else{
    saveRDS(list(),config.file)
    return(list())
  }
  
}