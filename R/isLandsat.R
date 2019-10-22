#' @title Check if a String is a Landsat file
#' 
#' @description Parses a vector of characters to check if they conform to Landsat naming convention
#'  
#' @param x Character 
#' 
#' @return \code{TRUE} if all elements of \code{x} match Landsat scene ID criteria, or \code{FALSE} otherwise.
#' 
#' @author Hugo Costa
#' 
#' @examples 
#' sourcefile<-c("LC08_L1TP_204032_20190925_20190925_01_RT_sr_band1", "LC08_L1TP_204032_20160425_20170326_01_T1_sr_band1.tif")
#' DGTbfast:::.isLandsat(sourcefile)
#' 
#' sourcefile<-c("LC08_L1TP_204032_20190925_20190925_01_RT_sr_band1", "LC08_L1TP_204032")
#' DGTbfast:::.isLandsat(sourcefile)
#' 
.isLandsat<- function(x){
    
    if(inherits(x, 'RasterStackBrick')){
        x<-names(x)
    }
    
    x<-basename(x)
    
    test1<-all(grepl(pattern='^(LT4|LT5|LE7|LC08)', x))
    test2<-all(nchar(x)>=40)
    test3<-all(sapply(strsplit(x,"_"), length)>=7)
    
    test1 & test2 & test3
    
}
