#' Retrieve Landsat info from filenames
#' 
#' Retrieves information on sensor, pathrow, date, collection and tier from the 
#' \href{https://www.usgs.gov/faqs/what-naming-convention-landsat-collections-level-1-scenes?qt-news_science_products=0#qt-news_science_products}{Landsat Collections Level-1 scenes naming convention}.
#' 
#' @param x Character vector or Raster* objec. If Raster* object, its \code{\link[raster]{names}} will be used.
#' @param ... Additional arguments to pass to \code{\link[utils]{write.csv}}.
#' 
#' @author Hugo Costa
#' 
#' @return \code{data.frame}
#' 
#' @examples 
#' \dontrun{
#' sourcefile<-c("LC08_L1TP_204032_20190925_20190925_01_RT_sr_band1", "LC08_L1TP_204032_20160425_20170326_01_T1_sr_band1.tif")
#' Lmetadata(sourcefile)
#' }
#' 
#' @export 
Lmetadata <- function(x, ...){
    
    if(!.isLandsat(x)) stop("Scene info not found.")
    
    if(inherits(x, 'RasterStackBrick')){
        x<-names(x)
    }
    
    # get scene info
    # What is the naming convention for Landsat Collections Level-1 scenes?
    # https://www.usgs.gov/faqs/what-naming-convention-landsat-collections-level-1-scenes?qt-news_science_products=0#qt-news_science_products
    x<-basename(x)
    y<-strsplit(x,"_")
    
    LXSS<-sapply(y, '[[', 1)      # Sensor
    LLLL<-sapply(y, '[[', 2)      # Processing correction level (L1TP/L1GT/L1GS)
    PPPRRR<-sapply(y, '[[', 3)    # WRS path row
    YYYYMMDD<-sapply(y, '[[', 4)  # Acquisition date
    yyyymmdd<-sapply(y, '[[', 5)  # Processing date
    CC<-sapply(y, '[[', 6)        # Collection number
    TX<-sapply(y, '[[', 7)        # Collection category
    
    info <- data.frame(Sensor = LXSS, 
                       Pathrow = PPPRRR, 
                       Acquisition_date = as.Date(YYYYMMDD, format="%Y%m%d"), 
                       Processing_date = as.Date(yyyymmdd, format="%Y%m%d"),
                       Collection = CC,
                       Tier = TX)
    rownames(info)<-substr(x, 1, 40)
  
    
    # optional: print to .csv for future reference
    if(methods::hasArg(file)) 
        utils::write.csv(info, ...)
  
    return(info)
}
