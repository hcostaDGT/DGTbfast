#' Apply BFAST on a single pixel
#' 
#' Apply \code{\link[bfast]{bfast}} on a single pixel of known cell index or coordinates, or 
#' interactively by clicking on a plot.
#' 
#' \strong{Overview}
#' 
#' \code{bfastPixel} is theoretically designed to work on any generic raster time series, as long 
#' as a \code{dates} vector is provided. In the absence of a \code{dates} vector, \code{x} or 
#' \code{names(x)} should correspond to respective Landsat filenames (see the
#' \href{https://www.usgs.gov/faqs/what-naming-convention-landsat-collections-level-1-scenes?qt-news_science_products=0#qt-news_science_products}{Landsat Collections Level-1 scenes naming convention}).
#' In this case, \code{\link{Lmetadata}} is used to extract a dates vector.
#' 
#' \strong{Arguments}
#' 
#' Arguments \code{start} and \code{end} adjust the period analysed. The time series is trimmed, 
#' but only after interpolation (see argument \code{interpolate}).
#' 
#' Argument \code{interpolate} defines the method used to interpolate missing data. It
#' can be either 'linear' or 'periodic'. These methods use \code{\link[zoo]{na.approx}} and 
#' \code{\link[forecast]{na.interp}}, respectively (inspired in 
#' \href{https://philippgaertner.github.io/2018/04/bfast-preparation/}{this blog}).
#' 
#' Argument \code{aggregate} defines how the interpolated data (see argument \code{interpolate}),
#' which are daily time series, are aggregated to a different time frequency: "weekly" , "biweekly"
#' (also "fortnightly"), "monthly", "quarterly", or "yearly" (also "annually"). The funtions used
#' for aggregate were copied from 
#' \href{https://philippgaertner.github.io/2018/04/bfast-preparation/}{this same blog}).
#' 
#' Arguments \code{h}, \code{season}, \code{max.iter}, \code{breaks}, \code{hpc}, \code{level}, and
#' \code{type} are inherited from \code{\link[bfast]{bfast}}. Note, however, that although argument
#' \code{h} is expected to be <= 0.5 in \code{\link[bfast]{bfast}}, values >=1 are accepted here. 
#' When so, \code{h} is interpreted as the minimal segment size in \strong{absolute} value. For 
#' example, in a monthly analysis, \code{h=12} means that the minimal segment size is 12 months.
#' Then, \code{h} is divided by the total length of the time series so to produce a relative 
#' fraction as expected by \code{\link[bfast]{bfast}}. This relative fraction is printed out as a 
#' warning.
#' 
#' @param x Character vector. Full path names to images or to parent folders containing images.
#' @param start (Optional) Numeric or character. The start of the period of interest in the format \code{c(year, julian day)}.
#' @param end (Optional) Numeric or character. The end of the period of interest in the format \code{c(year, julian day)}.
#' @param interpolate Character. Either "linear" or "periodic" (see details).
#' @param aggregate Character. Either "weekly" , "biweekly" (also "fortnightly"), "monthly", "quarterly", or "yearly" (also "annually").
#' @param cell Numeric. Either (i) a numeric of length 1 indicating the cell ID to be observed, or (ii) a numeric of length 2 representing the (x,y) coordinate of the raster cell to be observed. If \code{NULL}, argument \code{interactive} must be set to \code{TRUE}.
#' @param interactive Logical. Select cell by clicking on an already plotted map? (see \code{\link[raster]{click}}). If \code{FALSE}, \code{cell} cannot be \code{NULL}.
#' @param f Numeric. Factor by which to rescale values before running \code{bfast}. Defaults to 1 (no rescaling)
#' @param plot Logical. Plot the result? Defaults to \code{FALSE}.
#' @inheritParams bfast::bfast
#' 
#' @return A list with the following components:
#' \enumerate{
#'   \item $bfast - an object of class 'bfast' (see \code{\link[bfast]{bfast}})
#'   \item $cell - the cell index (an integer of length 1). This can be used to run 
#'   \code{bfastPixel} again on the same pixel (with different parameters) without having to click
#'   on a plot again to find the same pixel (in that case, be sure to set 
#'   \code{interactive=FALSE} for subsequent trials!).
#' }
#' 
#' @author Hugo Costa and Ben DeVries
#' 
#' @examples
#' \dontrun{
#' list.files("//DGT-692/Dgt-692- externo/Landsat8/Preprocessed/204032", "_preprocessed_NDVI.tif$",  full.names = TRUE, recursive = TRUE)
#' data(fire3)
#' 
#' bfastPixel(fire3, interpolate="linear", cell=8761217, plot=TRUE, h=52, level=1, max.iter=5)
#' bfastPixel(fire3, interpolate="periodic", cell=8761217, plot=TRUE, h=52, level=1, max.iter=5)
#' bfastPixel(fire3, interpolate="periodic", aggregate="weekly", cell=8761217, plot=TRUE, h=0.30952380952381, level=1, max.iter=5)
#' bfastPixel(fire3, interpolate="periodic", aggregate="monthly", cell=8761217, plot=TRUE, h=0.30952380952381, level=1, max.iter=5)
#' bfastPixel(fire3, interpolate="periodic", aggregate="quarterly", cell=8761217, plot=TRUE, h=0.30952380952381, level=1, max.iter=5)
#' bfastPixel(fire3, start=c(2015,1), end=c(2018,365), interpolate="periodic", cell=8761217, plot=TRUE, h=0.30952380952381, level=1, max.iter=5)
#' }
#' 
#' @seealso 
#' \code{\link{bfastRaster}}
#' @export

bfastPixel <- function (x, dates=NULL, start=NULL, end=NULL, interpolate, aggregate="biweekly", cell=NULL, interactive=FALSE, f=1, plot=TRUE,
                        h = 0.15, season = c('dummy', 'harmonic', 'none'), max.iter = NULL, 
                        breaks = NULL, hpc = 'none', level = 0.05, type= 'OLS-MOSUM'){
    
    # check arguments
    stopifnot((is.null(cell) & interactive) | (!is.null(cell) & !interactive))
    
    
    # get dates
    if(is.null(dates)) {
        info<-Lmetadata(x)
        dates <- info$Acquisition_date
        
    }else{
        if(inherits(x, 'RasterStackBrick')){
            stopifnot(length(dates)==nlayers(x))
        }else{
            stopifnot(length(dates)==length(x))
        }
        
    }
    
    
    # import rasters, intersect extents and stack
    if(class(x)=="character"){
        r<-lapply(x, raster::raster)
        e<-raster::extent(r[[1]])
        for(i in 2:length(r)) e<-raster::intersect(e, extent(r[[i]]))
        r<-lapply(r, raster::crop, y=e)
        x<-raster::stack(r)
    }
    
    # select pixel
    if (interactive) {
        cell <- as.data.frame(raster::click(x, n=1, id=TRUE, cell=TRUE, show=FALSE))$cell
    } else {
        cell <- ifelse(length(cell)==2, raster::cellFromXY(x, t(as.matrix(cell))), cell)
    }
    
    
    # extract pixel time series
    pixelts <- as.vector(x[cell])
    
    
    # rescale values
    if(f != 1) pixelts <- pixelts * f
    
    # convert time series to to a bfast ts object
    pixelts <- bfastts(pixelts, dates, type=c("irregular"))
    
    
    # interpolate time series
    if(interpolate=="linear")   s.d. <- round(zoo::na.approx(pixelts),0)
    if(interpolate=="periodic") s.d. <- round(forecast::na.interp(pixelts),0)
    
    
    # redefine period of analysis
    if(!is.null(start)) s.d. <- window(s.d., start=as.numeric(start))
    if(!is.null(end))   s.d. <- window(s.d., end=as.numeric(end))

    
    
    if(!all(is.na(s.d.))){
        
        # Aggregate Time Series
        s.f. <-.aggregate.time.series(s.d., aggregate=casefold(aggregate))
        
        # run bfast
        if(h>=1){
            h<-h/length(s.f.)
            warning("Argument 'h' >= 1; converted to the relative fraction ", h, call. = FALSE, immediate. = TRUE)
        }
        
        bfast<-bfast(s.f., h = h, season = season, max.iter = max.iter,
                         breaks = breaks, hpc = hpc, level = level, type = type)
        
        
        # return
        if(plot) plot(bfast)
        return(list(bfast=bfast, cell=cell))
        
    }else{
        warning("The pixel selected has no values in the time series. Pick a different pixel.")
    }
}
