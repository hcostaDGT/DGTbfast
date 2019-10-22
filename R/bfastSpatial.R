#' Run BFAST on an entire raster grid
#' 
#' Apply \code{\link[bfast]{bfast}} on a entire raster (Landsat time series).
#' 
#' \strong{Overview}
#' 
#' \code{bfastSpatial} is theoretically designed to work on any generic raster time series, as long 
#' as a \code{dates} vector is provided. In the absence of a \code{dates} vector, \code{names(x)} 
#' should correspond exactly to respective Landsat scene ID's. In this case, 
#' \code{\link{Lmetadata}} is used to extract a dates vector, and subset by sensor if desired.
#' 
#' \strong{Arguments}
#' 
#' Arguments \code{start} and \code{end} adjust the period analysed. The time series is trimmed, 
#' butonly after interpolation (see argument \code{interpolate}).
#' 
#' Argument \code{interpolate} defines the method used to interpolate missing data. It
#' can be either 'linear' or 'periodic'. These methods use \code{\link[zoo]{na.approx}} and 
#' \code{\link[forecast]{na.interp}}, respectively (inspired in 
#' \href{https://philippgaertner.github.io/2018/04/bfast-preparation/}{this blog}).
#' 
#' Argument \code{aggregate} defines how the interpolated data (see argument \code{interpolate}),
#' which are daily time series, are aggregated to a different time frequency: weekly" , "biweekly"
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
#' fraction as expected by \code{\link[bfast]{bfast}}.
#' 
#' @param fast_bfast Logical. Enabling performance optimisations. This requires package \pkg{bfast} to be installed from a specific fork on \href{https://github.com/bfast2/bfast}{github}. If \code{FALSE}, the default implementation of \pkg{bfast} is used.
#' @inheritParams bfastPixel
#' @inheritParams bfast::bfast
#' 
#' @return 
#' A rasterBrick with layers depending on what has been supplied to \code{returnLayers}. 
#' By default, 3 layers are returned: 
#' (1) breakpoint: timing of breakpoints detected for each pixel; 
#' (2) magnitude: the median of the residuals within the monitoring period; 
#' (3) error: a value of 1 for pixels where an error was encountered by the algorithm (see \code{\link{try}}), and NA where the method was successfully run. See \code{\link{bfastmonitor}} for more information on the other possible layers.
#' 
#' @author Hugo Costa, Loic Dutrieux, Ben DeVries
#' 
#' @examples
#' 
#' \dontrun{
#' x<-list.files("C:/Landsat/Preprocessed/Landsat204032", "_preprocessed_NDVI.tif$",  full.names = TRUE, recursive = TRUE)
#' data(fire)
#' x<-fire
#' dates=NULL
#' #cell=8761217 # not used in this function
#' #start=c(2015,1)
#' #end=c(2018,365)
#' interpolate="periodic"
#' aggregate="weekly"
#' f=1
#' sensor=NULL
#' fast_bfast=TRUE
#' cores=2
#' #min.thresh=NULL # used in bfastPixel only
#' h=0.30952380952381
#' season = c('harmonic')
#' max.iter=10
#' breaks = NULL
#' hpc = 'none'
#' level=1
#' type= 'OLS-MOSUM'
#' 
#' 
#' data(fire)
#' system.time(
#'   out<-bfastSpatial(fire, interpolate="linear", h=0.30952380952381, season = c('harmonic'), level=1, max.iter=5, fast_bfast=TRUE, cores=2)
#' )
#' plot(out)
#' 
#' bfastSpatial(x, interpolate="periodic", cell=8761217, plot=TRUE, h=52, level=1, max.iter=5)
#' bfastSpatial(x, interpolate="periodic", aggregate="weekly", cell=8761217, plot=TRUE, h=0.30952380952381, level=1, max.iter=5)
#' bfastSpatial(x, interpolate="periodic", aggregate="monthly", cell=8761217, plot=TRUE, h=0.30952380952381, level=1, max.iter=5)
#' bfastSpatial(x, interpolate="periodic", aggregate="quarterly", cell=8761217, plot=TRUE, h=0.30952380952381, level=1, max.iter=5)
#' bfastSpatial(x, start=c(2015,1), end=c(2018,365), interpolate="periodic", cell=8761217, plot=TRUE, h=0.30952380952381, level=1, max.iter=5)
#' }
#' 
#' @seealso 
#' \code{\link{bfastPixel}}
#' 
#' @export
bfastSpatial <- function (x, dates=NULL, start=NULL, end=NULL, interpolate, aggregate="biweekly", f=1, sensor=NULL, fast_bfast=TRUE, cores,
                          h = 0.15, season = c('dummy', 'harmonic', 'none'), max.iter = NULL, 
                          breaks = NULL, hpc = 'none', level = 0.05, type= 'OLS-MOSUM'){
    
    # library(strucchange)
    
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
    
    # define cores
    DGTdiffNDVI:::.define_cores(cores)
    
    
    # import rasters, intersect extents and stack
    if(class(x)=="character"){
        r<-lapply(x, raster)
        e<-extent(r[[1]])
        for(i in 2:length(r)) e<-intersect(e, extent(r[[i]]))
        r<-lapply(r, crop, y=e)
        x<-stack(r)
    }
    
    
    
    # subset by sensor
    if(!is.null(sensor) & !DGTbfast:::.isLandsat(x)){
        warning("Scene info not available for subsetting based on sensor. Ignoring argument 'sensor'...")
        sensor <- NULL
    }
    if(!is.null(sensor)){
        if("ETM+" %in% sensor){
            sensor <- c(sensor, "ETM+ SLC-on", "ETM+ SLC-off")
        }
        # s <- Lmetadata(names(names))
        pixelts <- pixelts[which(info$Sensor %in% sensor)]
        s <- info[which(info$Sensor %in% sensor), ]
        dates <- s$Acquisition_date
    }
    
    # # optional: apply a threshold (if supplied)
    # if (!is.null(min.thresh))
    #     pixelts[pixelts <= min.thresh] <- NA
    
    # extract pixel time series JUST FOR TESTING
    # pixelts<- as.vector(x2[1])
    
    
    # define period of analysis
    startend<-format(as.Date(dates), "%Y%j")
    if(is.null(start)) start <- c(substr(min(startend), 1, 4), substr(min(startend), 5, nchar(min(startend))))
    if(is.null(end))   end <- c(substr(max(startend), 1, 4), substr(max(startend), 5, nchar(max(startend))))
    
    
    fun<-function(x)#, h2 = h, season2 = season, max.iter2 = max.iter, 
                  # breaks2 = breaks, hpc2 = hpc, level2 = level, type2 = type)
        {
        
        
        # optional: rescale values
        if(f != 1)
            x <- x * f
        
        # convert to a bfast ts object
        # TO DO: first, test if wheter pixelts is irregular or not
        pixelts <- bfastts(x, dates, type=c("irregular"))

        # interpolate
        if(interpolate=="linear"){
            s.d. <- round(zoo::na.approx(pixelts),0)
        }else if(interpolate=="periodic"){
            s.d. <- round(forecast::na.interp(pixelts),0)
        }

        # # optional: trim ts if monend is supplied
        # if(!is.null(monend))
        #     pixelts <- window(pixelts, end=monend)

        s.d. <- window(s.d., start=as.numeric(start))
        s.d. <- window(s.d., end=as.numeric(end))

        if(all(is.na(s.d.))){ #----------------------------- IF 1: time series is all NA
            
            t_magn <- t_time <- s_time <- NA
            # bkpt <- NA
            # magn <- NA
            # err <- NA
            # history <- NA
            # rsq <- NA
            # adj_rsq <- NA
            # coefficients <- rep(NA, coef_len)
            
            
        }else{
            
            # Aggregate Time Series
            s.f. <-.aggregate.time.series(s.d., aggregate=casefold(aggregate))
            
            # run bfast on the pixel time series
            if(h>=1){
                h<-h/length(s.f.)
                warning("Argument 'h' >= 1; converted to the relative fraction ", h, call. = FALSE, immediate. = TRUE)
            }
            
            if(fast_bfast){
                set_fast_options()
            }else{
                set_default_options()
            }
            bfast<-try(bfast(s.f., h = h, season = season, max.iter = max.iter,
                             breaks = breaks, hpc = hpc, level = level, type = type))
            
            
            
            # assign 1 to error and NA to all other fields if an error is encountered
            if(class(bfast) == "try-error") {
                t_magn <- t_time <- s_time <- -9999
                
            } else {
                # Yt,
                # output, bfast$output[[1]]$Tt,
                #                           St,
                #                           Nt,
                #                           Vt,
                #                           bp.Vt:List of 12,
                #                           Vt.bp,
                #                           ci.Vt:List of 5,
                #                           Wt,
                #                           bp.Wt:List of 12,
                #                           Wt.bp,
                #                           ci.Wt
                # nobp
                # Magnitude
                # Mags
                # Time
                # jump
                
                if(all(bfast$output[[1]]$Vt.bp==0)){
                    t_magn <- t_time <- NA
                }else{
                    t_magn   <-bfast$Magnitude                                 # magnitude of the biggest change detected in the trend component
                    t_time   <-bfast$Time                                      # timing of the biggest change detected in the trend component
                    # t_bpci2.5<-bfast$output[[1]]$ci.Vt$confint[,"2.5 %"]
                    # t_pb     <-bfast$output[[1]]$ci.Vt$confint[,"breakpoints"]
                    # t_bpci2.5<-bfast$output[[1]]$ci.Vt$confint[,"97.5 %"]
                    # t_mags   <-bfast$Mags[,3]
                    # s_bpci2.5<-bfast$output[[1]]$ci.Wt$confint[,"2.5 %"]
                    # s_pb     <-bfast$output[[1]]$ci.Wt$confint[,"breakpoints"]
                    # s_bpci2.5<-bfast$output[[1]]$ci.Wt$confint[,"97.5 %"]
                }
                
                if(all(bfast$output[[1]]$Wt.bp==0)){
                    s_time <- NA
                }else{
                    s_time   <-bfast$output[[1]]$Wt.bp
                    # t_bpci2.5<-bfast$output[[1]]$ci.Vt$confint[,"2.5 %"]
                    # t_pb     <-bfast$output[[1]]$ci.Vt$confint[,"breakpoints"]
                    # t_bpci2.5<-bfast$output[[1]]$ci.Vt$confint[,"97.5 %"]
                    # t_mags   <-bfast$Mags[,3]
                    # s_bpci2.5<-bfast$output[[1]]$ci.Wt$confint[,"2.5 %"]
                    # s_pb     <-bfast$output[[1]]$ci.Wt$confint[,"breakpoints"]
                    # s_bpci2.5<-bfast$output[[1]]$ci.Wt$confint[,"97.5 %"]
                }
                
                # bfast$output[[1]]$bp.Vt$datatsp
                
                # bkpt <- bfm$breakpoint
                # magn <- bfm$magnitude
                # err <- NA
                # history <- bfm$history[2] - bfm$history[1]
                # rsq <- summary(bfm$model)$r.squared
                # adj_rsq <- summary(bfm$model)$adj.r.squared
                # coefficients <- coef(bfm$model)
            }
            
        }         #----------------------------- END IF 1: time series is all NA

        res <- c(t_magn, t_time, s_time[1])
        # names(res) <- c("Mag Bp trend", "Time Bp trend", "Time Bp season")
        return(res)
        # return(x[1])
        
        # res <- c(bkpt, magn, err, history, rsq, adj_rsq)
        # names(res) <- c("breakpoint", "magnitude", "error", "history", "r.squared", "adj.r.squared")
        # res <- res[which(names(res) %in% returnLayers)]
        # if("coefficients" %in% returnLayers)
        #     res <- c(res, coefficients)
        # return(res)
        
        
    }

    # run function in parallel
    beginCluster(cores)
    out <- clusterR(x, calc, args=list(fun=fun))
    endCluster()
    names(out) <- c("Mag Bp trend", "Time Bp trend", "Time Bp season")
    
    return(out)
}
