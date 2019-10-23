#' Internal functions
#'
#' Functions to aggregate daily time series (copied from
#' \href{https://philippgaertner.github.io/2018/04/bfast-preparation/}{here})
#'
#' @keywords internal
#'
#' @usage
#' .aggregate.time.series(x, aggregate)
#'
#' @param x Regular time series.
#' @param aggregate Character. See \code{\link{bfastPixel}}.
#'
#' @return Aggregated time series
#'
#' @author Hugo Costa
#'
#' @name .DGTaggregate
NULL

#' @rdname .DGTaggregate
.aggregate.time.series <- function(x, aggregate) {
    
    if(casefold(aggregate)=="weekly"){
        x   <- DGTbfast:::.aggregate.daily.to.weekly(x)
    }else if(casefold(aggregate)=="fortnightly" | casefold(aggregate)=="biweekly"){
        x   <- DGTbfast:::.aggregate.daily.to.fortnight(x)
    }else if(casefold(aggregate)=="monthly"){
        x   <- DGTbfast:::.aggregate.daily.to.monthly(x)
    }else if(casefold(aggregate)=="quarterly"){
        x   <- DGTbfast:::.aggregate.daily.to.quarterly(x)
    }else if(casefold(aggregate)=="yearly" | casefold(aggregate)=="annually"){
        x   <- DGTbfast:::.aggregate.daily.to.yearly(x)
    }else{
        stop("Argument 'aggregate' wrongly defined.", call. = FALSE)
    }
    x
}

#' @rdname .DGTaggregate
.aggregate.daily.to.weekly <- function(daily.ts) {

  dates      <- as.Date(lubridate::date_decimal(as.numeric(time(daily.ts))))

  xts.daily  <- xts::xts(daily.ts, order.by = dates)

  xts.weekly <- round(xts::apply.weekly(xts.daily, median),0)  # xts

  start(xts.weekly)
  ts.weekly <- ts(data = as.numeric(xts.weekly),
                  # define the start and end (Year, Week)
                  start = c(as.numeric(format(start(xts.weekly), "%Y")),
                            as.numeric(format(start(xts.weekly), "%W"))),
                  end   = c(as.numeric(format(end(xts.weekly), "%Y")),
                            as.numeric(format(end(xts.weekly), "%W"))),
                  frequency = 52)

  return(ts.weekly) # original code from the blog
  # return(ts.weekly[1:length(ts.weekly)])
}

#' @rdname .DGTaggregate
.aggregate.daily.to.fortnight <- function(daily.ts) {

  daily.ts.df           <- data.frame(as.numeric(daily.ts))
  colnames(daily.ts.df) <- "NDVI"
  daily.ts.df$date      <- as.Date(lubridate::date_decimal(as.numeric(time(daily.ts))))


  ## bi-weekly / fortnightly averages with function timeAverage {openair}
  fortnight      <- openair::timeAverage(daily.ts.df, avg.time = "2 week")
  fortnight$time <- as.Date(fortnight$date)

  fortnight.ts <- ts(fortnight$NDVI,
                     # define the start and end (Year, Week)
                     start = c(as.numeric(format(min(fortnight$date), "%Y")),
                               round(lubridate::yday(min(fortnight$date)) / 14,0)), # DOY / 14
                     end   = c(as.numeric(format(max(fortnight$date), "%Y")),
                               round(lubridate::yday(max(fortnight$date)) / 14,0)), # DOY / 14
                     frequency = 26)

  return(fortnight.ts)
}

#' @rdname .DGTaggregate
.aggregate.daily.to.monthly <- function(daily.ts) {
    
    s.month <- round(aggregate(as.zoo(daily.ts), as.yearmon, median), 0)   # zoo
    s.month <- as.ts(s.month)
    
    return(s.month)
    
}

#' @rdname .DGTaggregate
.aggregate.daily.to.quarterly <- function(daily.ts) {
    
    s.qtr <- round(aggregate(zoo::as.zoo(daily.ts), as.yearqtr, median), 0)   # zoo
    s.qtr <- as.ts(s.qtr)
    
    return(s.qtr)
    
}

#' @rdname .DGTaggregate
.aggregate.daily.to.yearly <- function(daily.ts, interpolate) {
    
    # function to create dataframe
    time.series.to.dataframe <- function(time_series, interpolate) {
        
        s.df           <- data.frame(as.numeric(time_series))
        colnames(s.df) <- "NDVI"
        s.df$Time      <- as.Date(lubridate::date_decimal(as.numeric(time(time_series))))
        s.df$Type      <- interpolate
        return(s.df)
        
    }
    
    # create data.frame
    daily.df <- time.series.to.dataframe(daily.ts, "not needed here")
    
    # add year
    daily.df$year <- format(daily.df$Time, "%Y")
    
    # add season
    daily.df$seas <- seas::mkseas(x = daily.df, width = "DJF")
    
    # calculate median per season within each year
    year.df <- aggregate(NDVI ~ seas + year, data = daily.df, median)
    
    # based on variable values
    year.df.JJA <- year.df[ which(year.df$seas == 'JJA'), ]
    
    # create ts
    year.ts.JJA <- ts(year.df.JJA$NDVI, 
                      start = c(as.numeric(min(year.df.JJA$year)), 1), # freq 1
                      end   = c(as.numeric(max(year.df.JJA$year)), 1), # freq 1
                      frequency = 1)
    
    return(year.ts.JJA)
}


