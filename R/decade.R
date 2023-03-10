## Multiple plot function

multiplot <- function(..., plotlist=NULL, cols=1, layout=NULL) {
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots <- length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots == 1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), 
                                                                 ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
    }
  }
}


###############################################################################


#' Plot decimals
#' 
#' Plot year-by-year distribution of the decimals
#' in order to investigate the actual reporting resolution.
#' 
#' @param Data A character string giving the path of the input file, or
#' a 5 or 7-column matrix (depending on data type) with following columns: 
#' variable code, year, month, day, (hour), (minute), value. 
#' @param outfile Character string giving the path of the output pdf file.
#' @param startyear First year to plot. If not indicated, all available years
#' until \code{endyear} will be plotted.
#' @param endyear Last year to plot. If not indicated, all available years
#' since \code{startyear} will be plotted.
#'
#' @details
#' The input file must follow the C3S Station Exchange Format (SEF).
#' 
#' Only the first digit after the decimal point is analysed. If there is more
#' than one digit, the data will be rounded to the first decimal place.
#' 
#' @note 
#' For precipitation and other bounded variables one needs to remove the values
#' at the boundaries from the input data (e.g., zeros for precipitation).
#' 
#' @author Stefan Hunziker, Yuri Brugnara
#' 
#' @references
#' Hunziker et al., 2017: Identifying, attributing, and overcoming common data quality 
#' issues of manned station observations. Int. J. Climatol, 37: 4131-4145.
#' 
#' Hunziker et al., 2018: Effects of undetected data quality issues on climatological
#' analyses. Clim. Past, 14: 1-20.
#' 
#' @examples 
#' plot_decimals(Rosario$Tx, outfile = paste0(tempdir(),"/test.pdf"))
#' 
#' @import grDevices
#' @export

plot_decimals <- function(Data, outfile, startyear = NA, endyear = NA) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("Package \"grid\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  if (is.null(dim(Data))) {
    Data <- read_sef(Data)
    if (dim(Data)[2] == 7) Data <- Data[, c(1:4,7)]
  }
  
  if (!is.na(startyear)) {
    Data <- Data[which(Data[,2] >= startyear), ]
  }
  
  if (!is.na(endyear)) {
    Data <- Data[which(Data[,2] <= endyear), ]
  }
  
  if (substr(outfile, nchar(outfile)-2, nchar(outfile)) != "pdf") {
    outfile <- paste(outfile, "pdf", sep = ".")
  }
  
  dec.labels <- c("x.0","x.1","x.2","x.3","x.4","x.5","x.6","x.7","x.8","x.9")
  
  ## Remove missing values
  Data <- Data[which(!is.na(Data[,dim(Data)[2]])), ]
  if (nrow(Data) == 0) warning("No valid data found in the selected time interval")
  
  ## Get decimals in precision of 0.1
  data.vec <- round(Data[,dim(Data)[2]], 1)
  decimal.vec <- as.integer(10 * abs(data.vec - trunc(data.vec))) 
  
  ## Define data frame with date and decimals
  decimal.frame <- data.frame(Data[,2], Data[,3], Data[,4], decimal.vec)
  
  ## Count decimals per year
  start <- min(decimal.frame[,1])
  end <- max(decimal.frame[,1])
  year.vec <- start:end
  decimal.frame.count <- data.frame(matrix(rep(year.vec,11), ncol=11))
  colnames(decimal.frame.count) <- c("year", dec.labels)  
  for (year in start:end) {
    for (i in 0:9) {
      decimal.frame.count[year-start+1, i+2] <- 
        sum(decimal.frame[which(decimal.frame[,1]==year),4] == i)
    }
  }
  
  ## Replace rows (i.e. years) without values by NAs
  for (i in 1:length(year.vec)){
    if (sum(decimal.frame.count[i,2:11]) <= 0) {
      decimal.frame.count[i, 2:11] <- NA
    }
  }
  
  ## Get maximum number of annual decimals
  max.no.decimals <- max(rowSums(decimal.frame.count[,2:11]), na.rm=TRUE)
  
  ## Prepare plot (max 30 years per segment)
  max.years <- 30
  segment.amounts.rounded <- ceiling(length(min(Data[,2]):max(Data[,2])) / max.years)
  start.year <- min(Data[,2])
  ## Loop on segments
  for (k in 1:segment.amounts.rounded) {
    seg <- decimal.frame.count[which(decimal.frame.count[,1] > start.year-1 & 
                                       decimal.frame.count[,1] < start.year+max.years),]
    ## Fill segment to length of max.years
    if (length(min(seg[,1]):max(seg[,1])) < max.years) {
      add.dates <- max.years - length(min(seg[,1]):max(seg[,1]))
      add.matrix <- matrix(nrow=add.dates, ncol=ncol(seg))
      add.matrix[, 1] <- (max(seg[,1])+1):(start.year+max.years-1)
      colnames(add.matrix) <- colnames(seg)
      seg <- rbind(seg, add.matrix)
    }
    # Name segment
    assign(paste0("segment.no.", k), seg)
    assign(paste0("pdf.width.segment.no.", k), 
           12*(length(start.year:min(start.year+max.years-1,max(Data[,2]))) / max.years))
    start.year <- start.year + max.years
  }
  
  ## Plot single or multiple time series segments in one pdf
  pdf(file = outfile, width = 12, height = 7*segment.amounts.rounded)
  plots <- list()
  rhg_cols <- c("black", "yellow", "orange", "red", "darkslateblue", "darkgray", 
                "magenta","blue", "cyan", "darkgreen")
  plot_title <- strsplit(outfile, "\\.pdf")[[1]]
  plot_title <- rev(strsplit(plot_title, "/")[[1]])[1]
  ## Loop on segments
  for (k in 1:segment.amounts.rounded) {
#    decimal.frame.count <- melt(get(paste0("segment.no.", k)), id.vars = "year")
    x <- get(paste0("segment.no.", k))
    decimal.frame.count <- data.frame(year = rep(x$year, 10),
                                      variable = rep(dec.labels, each=length(x$year)),
                                      value = as.vector(as.matrix(x[,2:11])))
    segment.plot <- (ggplot2::ggplot(decimal.frame.count, 
                                     ggplot2::aes(x=year, y=value, fill=variable))
                    + ggplot2::geom_bar(stat="identity", width=.7, 
                                        position=ggplot2::position_stack(reverse = TRUE)) 
                    + ggplot2::scale_fill_manual(values = rhg_cols)
                    + ggplot2::xlim(min(decimal.frame.count[,1])-0.5, 
                           max(decimal.frame.count[,1])+0.5)
                    + ggplot2::ylim(0, max.no.decimals)
                    + ggplot2::theme(plot.title = 
                                       ggplot2::element_text(size=22, face="bold", 
                                                             hjust=0.5), 
                            axis.text = ggplot2::element_text(size=18), 
                            axis.title = ggplot2::element_text(size=18, colour="black", 
                                                               face=NULL), 
                            axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r=15), 
                                                                 vjust=1.3), 
                            axis.title.x = ggplot2::element_blank(), 
                            legend.title = ggplot2::element_blank(), 
                            legend.text = ggplot2::element_text(colour="black", size=15), 
                            plot.margin = ggplot2::unit(c(0.6,0.2,0.5,0.5), "cm"))
                    + ggplot2::guides(fill = ggplot2::guide_legend(reverse=TRUE)) 
                    + ggplot2::ylab("observations / year")
                    + ggplot2::ggtitle(plot_title," \n")) 
    plots[[k]] <- segment.plot
  }
  
  multiplot(plotlist = plots)
  dev.off()
  
  message(paste("Plot saved in",outfile))
  
}


###############################################################################


#' Plot daily data points
#' 
#' Plot daily data points for custom intervals.
#'
#' @param dailydata A character string giving the path of the input file, or
#' a 5-column matrix with following columns: variable code, year, month, day,  
#' and the daily value.
#' @param len Integer indicating the number of years shown in 
#' each panel.
#' @param outfile Character string giving the path of the output pdf file.
#' @param startyear First year to plot. If not indicated, all available years
#' until \code{endyear} will be plotted.
#' @param endyear Last year to plot. If not indicated, all available years
#' since \code{startyear} will be plotted.
#' @param miss If TRUE (the default), missing data are plotted as red crosses at the bottom 
#' of the plot.
#' @param units Character string giving the units (will be printed in the y-axis). If
#' \code{dailydata} is a path to a file, then the units are read from the SEF header.
#' @param ...  Graphical parameters passed to the function \code{\link{plot}},
#' such as \code{cex}, \code{lwd}, \code{pch}, \code{col}, etc.
#' (see \code{\link{par}}).
#'
#' @details
#' The input file must follow the C3S Station Exchange Format (SEF).
#' 
#' Missing data are shown as red dots at the bottom of the plot.
#' 
#' @author Stefan Hunziker, Yuri Brugnara
#' 
#' @references
#' Hunziker et al., 2017: Identifying, attributing, and overcoming common data quality 
#' issues of manned station observations. Int. J. Climatol, 37: 4131-4145.
#' 
#' Hunziker et al., 2018: Effects of undetected data quality issues on climatological
#' analyses. Clim. Past, 14: 1-20.
#' 
#' @examples
#' plot_daily(Rosario$Tx, len = 2, outfile = paste0(tempdir(),"/test.pdf"))
#' 
#' @import graphics
#' @import grDevices
#' @import utils
#' @export


plot_daily <- function(dailydata, len = 1, outfile,                             
                       startyear = NA, endyear = NA, 
                       miss = TRUE, units = NA, ...) {
  
  if (is.null(dim(dailydata))) {
    units <- read_meta(dailydata, "units")
    dailydata <- read_sef(dailydata)[, 1:5]
  }
  y_lab <- paste(dailydata[1,1], ifelse(is.na(units), "", paste0("[",units,"]")))
  
  if (!is.na(startyear)) {
    dailydata <- dailydata[which(dailydata[,2] >= startyear), ]
  }
  
  if (!is.na(endyear)) {
    dailydata <- dailydata[which(dailydata[,2] <= endyear), ]
  }
  
  if (substr(outfile, nchar(outfile)-2, nchar(outfile)) != "pdf") {
    outfile <- paste(outfile, "pdf", sep = ".")
  }
  
  ## define grid
  abline_h <- as.numeric(pretty(round(min(dailydata[,5], na.rm = T)):
                           round(max(dailydata[,5], na.rm = T))))
  abline_v <- seq(ISOdate(min(dailydata[,2]), 1, 1, 0), 
                  ISOdate(max(dailydata[,2]) + 1, 1, 1, 0), "months")
  
  ## get segment lengths
  segment.amounts.rounded <- ceiling(length(min(dailydata[,2]):max(dailydata[,2])) / len)
  
  ## loop on the segments
  n <- 0
  start.year <- min(dailydata[,2])
  for (k in 1:segment.amounts.rounded) {
    seg <- dailydata[which(dailydata[,2] >= start.year & 
                             dailydata[,2] < start.year + 
                             len),]
    
    if (dim(seg)[1] > 0) {
      seg$Date <- ISOdate(seg[,2], seg[,3], seg[,4])
      n <- n + 1
      all.dates <- seq(ISOdate(min(seg[,2]), 1, 1), 
                       ISOdate(min(seg[,2]) + len - 1, 
                         12, 31), "days")
      seg.missing <- all.dates[which(!all.dates %in% seg$Date)]
      
      ## name segment
      assign(paste0("segment.no.", k), seg)
      assign(paste0("segment.missing.no.", k), seg.missing)
    }
    start.year <- start.year + len
  }
  
  ## plot single or multiple time series segments
  pdf(outfile, width = 7, height = n * 2)
  old_par <- list(mfrow = par()$mfrow, mar = par()$mar)
  on.exit(par(old_par))
  par(mar = c(3, 4, 2, 2))
  par(mfrow = c(n, 1))
  for (k in 1:segment.amounts.rounded) {
    if (exists(paste0("segment.no.", k))) {
      year1 <- min(get(paste0("segment.no.", k))[, 2])
      year2 <- max(get(paste0("segment.no.", k))[, 2])
      years <- ifelse(year1 != year2, paste0(year1, " - ", year2), year1)
      plot(get(paste0("segment.no.", k))$Date, 
           as.numeric(get(paste0("segment.no.", k))[,5]),
           panel.first = abline(h = abline_h, v = abline_v, col = "lightgray"),
           xlim = c(ISOdate(year1, 1, 1), 
                    ISOdate(year1 + len - 1, 12, 31)),
           ylim = c(min(dailydata[,5]-1, na.rm=T), max(dailydata[,5], na.rm=T)),
           xlab = "", ylab = y_lab, main = years, ...)
      if (length(get(paste0("segment.missing.no.", k))) > 0 & miss) {
        points(get(paste0("segment.missing.no.", k)), 
               rep(min(dailydata[,5], na.rm=T) - 1, 
                   length(get(paste0("segment.missing.no.", k)))), 
               pch = 4, col = "red", cex = 0.5)
      }
    }
  }
  
  dev.off()

  message(paste("Plot saved in",outfile))
  
}


###############################################################################


## Gives the number of days in a certain month

days_of_month <- function(year, month) {
  
  days <- 31
  if (month %in% c(4,6,9,11)) days <- 30
  if (month == 2) {
    days <- 28
    if (year %in% seq(1600,2200,4) & !(year %in% c(1700,1800,1900,2100,2200))) {
      days <- 29
    }
  }
  return(days)
  
}

###############################################################################


#' Plot sub-daily data points
#' 
#' Plot sub-daily data points divided by month.
#'
#' @param subdailydata A character string giving the path of the input file, or
#' a 7-column matrix with following columns: variable code, year, month, day,  
#' hour, minute, value.
#' @param year Integer vector giving the year(s) to plot. If not specified (NA), 
#' all available years will be plotted. One pdf per year will be created.
#' @param outfile Character string giving the path of the output pdf file. If  
#' \code{year} has more than one element or is NA, then this is a root to the 
#' filenames to which the year will be added ('root.year.pdf').
#' @param fixed If TRUE (default), use the same y axis for all months. If FALSE,
#' the axis limits are set based on the data range of each month.
#' @param units Character string giving the units (will be printed in the y-axis). 
#' If \code{subdailydata} is a path to a file, then the units are read from the 
#' SEF header.
#' @param time_offset Numeric vector of offsets in hours to be applied to the observation
#' times. Recycled if only one value is given. The default is no offset, 
#' i.e. UTC times for SEF input.
#' @param ...  Graphical parameters passed to the function \code{\link{plot}},
#' such as \code{cex}, \code{lwd}, \code{pch}, \code{col}, etc.
#' (see \code{\link{par}}).
#'
#' @details
#' Creates one pdf for each year plotted.
#' 
#' The input file must follow the C3S Station Exchange Format (SEF).
#' 
#' The parameter \code{time_offset} can be used to plot observations in local time
#' when reading the data in SEF.
#' 
#' @author Stefan Hunziker, Yuri Brugnara
#' 
#' @references
#' Hunziker et al., 2017: Identifying, attributing, and overcoming common data quality 
#' issues of manned station observations. Int. J. Climatol, 37: 4131-4145.
#' 
#' Hunziker et al., 2018: Effects of undetected data quality issues on climatological
#' analyses. Clim. Past, 14: 1-20.
#' 
#' @examples
#' plot_subdaily(Bern$p, year = 1803:1804, outfile = paste0(tempdir(),"/test"))
#' 
#' @import graphics
#' @import grDevices
#' @import utils
#' @export


plot_subdaily <- function(subdailydata, year = NA, outfile,
                          fixed = TRUE, units = NA, time_offset = 0, ...) {
  
  if (is.null(dim(subdailydata))) {
    units <- read_meta(subdailydata, "units")
    subdailydata <- read_sef(subdailydata)
  }
  v <- subdailydata[1, 1]
  y_lab <- paste(v, ifelse(is.na(units), "", paste0("[",units,"]")))
  
  if (is.na(year[1])) {
    year <- sort(unique(subdailydata[,2]))
  }
  subdailydata <- subdailydata[which(subdailydata[,2] %in% year), ]
  
  if (dim(subdailydata)[1] == 0) stop("No data to plot")
  
  ## define filename(s)
  if (length(year) > 1) {
    outfile <- paste(outfile, year, "pdf", sep = ".")
  } else if (substr(outfile, nchar(outfile)-2, nchar(outfile)) != "pdf") {
    outfile <- paste(outfile, "pdf", sep = ".")
  }
  names(outfile) <- year
  
  old_par <- list(mfrow = par()$mfrow, mar = par()$mar, xaxs = par()$xaxs)
  on.exit(par(old_par))
  
  subdailydata$Date <- ISOdate(subdailydata[,2], subdailydata[,3], subdailydata[,4], 
                               subdailydata[,5], subdailydata[,6]) 
  ## Apply offset
  if (!all(time_offset == 0)) {
    subdailydata$Date <- subdailydata$Date + time_offset * 3600
    subdailydata[,2] <- as.integer(format(subdailydata$Date,"%Y"))
    subdailydata[,3] <- as.integer(format(subdailydata$Date,"%m"))
  }
  
  for (y in year) {
    
    ## loop on the months
    for (k in 1:12) {
      seg <- subdailydata[which(subdailydata[,2] == y & subdailydata[,3] == k), ]
      
      if (dim(seg)[1] == 0) {
        seg <- as.data.frame(matrix(nrow = days_of_month(y,k), ncol = 7))
        names(seg) <- names(subdailydata)[1:7]
        seg$Date <- ISOdate(y, k, 1:days_of_month(y,k))
      }
        
      ## name segment
      assign(paste0("segment.no.", k), seg)
    }
    
    ## plot year
    pdf(outfile[as.character(y)], width = 7, height = 24)
    par(mar = c(3, 4, 2, 2))
    par(mfrow = c(12, 1))
  
    ## loop on the months
    for (k in 1:12) {
      
      ## define grid
      if (k == 12) {
        abline_v <- seq(ISOdate(y, k, 1, 0), ISOdate(y+1, 1, 1, 0), "days")
      } else {
        abline_v <- seq(ISOdate(y, k, 1, 0), ISOdate(y, k+1, 1, 0), "days")
      }
      if (fixed | sum(!is.na(get(paste0("segment.no.", k))[,7])) == 0) {
        abline_h <- as.numeric(pretty(min(subdailydata[,7], na.rm = T):
                                        max(subdailydata[,7], na.rm = T))) 
      } else {
        abline_h <- as.numeric(pretty(min(get(paste0("segment.no.", k))[,7], 
                                                na.rm = T):
                                        max(get(paste0("segment.no.", k))[,7], 
                                                  na.rm = T)))
      }
      
      plot(get(paste0("segment.no.", k))$Date, 
           as.numeric(get(paste0("segment.no.", k))[,7]),
           panel.first = abline(h = abline_h, v = abline_v, col = "lightgray"),
           xlim = c(abline_v[1]-3600*12, rev(abline_v)[1]+3600*12),
           ylim = c(min(abline_h), max(abline_h)),
           xlab = "", ylab = y_lab, xaxs = "i",
           main = paste(format(as.Date(paste(y,k,1,sep="-")),format="%B"), y), ...)
      
    }
    
    dev.off()
    par(old_par)
  
  }
  
  message(paste("Plot saved in",outfile))
  
}


###############################################################################


## Binomial test for the weekly cycle

weekly_test <- function(x, p) {
  
  wetdays <- as.integer(x$ramount >= 1)
  raindays.count <- aggregate(wetdays, list(x[,"weekday"]),
                              FUN = sum, na.rm = TRUE) # no. of wetdays per weekday
  alldays <- replace(wetdays, wetdays == 0, 1)
  weekdays.count <- aggregate(alldays, list(x[,"weekday"]),
                              FUN = sum, na.rm = TRUE) # total no. of days per weekday
  day.names <- weekdays(ISOdate(2019, 3, 11:17))
  
  mon_wd <- raindays.count[which(raindays.count[,1] == day.names[1]), 2]
  tue_wd <- raindays.count[which(raindays.count[,1] == day.names[2]), 2]
  wed_wd <- raindays.count[which(raindays.count[,1] == day.names[3]), 2]
  thu_wd <- raindays.count[which(raindays.count[,1] == day.names[4]), 2]
  fri_wd <- raindays.count[which(raindays.count[,1] == day.names[5]), 2]
  sat_wd <- raindays.count[which(raindays.count[,1] == day.names[6]), 2]
  sun_wd <- raindays.count[which(raindays.count[,1] == day.names[7]), 2]
  
  mon_all <- weekdays.count[which(raindays.count[,1] == day.names[1]), 2]
  tue_all <- weekdays.count[which(raindays.count[,1] == day.names[2]), 2]
  wed_all <- weekdays.count[which(raindays.count[,1] == day.names[3]), 2]
  thu_all <- weekdays.count[which(raindays.count[,1] == day.names[4]), 2]
  fri_all <- weekdays.count[which(raindays.count[,1] == day.names[5]), 2]
  sat_all <- weekdays.count[which(raindays.count[,1] == day.names[6]), 2]
  sun_all <- weekdays.count[which(raindays.count[,1] == day.names[7]), 2]
  
  wetall <- mon_wd + tue_wd + wed_wd + thu_wd + fri_wd + sat_wd + sun_wd
  daysall <- mon_all + tue_all + wed_all + thu_all + fri_all + sat_all + sun_all
  day_all <- c(mon_all, tue_all, wed_all, thu_all, fri_all, sat_all, sun_all)
  
  phi <- wetall / daysall # probability of a day to be a wd
  
  ## Binomial test
  ci_mon <- binom.test(mon_wd, mon_all, p = phi, alternative ="two.sided", 
                       conf.level = p)
  ciind_mon <- ifelse(!(ci_mon$conf.int[1] < phi & phi < ci_mon$conf.int[2]), 1, 0)
  ci_tue <- binom.test(tue_wd, tue_all, p = phi, alternative ="two.sided", 
                       conf.level = p)
  ciind_tue <- ifelse(!(ci_tue$conf.int[1] < phi & phi < ci_tue$conf.int[2]), 2, 0)
  ci_wed <- binom.test(wed_wd, wed_all, p = phi, alternative ="two.sided", 
                       conf.level = p) 
  ciind_wed <- ifelse(!(ci_wed$conf.int[1] < phi & phi < ci_wed$conf.int[2]), 3, 0)
  ci_thu <- binom.test(thu_wd, thu_all, p = phi, alternative ="two.sided", 
                       conf.level = p) 
  ciind_thu <- ifelse(!(ci_thu$conf.int[1] < phi & phi < ci_thu$conf.int[2]), 4, 0)
  ci_fri <- binom.test(fri_wd, fri_all, p = phi, alternative ="two.sided", 
                       conf.level = p) 
  ciind_fri <- ifelse(!(ci_fri$conf.int[1] < phi & phi < ci_fri$conf.int[2]), 5, 0)
  ci_sat <- binom.test(sat_wd, sat_all, p = phi, alternative ="two.sided", 
                       conf.level = p) 
  ciind_sat <- ifelse(!(ci_sat$conf.int[1] < phi & phi < ci_sat$conf.int[2]), 6, 0)
  ci_sun <- binom.test(sun_wd, sun_all, p = phi, alternative ="two.sided", 
                       conf.level = p) 
  ciind_sun <- ifelse(!(ci_sun$conf.int[1] < phi & phi < ci_sun$conf.int[2]), 7, 0)
  ci_out <- c(ciind_mon, ciind_tue, ciind_wed, ciind_thu, ciind_fri, ciind_sat, ciind_sun)

  mon_prob <- mon_wd / mon_all
  tue_prob <- tue_wd / tue_all
  wed_prob <- wed_wd / wed_all
  thu_prob <- thu_wd / thu_all
  fri_prob <- fri_wd / fri_all
  sat_prob <- sat_wd / sat_all
  sun_prob <- sun_wd / sun_all
  d.count <- c(mon_prob, tue_prob, wed_prob, thu_prob, fri_prob, sat_prob, sun_prob)
  output <- data.frame(1:7, d.count)
  
  return(list(phi=phi, day_all=day_all, ci_out=ci_out, output=output))
  
}

###############################################################################


#' Plot weekly cycle
#' 
#' Check if there is a significant weekly cycle in daily precipitation data by
#' means of a binomial test.
#'
#' @param dailypcp A character vector giving the paths of the input files, or
#' a list of 5-column matrices with following columns: variable code (must be 'rr'),   
#' year, month, day, value. The names of the list elements are assumed to be
#' the station IDs.
#' @param outpath Character string giving the path for the output files.
#' @param p Probability threshold for the binomial test (default is 0.95).
#'
#' @details
#' The input files must follow the C3S Station Exchange Format (SEF).
#' 
#' Creates one pdf for each station ('weekly.ID.pdf') plus one pdf with an overview 
#' of the entire dataset ('weekly.pdf').
#' 
#' @author Stefan Hunziker, Yuri Brugnara
#' 
#' @references
#' Hunziker et al., 2017: Identifying, attributing, and overcoming common data quality 
#' issues of manned station observations. Int. J. Climatol, 37: 4131-4145.
#' 
#' Hunziker et al., 2018: Effects of undetected data quality issues on climatological
#' analyses. Clim. Past, 14: 1-20.
#' 
#' @examples
#' plot_weekly_cycle(list(Rosario = Rosario$rr), outpath = tempdir())
#' 
#' @import graphics
#' @import grDevices
#' @import utils
#' @import stats
#' @export

plot_weekly_cycle <- function(dailypcp, outpath, p = 0.95) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  if (!requireNamespace("grid", quietly = TRUE)) {
    stop("Package \"grid\" needed for this function to work. Please install it.",
         call. = FALSE)
  }  
  
  if (nchar(outpath) > 0) {
    if (substr(outpath, nchar(outpath), nchar(outpath)) != "/") {
      outpath <- paste0(outpath, "/")
    }
  }
  
  ## Read data
  if (is.data.frame(dailypcp) | is.matrix(dailypcp)) {
    stop("dailypcp must be a list or a character vector") 
  }
  if (!is.list(dailypcp)) {
    tmp <- list()
    for (i in 1:length(dailypcp)) {
      station <- read_meta(dailypcp[i], "id")
      tmp[[station]] <- read_sef(dailypcp[i])[, 1:5]
    }
    dailypcp <- tmp
    rm(tmp)
  } 
  
  old_par <- list(mar = par()$mar)
  on.exit(par(old_par))
  
  ## Loop on stations
  for (st in names(dailypcp)) {
    
    if (dailypcp[[st]][1,1] != "rr") stop(paste0(st, ": not daily precipitation"))
    message(paste("Working on", st))
    
    tab.data <- data.frame(yr = dailypcp[[st]][, 2],
                           m = dailypcp[[st]][, 3],
                           d = dailypcp[[st]][, 4],
                           ramount = dailypcp[[st]][, 5])
    tab.data$weekday <- weekdays(as.Date(paste(tab.data$yr, tab.data$m, tab.data$d, sep="-")))
    
    bintest <- weekly_test(tab.data, p)
    output <- bintest$output

    ## Create and save plots
    pdf(file = paste0(outpath, "weekly.", st, ".pdf"), width = 7, height = 5)
    par(mar = c(2.5, 5, 3, 2))
    
    ## Create dummy-plot
    day.names <- weekdays(ISOdate(2019, 3, 11:17), abbreviate = TRUE)
    w.day <- 1:7
    plot(0, t = "p", pch = 16, col = "white",
         ylim = c(min(output$d.count) - 0.01, max(output$d.count) + 0.01),
         xlim = c(1, 7), axes = FALSE,
         xlab = "", ylab = "WD fraction (rr >= 1 mm)",
         main = st, cex.main = 1.8, cex.lab = 1.3)
    axis(1, w.day, day.names, cex.axis = 1.5)
    axis(2, cex.axis = 1.5,las = 0.05)
    abline(h = bintest$phi, lty = 2, col = "black")
    
    par(new = T)
    plot(output, type = "h", lwd = 8, pch = 4, 
         col = ifelse(w.day == bintest$ci_out, "red", "black"), 
         ylim = c(min(output$d.count) - 0.01, max(output$d.count) + 0.01),
         axes = F, xlab = "", ylab = "", cex = 1)
    text(output, labels = bintest$day_all, pos = 3)
    
    dev.off()
    
  }
  
  ## Prepare data for overview plot
  message("Working on overview plot...")
  ## Find start and end years
  n <- length(dailypcp)
  start <- array(dim = n)
  end <- start
  for (i in 1:n) {
    start[i] <- min(dailypcp[[names(dailypcp)[i]]][, 2])
    end[i] <- max(dailypcp[[names(dailypcp)[i]]][, 2])
  }
  
  ## Create dummy data frame to fill in if a year of a station is affected
  ci_table <- as.data.frame(matrix((c(1:(n+1) * (max(end)-min(start)+2)) * NA), 
                                   nrow = (max(end)-min(start)+2), 
                                   ncol = (n+1)))
  ci_table[2:(max(end)-min(start)+2), 1] <- min(start):max(end)
  
  ## Loop on stations
  for(i in 1:n) {
    
    station.name <- names(dailypcp)[i]
    tab.dailydata <- dailypcp[[station.name]]
    
    ## Create dummie vector for yearly - ci indication
    ci_vec <- array(dim = (max(tab.dailydata[,2]) - min(tab.dailydata[,2]) + 1)) 
    ci_vec_pos <- 1
    
    ## Loop over the years
    for (year in min(tab.dailydata[,2]):max(tab.dailydata[,2])) { 
      single_yr <- tab.dailydata[which(tab.dailydata[,2] == year), ]
      
      if (sum(!is.na(single_yr[,5])) > 180) { # exclude years that contain too many NAs
        tab.data <- data.frame(yr = single_yr[, 2],
                               m = single_yr[, 3],
                               d = single_yr[, 4],
                               ramount = single_yr[, 5])
        tab.data$weekday <- weekdays(as.Date(paste(tab.data$yr, tab.data$m, tab.data$d, sep="-")))
        
        bintest <- weekly_test(tab.data, p)
        output <- bintest$output
        ## ci vector: 0 = in ci, 1 = outside ci
        ci_vec[ci_vec_pos] <- ifelse(sum(bintest$ci_out, na.rm=TRUE) == 0, 0, 1) 
      }
      
      ci_vec_pos_new <- ci_vec_pos + 1 #position in ci-vec (changes with each year)
      ci_vec_pos <- ci_vec_pos_new
    }
    ci_table[1, i+1] <- station.name #write station name in table
    ci_table[(min(tab.dailydata[,2])-min(start)+2):(max(tab.dailydata[,2])-min(start)+2), i+1] <- 
      ci_vec #fill data in ci_table
  }

  ## Print data table in ggplot
  ## structure data to plot in ggplot (column station.name, column year, column ci_index)
  year_vec <- rep(ci_table[2:nrow(ci_table),1], ncol(ci_table)-1)
  station_vec <- array(dim = length(year_vec))
  pos_first <- 1
  pos_last <- nrow(ci_table) - 1
  for(name in 2:ncol(ci_table)) {
    station_vec[pos_first:pos_last] <- rep(ci_table[1,name], nrow(ci_table)-1)
    pos_first <- pos_first + nrow(ci_table) - 1
    pos_last <- pos_last + nrow(ci_table) - 1
  }
  data_vec <- array(dim = length(year_vec))
  pos_first <- 1
  pos_last <- nrow(ci_table) - 1
  for(name in 2:ncol(ci_table)) {
    data_vec[pos_first:pos_last] <- ci_table[2:nrow(ci_table),name]
    pos_first <- pos_first + nrow(ci_table) - 1
    pos_last <- pos_last + nrow(ci_table) - 1
  }
  ggplot_table <- data.frame(station_vec, year_vec, data_vec)
  name.vec <- as.vector(unique(ggplot_table[,1][!is.na(ggplot_table[,1])]))
  stations.page <- min(length(name.vec)+1, 20) # max 20 stations per page
  plot.no <- ceiling(length(name.vec)/stations.page) 
  
  #plot in ggplot
  pdf(file = paste0(outpath, "weekly.pdf"), width = max(9*(max(end)-min(start)+1)/100, 3), 
      height = 10.5*stations.page/20)
  for(station.no in 1:plot.no) {
    start.no <- station.no * stations.page - stations.page + 1
    end.no <- min(station.no * stations.page, length(name.vec))
    table <- ggplot_table[which(ggplot_table[,1]  %in% name.vec[start.no:end.no]), ]
    colnames(table) <- (c("station", "year", "ind"))
    plot <- ggplot2::ggplot(data=table, ggplot2::aes(x=year, y=station, colour=ind)) + 
      ggplot2::geom_point(size=1.5) + 
      ggplot2::scale_color_manual(breaks = c("0","1"),values=c("black","red")) + 
      ggplot2::facet_wrap(~ station, scales = 'free_y', ncol = 1) + 
      ggplot2::theme(legend.position="none", axis.title.y=ggplot2::element_blank(), 
                     axis.text.y=ggplot2::element_blank(), 
            axis.title.x=ggplot2::element_blank(), 
            panel.grid.major.y = ggplot2::element_blank(), 
            panel.grid.major = ggplot2::element_line(colour = "white",size=0.5), 
            axis.text=ggplot2::element_text(size=14), 
            plot.margin = ggplot2::unit(c(0.5,1,0.5,0.5), "cm"), 
            axis.ticks.length = ggplot2::unit(0.2,"cm")) + 
      ggplot2::xlim(min(ggplot_table[,2]), max(ggplot_table[,2]))
    print(plot + ggplot2::theme(axis.text.y = ggplot2::element_blank(), 
                                axis.ticks = ggplot2::element_blank()))
  }
  dev.off()
  
  message(paste("Plots saved in", outpath))
  
}