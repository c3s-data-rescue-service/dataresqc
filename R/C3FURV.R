#' Climatic outliers test
#' 
#' Considers as outliers all values falling outside a range between,
#' for example, p25 - 3 interquartile ranges and p75 + 3 interquartile 
#' The number of interquantile ranges can be modified through the parameter
#' \code{IQR}. 
#' 
#' @param Data A character string giving the path of the input file,
#' or a matrix with 5 (7) columns for daily (sub-daily) data: variable code, year, 
#' month, day, (hour), (minute), value.
#' @param meta A character vector with 6 elements: station ID, latitude, longitude,
#' altitude, variable code, units. If \code{Data} is a path, \code{meta} is
#' ignored.
#' @param outpath Character string giving the path for the QC results. By default
#' this is the working directory.
#' @param IQR Interquantile range used to define outliers. By default it is 5 for
#' precipitation, 3 for air temperature, and 4 for any other variable.
#' @param bplot If TRUE, create a boxplot and print it into a PDF.
#' @param outfile Filename for the plot. Ignored if \code{bplot} is FALSE.
#' @param ...  Graphical parameters passed to the function \code{\link{boxplot}}.
#' 
#' @details The input file must follow the Copernicus Station Exchange 
#' Format (SEF). This function works with any numerical variable.
#' 
#' Zeroes are automatically excluded in bounded variables such as precipitation.
#' 
#' @author Alba Gilabert, Yuri Brugnara
#'
#' @examples 
#' climatic_outliers(Rosario$Tn, Meta$Tn, IQR = 4) 
#'
#' @import graphics
#' @import grDevices
#'
#' @export

climatic_outliers <- function(Data, meta = NULL, outpath = getwd(), 
                              IQR = NA, bplot = FALSE, outfile = NA, ...) {
  
  #Read data and metadata
  if (is.null(dim(Data))) {
    meta <- read_meta(Data, c("id","lat","lon","alt","var","units"))
    Data <- read_sef(Data)
  }
  if (length(meta) != 6) stop("Incorrect metadata (must be a vector of 6 elements)")
  n <- dim(Data)[2]
  if (!n %in% c(5,7)) stop("Incorrect dimension of Data")
  meta <- as.character(meta)
  
  #Define IQR and flag
  if (is.na(IQR)) {
    if (meta[5] == "rr") {
      outrange <- 5
    } else if (meta[5] %in% c("Tx","Tn","ta")) {
      outrange <- 3
    } else {
      outrange <- 4
    }
  } else {
    outrange <- IQR
  }
  flag <- "climatic_outliers"
  
  #Enough data to compute outliers
  if (sum(!is.na(Data[,n])) > 5*365) {
    
    #Remove zeroes
    if (meta[5] %in% c("rr","sd","fs","sc","sw")) {
      Data <- Data[which(Data[,n] != 0), ]
    }
    
    #Make boxplot
    if (bplot) {
      if (is.na(outfile)) {
        namefile <- paste0('climatic_outliers_boxplot_', 
                           meta[1], "_", meta[5], '.pdf')
      } else {
        namefile <- outfile
      }
      pdf(file = namefile)
    }
    x <- boxplot(Data[,n] ~ Data[,3], main = meta[5],
                 range = outrange, plot = bplot,
                 xlab = "Months", ylab = meta[6], ...)
    if (bplot) dev.off()
    
    #Flag data
    out <- data.frame(Var=character(), Year=numeric(), Month=numeric(),
                      Day=numeric(), Value=numeric(), Meta=character())
    for (a in 1:12) {
      out <- rbind(out, subset(Data, Data[,3] == a &
                                 ((Data[,n] > x$stats[5,a]) | (Data[,n] < x$stats[1,a]))))
    }
    
    #Write flags to file
    if (nrow(out) != 0) {
      out$Test <- flag
      res <- ifelse(n == 5, "daily", "subdaily")
      namefile <- paste0("qc_", meta[1], "_", meta[5], "_", res, ".txt")
      namefile <- paste(file.path(outpath), namefile, sep="/")
      if (n == 5) {
        write_intermediate_daily(namefile, out)
      } else {
        write_intermediate_subdaily(namefile, out)
      }
    }
    outMsg <- "Climatic outliers test completed"
    
  } else {
    outMsg <- paste("Not enough data for outliers test", 
                    "(minimum 5 years of non-zero values for daily data)")
  }
  
  return(print(outMsg, quote=FALSE))
  
}


################################################################################


#' Daily internal consistency test
#' 
#' Determines the coherence between daily maximum temperature (Tx) values and
#' daily minimum temperature (Tn) values; daily wind speed (w) and wind
#' direction (dd); daily snow cover (sc) and snow depth (sd); daily fresh snow 
#' (fs) and snow depth (sd); daily fresh snow (fs) and minimum temperature (Tn);
#' daily snow depth (sd) and minimum temperature (Tn). 
#'
#' @param dailydata A character vector giving the paths of two input files,
#' or a 5-column matrix with following columns: variable code, year, month,
#' day, and the daily value.
#' @param meta A data frame with 2 rows and 6 columns: station ID, latitude, 
#' longitude, altitude, variable code, units. If \code{dailydata} is a vector, 
#' \code{meta} is ignored.
#' @param outpath Character string giving the path for the QC results. By default
#' this is the working directory.
#' 
#' @details The input file must follow the Copernicus Station Exchange 
#' Format (SEF).
#' 
#' The daily minimum temperature is assumed to be observed at the same time of 
#' the snow depth / fresh snow, and to refer to the same 24-hour period. Snow 
#' accumulation is flagged if the minimum temperature is
#' higher than 3 degrees Celsius.
#'
#' @author Alba Gilabert, Yuri Brugnara
#' 
#' @examples 
#' internal_consistency(rbind(Rosario$Tx, Rosario$Tn), 
#'                      rbind(Meta$Tx, Meta$Tn))
#'
#' @export

internal_consistency <- function(dailydata, meta = NULL, outpath = getwd()) {
  
  #Read data and metadata
  if (is.null(dim(dailydata))) {
    if (length(dailydata) != 2) {
      stop("Incorrect dimension of files (2 input files required)")
    }
    meta <- rbind(read_meta(dailydata[1], c("id","lat","lon","alt","var","units")), 
                  read_meta(dailydata[2], c("id","lat","lon","alt","var","units")))
    dailydata <- rbind(read_sef(dailydata[1]), read_sef(dailydata[2]))
  }
  if (dim(meta)[1] != 2 | dim(meta)[2] != 6) {
    stop("Incorrect metadata (must be a data frame with 6 columns and 2 rows)")
  }
  if (dim(dailydata)[2] != 5) stop("Incorrect dimension of dailydata")
  names(dailydata) <- c('Var', 'Year', 'Month', 'Day', 'Value')
  
  #Check units
  dailydata$Value[which(dailydata$Var==meta[1,5])] <- 
    check_units(dailydata$Value[which(dailydata$Var==meta[1,5])], meta[1,5], meta[1,6])
  dailydata$Value[which(dailydata$Var==meta[2,5])] <- 
    check_units(dailydata$Value[which(dailydata$Var==meta[2,5])], meta[2,5], meta[2,6]) 
  
  #Define flag
  flag <- "internal_consistency"
  
  #Check that a compatible pair of variables is given
  if (all(sort(meta[,5]) == c("Tn","Tx")) |
      all(sort(meta[,5]) == c("dd","w")) |
      all(sort(meta[,5]) == c("sc","sd")) |
      all(sort(meta[,5]) == c("fs","sd")) |
      all(sort(meta[,5]) == c("fs","Tn")) |
      all(sort(meta[,5]) == c("sd","Tn"))) {
    
    #Internal consistency tmax-tmin
    if ("Tx" %in% meta[,5]) {
      datax <- dailydata[which(dailydata[,1]=="Tx"), ]
      datan <- dailydata[which(dailydata[,1]=="Tn"), ]
      datatxtn <- merge(datax, datan, by = c('Year', 'Month', 'Day'))
      consistency <- subset(datatxtn, datatxtn$Value.x < datatxtn$Value.y)
    }
    
    #Internal consistency wind speed and wind direction
    if ("w" %in% meta[,5]) {
      dataw <- dailydata[which(dailydata[,1]=="w"), ]
      datadd <- dailydata[which(dailydata[,1]=="dd"), ]
      datawdd <- merge(dataw, datadd, by = c('Year', 'Month', 'Day'))
      consistency <- subset(datawdd, datawdd$Value.x == 0 & !is.na(datawdd$Value.y))
    }
    
    #Internal consistency snow cover and snow depth
    if ("sc" %in% meta[,5]) {
      datasc <- dailydata[which(dailydata[,1]=="sc"), ]
      datasd <- dailydata[which(dailydata[,1]=="sd"), ]
      datascsd <- merge(datasc, datasd, by = c('Year', 'Month', 'Day'))
      consistency <- subset(datascsd, datascsd$Value.x == 0 & datascsd$Value.y > 0)
    }
    
    #Internal consistency fresh snow and snow depth / minimum temperature
    if ("fs" %in% meta[,5]) {
      datafs <- dailydata[which(dailydata[,1]=="fs"), ]
      if ("sd" %in% meta[,5]) {
        datasd <- dailydata[which(dailydata[,1]=="sd"), ]
        datafssd <- merge(datafs, datasd, by = c('Year', 'Month', 'Day'))
        n <- dim(datafssd)[1]
        diffsd <- c(NA, datafssd$Value.y[2:n] - datafssd$Value.y[1:(n-1)])
        dates <- ISOdate(datafssd$Year, datafssd$Month, datafssd$Day)
        diffdate <- c(NA, as.integer(dates[2:n] - dates[1:(n-1)]))
        consistency <- subset(datafssd, diffsd > 0 & datafssd$Value.x == 0 & diffdate == 1)
      } else {
        datatn <- dailydata[which(dailydata[,1]=="Tn"), ]
        datafstn <- merge(datafs, datatn, by = c('Year', 'Month', 'Day'))
        consistency <- subset(datafstn, datafstn$Value.x > 0 & datafstn$Value.y > 3)
      }
    }
    
    #Internal consistency snow depth and minimum temperature
    if ("sd" %in% meta[,5] & "Tn" %in% meta[,5]) {
      datasd <- dailydata[which(dailydata[,1]=="sd"), ]
      datatn <- dailydata[which(dailydata[,1]=="Tn"), ]
      datasdtn <- merge(datasd, datatn, by = c('Year', 'Month', 'Day'))
      n <- dim(datasdtn)[1]
      diffsd <- c(NA, datasdtn$Value.y[2:n] - datasdtn$Value.y[1:(n-1)])
      dates <- ISOdate(datasdtn$Year, datasdtn$Month, datasdtn$Day)
      diffdate <- c(NA, as.integer(dates[2:n] - dates[1:(n-1)]))
      consistency <- subset(datasdtn, diffsd > 0 & datafssd$Value.y > 2.5 & diffdate == 1)
    }
    
    #Write flags to file
    if (nrow(consistency) != 0) {
      consistency$Test <- flag
      namefile <- paste0("qc_", meta[which(meta[,5]==consistency$Var.x[1]),1], "_", 
                         meta[which(meta[,5]==consistency$Var.x[1]),5], "_daily.txt")
      namefile <- paste(file.path(outpath), namefile, sep="/")
      out <- consistency[,c("Var.x","Year","Month","Day","Value.x","Test")]
      names(out) <- c("Var", "Year", "Month", "Day", "Value", "Test")
      write_intermediate_daily(namefile, out)
      namefile <- paste0("qc_", meta[which(meta[,5]==consistency$Var.y[1]),1], "_", 
                         meta[which(meta[,5]==consistency$Var.y[1]),5], "_daily.txt")
      namefile <- paste(file.path(outpath), namefile, sep="/")
      out <- consistency[,c("Var.y","Year","Month","Day","Value.y","Test")]
      names(out) <- c("Var", "Year", "Month", "Day", "Value", "Test")
      write_intermediate_daily(namefile, out)
    }
    
    outMsg <- "Internal consistency test completed"
    
  } else {
    outMsg <- "The variables provided are incompatible with this test"
  }
  
  return(print(outMsg, quote=FALSE))
}


################################################################################


#' Daily temporal coherence test
#' 
#' Find those records where daily maximum or minimum temperature, 
#' mean wind speed, 
#' snow depth, snow cover, or fresh snow differences with previous day are too
#' large.
#'
#' @param dailydata A character string giving the path of the input file,
#' or a 5-column matrix with following columns: variable code, year, month,
#' day, and the daily value.
#' @param meta A character vector with 6 elements: station ID, latitude, longitude,
#' altitude, variable code, units. If \code{dailydata} is a path, \code{meta} is
#' ignored.
#' @param outpath Character string giving the path for the QC results. By default
#' this is the working directory.
#' @param temp_jumps given a daily maximum or minimum temperature values of two
#' consecutive days, maximum difference in degrees Celsius. By default, 
#' temp_jumps = 20 C.
#' @param windspeed_jumps given a daily mean wind speed value of two consecutive days,
#' maximum difference in metres per second. By default, wind_jumps = 15 m/s.
#' @param snowdepth_jumps given a daily snow depth of two consecutive days, 
#' maximum difference in centimetres. By default, snowdepth_jumps = 50 cm.
#' 
#' @details The input file must follow the Copernicus Station Exchange 
#' Format (SEF).
#' 
#' @author Alba Gilabert, Yuri Brugnara
#'
#' @examples
#' temporal_coherence(Rosario$Tx, Meta$Tx, temp_jumps = 10)
#'
#' @export

temporal_coherence <- function(dailydata, meta = NULL, outpath = getwd(), 
                               temp_jumps = 20, windspeed_jumps = 15, 
                               snowdepth_jumps = 50) {
  
  #Read data and metadata
  if (is.null(dim(dailydata))) {
    meta <- read_meta(dailydata, c("id","lat","lon","alt","var","units"))
    dailydata <- read_sef(dailydata)
  }
  if (length(meta) != 6) stop("Incorrect metadata (must be a vector of 6 elements)")
  if (dim(dailydata)[2] != 5) stop("Incorrect dimension of dailydata")
  meta <- as.character(meta)

  #Check variable
  if (meta[5] %in% c("Tx","Tn","w","sd")) {
    
    #Check units
    dailydata[,5] <- check_units(dailydata[,5], meta[5], meta[6])
    
    #Define jumps and flag
    if (meta[5] %in% c("Tx","Tn")) {
      jumps <- temp_jumps
    }
    if (meta[5] == "w") {
      jumps <- windspeed_jumps
    }
    if (meta[5] == "sd") {
      jumps <- snowdepth_jumps
    }
    flag <- "temporal_coherence"
    
    #Flag jumps
    dailydata$Date <- as.Date(paste(dailydata[,2], dailydata[,3], 
                                    dailydata[,4], sep = "-"))
    n <- dim(dailydata)[1]
    dailydiff <- c(NA, dailydata[2:n,5] - dailydata[1:(n-1),5])
    datediff <- c(NA, as.integer(difftime(dailydata$Date[2:n],
                                          dailydata$Date[1:(n-1)],
                                          units="days")))
    flags <- which(abs(dailydiff) > jumps & datediff == 1)
    flags <- append(flags, flags+1)
    out <- dailydata[flags, 1:5]
    if (length(flags) != 0) {
      out <- out[!duplicated(out), ]
      out$Test <- flag
      namefile <- paste0("qc_", meta[1], "_", meta[5], "_daily.txt")
      namefile <- paste(file.path(outpath), namefile, sep="/")
      write_intermediate_daily(namefile, out)
    }
    
    outMsg <- "Temporal coherence test completed"
  } else {
    outMsg <- "Variable not supported by this test"
  }
  
  return(print(outMsg, quote=FALSE))
}


################################################################################


#' Daily repetition test
#' 
#' Report occurrences of equal consecutive values in daily data.
#' 
#' @param dailydata A character string giving the path of the input file,
#' or a 5-column matrix with following columns: variable code, year, month,
#' day, and the daily value.
#' @param meta A character vector with 6 elements: station ID, latitude, longitude,
#' altitude, variable code, units. If \code{dailydata} is a path, \code{meta} is
#' ignored.
#' @param outpath Character string giving the path for the QC results. By default
#' this is the working directory.
#' @param n Number of minimum equal consecutive values required for a flag. The
#' default is 4.
#' 
#' @details 
#' The input file must follow the Copernicus Station Exchange Format (SEF).
#' 
#' Zeroes are automatically excluded in bounded variables such as precipitation.
#' 
#' @author Alba Gilabert, Yuri Brugnara
#'
#' @examples
#' daily_repetition(Rosario$Tx, Meta$Tx, n = 3)
#'
#' @export

daily_repetition <- function(dailydata, meta = NULL, outpath = getwd(), n = 4) {
  
  #Read data and metadata
  if (is.null(dim(dailydata))) {
    meta <- read_meta(dailydata, c("id","lat","lon","alt","var","units"))
    dailydata <- read_sef(dailydata)
  }
  if (length(meta) != 6) stop("Incorrect metadata (must be a vector of 6 elements)")
  if (dim(dailydata)[2] != 5) stop("Incorrect dimension of dailydata")
  meta <- as.character(meta)
  
  #Define flag
  flag <- "daily_repetition"
  
  #Flag data
  vector_count <- rle(dailydata[,5])$lengths
  id_value <- cumsum(vector_count)
  value_flat <- which(vector_count > n-1)
  flat <- c()
  for (i in value_flat) flat <- append(flat, (id_value[i]-vector_count[i]+1):id_value[i])
  if (meta[5] %in% c("rr","sd","fs","sc","sw")) {
    flat <- flat[which(dailydata[flat,5] != 0)]
  }
  out <- dailydata[flat, ]
  if (nrow(out) != 0) {
    out$Test <- flag
    namefile <- paste0("qc_", meta[1], "_", meta[5], "_daily.txt")
    namefile <- paste(file.path(outpath), namefile, sep="/")
    write_intermediate_daily(namefile, out)
  }
  
  return(print("Repetition test completed", quote=FALSE))
}


################################################################################


#' Duplicate dates test
#' 
#' Flag dates that appear more than once in daily data.
#'
#' @param dailydata A character string giving the path of the input file,
#' or a 5-column matrix with following columns: variable code, year, month,
#' day, and the daily value.
#' @param meta A character vector with 6 elements: station ID, latitude, longitude,
#' altitude, variable code, units. If \code{dailydata} is a path, \code{meta} is
#' ignored.
#' @param outpath Character string giving the path for the QC results. By default
#' this is the working directory.
#' 
#' @details 
#' The input file must follow the Copernicus Station Exchange Format (SEF).
#' 
#' @author Alba Gilabert, Yuri Brugnara
#'
#' @examples 
#' duplicate_dates(Rosario$Tx, Meta$Tx)
#'
#' @export

duplicate_dates <- function(dailydata, meta = NULL, outpath = getwd()) {
  
  #Read data and metadata
  if (is.null(dim(dailydata))) {
    meta <- read_meta(dailydata, c("id","lat","lon","alt","var","units"))
    dailydata <- read_sef(dailydata)
  }
  if (length(meta) != 6) stop("Incorrect metadata (must be a vector of 6 elements)")
  n <- dim(dailydata)[2]
  if (n != 5) stop("Incorrect dimension of dailydata")
  meta <- as.character(meta)
  
  #Define flag
  flag <- "duplicate_dates"
  
  #Find duplicates
  dupli <- duplicated(dailydata[,2:4])
  dupli <- unique(c(which(dupli) - 1, which(dupli)))
  out <- dailydata[dupli, ]
  if (nrow(out) != 0) {
    out$Test <- flag
    namefile <- paste0("qc_", meta[1], "_", meta[5], "_daily.txt")
    namefile <- paste(file.path(outpath), namefile, sep="/")
    write_intermediate_daily(namefile, out)
  }
  
  return(print("Duplicate dates test completed", quote=FALSE))
}


################################################################################


#' Daily big errors test
#' 
#' Find the daily maximum, minimum, precipitation, mean wind direction,
#' mean wind speed, snow cover and snow depth that exceed thresholds selected by the
#' user. The output is a list with the days in which Tx, Tn, rr, dd, w, sc, sd
#' or fs exceeds some threshold.
#' 
#' @param dailydata A character string giving the path of the input file,
#' or a 5-column matrix with following columns: variable code, year, month,
#' day, and the daily value.
#' @param meta A character vector with 6 elements: station ID, latitude, longitude,
#' altitude, variable code, units. If \code{dailydata} is a path, \code{meta} is
#' ignored.
#' @param outpath Character string giving the path for the QC results. By default
#' this is the working directory.
#' @param tmax_upper is the tx maximum threshold in degrees Celsius.
#' By default, tmax_upper = 45 C.
#' @param tmax_lower is the tx minimum threshold in degrees Celsius.
#' By default, tmax_lower = -30 C.
#' @param tmin_upper is the tn maximum threshold in degrees Celsius.
#' By default, tmin_upper = 30 C.
#' @param tmin_lower is the tn minimum threshold in degrees Celsius.
#' By default, tmin_lower = -40 C.
#' @param rr_upper is the rr maximum threshold in millimetres.
#' By default, rr_upper = 200 mm.
#' @param rr_lower is the rr minimum threshold in millimetres.
#' By default, rr_lower = 0 mm. 
#' @param w_upper is the w maximum threshold in metres per second.
#' By default, w_upper = 30 m/s.
#' @param w_lower is the w mimumum threshold in metres per second.
#' By default, w_lower = 0 m/s.
#' @param dd_upper is the dd maximum threshold in degrees North.
#' By default, dd_upper = 360.
#' @param dd_lower is the dd minimum threshold in degrees North.
#' By default, dd_lower = 0.
#' @param sc_upper is the sc maximum threshold in percent.
#' By default, sc_upper = 100\%.
#' @param sc_lower is the sc minimum threshold in percent.
#' By default, sc_lower = 0\%.
#' @param sd_upper is the sd maximum threshold in centimetres.
#' By default, sd_upper = 200 cm.
#' @param sd_lower is the sd minimum threshold in centimetres.
#' By default, sd_lower = 0 cm.
#' @param fs_upper is the fs maximum threshold in centimetres.
#' By default, fs_upper = 100 cm.
#' @param fs_lower is the fs minimum threshold in centimetres.
#' By default, fs_lower = 0 cm.
#' 
#' @details 
#' The input file must follow the Copernicus Station Exchange Format (SEF).
#'
#' @author Alba Gilabert, Yuri Brugnara
#'
#' @examples 
#' daily_out_of_range(Rosario$Tn, Meta$Tn, tmin_upper = 25)
#'
#' @export

daily_out_of_range <- function(dailydata, meta = NULL, outpath = getwd(), tmax_upper = 45, 
                               tmax_lower = -30, tmin_upper = 30, tmin_lower = -40, 
                               rr_upper = 200, rr_lower = 0, w_upper = 30, 
                               w_lower = 0, dd_upper = 360, dd_lower =  0, 
                               sc_upper = 100, sc_lower = 0, sd_upper = 200, 
                               sd_lower = 0, fs_upper = 100, fs_lower = 0) {
  
  #Read data and metadata
  if (is.null(dim(dailydata))) {
    meta <- read_meta(dailydata, c("id","lat","lon","alt","var","units"))
    dailydata <- read_sef(dailydata)
  }
  if (length(meta) != 6) stop("Incorrect metadata (must be a vector of 6 elements)")
  n <- dim(dailydata)[2]
  if (n != 5) stop("Incorrect dimension of dailydata")
  names(dailydata) <- c("Var", "Year", "Month", "Day", "Value")
  meta <- as.character(meta)
  
  #Check variable
  if (meta[5] %in% c("Tx","Tn","rr","w","dd","sc","sd","fs")) {
    
    #Check units
    dailydata[,5] <- check_units(dailydata[,5], meta[5], meta[6])
    
    #Define thresholds and flag  
    if (meta[5] == "Tx") {
      upper <- tmax_upper
      lower <- tmax_lower
    }
    if (meta[5] == "Tn") {
      upper <- tmin_upper
      lower <- tmin_lower
    }
    if (meta[5] == "rr") {
      upper <- rr_upper
      lower <- rr_lower
    }
    if (meta[5] == "w") {
      upper <- w_upper
      lower <- w_lower
    }
    if (meta[5] == "dd") {
      upper <- dd_upper
      lower <- dd_lower
    }
    if (meta[5] == "sc") {
      upper <- sc_upper
      lower <- sc_lower
    }
    if (meta[5] == "sd") {
      upper <- sd_upper
      lower <- sd_lower
    }
    if (meta[5] == "fs") {
      upper <- fs_upper
      lower <- fs_lower
    }
    flag <- "daily_out_of_range"
    
    #Flag data
    out <- subset(dailydata, dailydata$Value > upper | dailydata$Value < lower)
    if (nrow(out) != 0) {
      out$Test <- flag
      namefile <- paste0("qc_", meta[1], "_", meta[5], "_daily.txt")
      namefile <- paste(file.path(outpath), namefile, sep="/")
      write_intermediate_daily(namefile, out)
    }
    outMsg <- "Big errors test completed"
    
  } else {
    outMsg <- "Variable not supported by this test"
  }
  
  return(print(outMsg, quote=FALSE))
}


################################################################################


#' Sub-daily big errors test
#' 
#' Find the subdaily temperature (ta), wind speed (w), 
#' wind direction (dd), snow cover (sc), snow depth (sd)
#' and fresh snow (fs) values that exceed thresholds selected by the
#' user. The output is a list with the days in which ta, rr, dd, w, sc, sd
#' or fs exceeds some threshold.
#' 
#' @param subdailydata A character string giving the path of the input file,
#' or a 7-column matrix with following columns: variable code, year, month,
#' day, hour, minute, value.
#' @param meta A character vector with 6 elements: station ID, latitude, longitude,
#' altitude, variable code, units. If \code{subdailydata} is a path, \code{meta} is
#' ignored.
#' @param outpath Character string giving the path for the QC results. By default
#' this is the working directory.
#' @param time_offset Offset in hours to add to the time to obtain local time.
#' By default, time_offset = 0.
#' @param ta_day_upper is the ta maximum day threshold in degrees Celsius.
#' By default, ta_day_upper = 45 C.
#' @param ta_day_lower is the ta minimum day threshold in degrees Celsius.
#' By default, ta_day_lower = -35 C.
#' @param ta_night_upper is the ta maximum night threshold in degrees Celsius.
#' By default, ta_night_upper = 40 C.
#' @param ta_night_lower is the ta minimum night threshold in degrees Celsius.
#' By default, ta_night_lower = -40 C.
#' @param rr_upper is the rr maximum threshold in millimetres.
#' By default, rr_upper = 100 mm.
#' @param rr_lower is the rr minimum threshold in millimetres.
#' By default, rr_lower = 0 mm. 
#' @param w_upper is the w maximum threshold in metres per second.
#' By default, w_upper = 50 m/s.
#' @param w_lower is the w mimumum threshold in metres per second.
#' By default, w_lower = 0 m/s.
#' @param dd_upper is the dd maximum threshold in degrees North.
#' By default, dd_upper = 360.
#' @param dd_lower is the dd minimum threshold in degrees North.
#' By default, dd_lower = 0.
#' @param sc_upper is the sc maximum threshold in percent.
#' By default, sc_upper = 100\%.
#' @param sc_lower is the sc minimum threshold in percent.
#' By default, sc_lower = 0\%.
#' @param sd_upper is the sd maximum threshold in centimetres.
#' By default, sd_upper = 200 cm.
#' @param sd_lower is the sd minimum threshold in centimetres.
#' By default, sd_lower = 0 cm.
#' @param fs_upper is the fs maximum threshold in centimetres.
#' By default, fs_upper = 100 cm.
#' @param fs_lower is the fs minimum threshold in centimetres.
#' By default, fs_lower = 0 cm.
#' 
#' @details 
#' The input file must follow the Copernicus Station Exchange Format (SEF).
#'
#' @author Alba Gilabert, Yuri Brugnara
#'
#' @examples 
#' subdaily_out_of_range(Rosario$ta, Meta$ta[which(Meta$ta$id=="Rosario"),], 
#'                       time_offset = -4.28, ta_day_upper = 35)
#'
#' @export

subdaily_out_of_range <- function(subdailydata, meta = NULL, outpath = getwd(), time_offset = 0,
                                  ta_day_upper = 45, ta_day_lower = -35,
                                  ta_night_upper = 40, ta_night_lower = -40,
                                  rr_upper = 100, rr_lower = 0, w_upper = 50, 
                                  w_lower = 0, dd_upper = 360, dd_lower =  0, 
                                  sc_upper = 100, sc_lower = 0, sd_upper = 200, 
                                  sd_lower = 0, fs_upper = 100, fs_lower = 0) {
  
  #Read data and metadata
  if (is.null(dim(subdailydata))) {
    meta <- read_meta(subdailydata, c("id","lat","lon","alt","var","units"))
    subdailydata <- read_sef(subdailydata)
  }
  if (length(meta) != 6) stop("Incorrect metadata (must be a vector of 6 elements)")
  n <- dim(subdailydata)[2]
  if (n != 7) stop("Incorrect dimension of subdailydata")
  names(subdailydata) <- c("Var", "Year", "Month", "Day", "Hour", "Minute", "Value")
  meta <- as.character(meta)
  
  #Check variable
  if (meta[5] %in% c("ta","rr","w","dd","sc","sd","fs")) {
    
    #Check units
    subdailydata[,7] <- check_units(subdailydata[,7], meta[5], meta[6])
    
    #Define thresholds and flag  
    if (meta[5] == "rr") {
      upper <- rr_upper
      lower <- rr_lower
    }
    if (meta[5] == "w") {
      upper <- w_upper
      lower <- w_lower
    }
    if (meta[5] == "dd") {
      upper <- dd_upper
      lower <- dd_lower
    }
    if (meta[5] == "sc") {
      upper <- sc_upper
      lower <- sc_lower
    }
    if (meta[5] == "sd") {
      upper <- sd_upper
      lower <- sd_lower
    }
    if (meta[5] == "fs") {
      upper <- fs_upper
      lower <- fs_lower
    }
    flag <- "subdaily_out_of_range"
    
    #Flag data
    if (meta[5] == "ta") {
      
      #Add time offset
      times <- ISOdate(subdailydata[,2], subdailydata[,3], subdailydata[,4],
                           subdailydata[,5], subdailydata[,6])
      times <- times + time_offset * 3600
      hours <- as.integer(substr(times,12,13))
      
      i_flags <- which((hours >= 8 & hours <= 19 & 
                          (subdailydata$Value < ta_day_lower | 
                             subdailydata$Value > ta_day_upper)) | 
                         ((hours < 8 | hours > 19) & 
                            (subdailydata$Value < ta_night_lower | 
                               subdailydata$Value > ta_night_upper)))
      out <- data.frame(Var = rep(as.character(meta[5]),length(i_flags)),
                        Year = subdailydata$Year[i_flags],
                        Month = subdailydata$Month[i_flags],
                        Day = subdailydata$Day[i_flags],
                        Hour = subdailydata$Hour[i_flags],
                        Minute = subdailydata$Minute[i_flags],
                        Value = subdailydata$Value[i_flags],
                        stringsAsFactors = FALSE)
      
    } else {
      out <- subset(subdailydata, subdailydata$Value > upper | subdailydata$Value < lower)
    }
    
    if (nrow(out) != 0) {
      out$Test <- flag
      namefile <- paste0("qc_", meta[1], "_", meta[5], "_subdaily.txt")
      namefile <- paste(file.path(outpath), namefile, sep="/")
      write_intermediate_subdaily(namefile, out)
    }
    outMsg <- "Big errors test completed"
    
  } else {
    outMsg <- "Variable not supported by this test"
  }
  
  return(print(outMsg, quote=FALSE))
}


################################################################################


#' Sub-daily repetition test
#' 
#' Report occurrences of equal consecutive values in subdaily data.
#' 
#' @param subdailydata A character string giving the path of the input file,
#' or a 7-column matrix with following columns: variable code, year, month,
#' day, hour, minute, value.
#' @param meta A character vector with 6 elements: station ID, latitude, longitude,
#' altitude, variable code, units. If \code{subdailydata} is a path, \code{meta} is
#' ignored.
#' @param outpath Character string giving the path for the QC results. By default
#' this is the working directory.
#' @param n Number of minimum equal consecutive values required for a flag. The
#' default is 6.
#' 
#' @details
#' The input file must follow the Copernicus Station Exchange Format (SEF).
#' 
#' Zeroes are automatically excluded in bounded variables such as precipitation.
#' 
#' @author Alba Gilabert, Yuri Brugnara
#'
#' @examples
#' subdaily_repetition(Rosario$ta, Meta$ta[which(Meta$ta$id=="Rosario"),], n = 3)
#'
#' @export

subdaily_repetition <- function(subdailydata = file.choose(), meta = NULL, 
                                outpath = getwd(), n = 6) {
  
  #Read data and metadata
  if (is.null(dim(subdailydata))) {
    meta <- read_meta(subdailydata, c("id","lat","lon","alt","var","units"))
    subdailydata <- read_sef(subdailydata)
  }
  if (length(meta) != 6) stop("Incorrect metadata (must be a vector of 6 elements)")
  if (dim(subdailydata)[2] != 7) stop("Incorrect dimension of subdailydata")
  meta <- as.character(meta)
  
  #Define flag
  flag <- "subdaily_repetition"
  
  #Flag data
  vector_count <- rle(subdailydata[,7])$lengths
  id_value <- cumsum(vector_count)
  value_flat <- which(vector_count > n-1)
  flat <- c()
  for (i in value_flat) flat <- append(flat, (id_value[i]-vector_count[i]+1):id_value[i])
  if (meta[5] %in% c("rr","sd","fs","sc","sw")) {
    flat <- flat[which(subdailydata[flat,7] != 0)]
  }
  out <- subdailydata[flat, ]
  if (nrow(out) != 0) {
    out$Test <- flag
    namefile <- paste0("qc_", meta[1], "_", meta[5], "_subdaily.txt")
    namefile <- paste(file.path(outpath), namefile, sep="/")
    write_intermediate_subdaily(namefile, out)
  }
  
  return(print("Repetition test completed", quote=FALSE))
}


################################################################################


#' Duplicate times test
#' 
#' Flag times that appear more than once.
#'
#' @param subdailydata A character string giving the path of the input file,
#' or a 7-column matrix with following columns: variable code, year, month,
#' day, hour, minute, value.
#' @param meta A character vector with 6 elements: station ID, latitude, longitude,
#' altitude, variable code, units. If \code{subdailydata} is a path, \code{meta} is
#' ignored.
#' @param outpath Character string giving the path for the QC results. By default
#' this is the working directory.
#' 
#' @details 
#' The input file must follow the Copernicus Station Exchange Format (SEF).
#' 
#' @author Alba Gilabert, Yuri Brugnara
#'
#' @examples 
#' duplicate_times(Bern$p, Meta$p[which(Meta$p$id=="Bern"),])
#'
#' @export

duplicate_times <- function(subdailydata, meta = NULL, outpath = getwd()) {
  
  #Read data and metadata
  if (is.null(dim(subdailydata))) {
    meta <- read_meta(subdailydata, c("id","lat","lon","alt","var","units"))
    subdailydata <- read_sef(subdailydata)
  }
  if (length(meta) != 6) stop("Incorrect metadata (must be a vector of 6 elements)")
  n <- dim(subdailydata)[2]
  if (n != 7) stop("Incorrect dimension of dailydata")
  meta <- as.character(meta)
  subdailydata <- subdailydata[order(subdailydata[,2], subdailydata[,3], subdailydata[,4],
                                     subdailydata[,5], subdailydata[,6]), ]
  
  #Define flag
  flag <- "duplicate_times"
  
  #Find duplicates
  dupli <- duplicated(subdailydata[,2:6]) & !is.na(subdailydata[,5])
  dupli <- unique(c(which(dupli) - 1, which(dupli)))
  out <- subdailydata[dupli, ]
  if (nrow(out) != 0) {
    out$Test <- flag
    namefile <- paste0("qc_", meta[1], "_", meta[5], "_subdaily.txt")
    namefile <- paste(file.path(outpath), namefile, sep="/")
    write_intermediate_subdaily(namefile, out)
  }
  
  return(print("Duplicate times test completed", quote=FALSE))
}
