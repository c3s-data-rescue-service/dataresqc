#' Convert pressure data to hPa
#' 
#' Converts pressure observations made with a mercury barometer to SI units.
#' If geographical coordinates are given, a gravity correction is applied.
#' If attached temperature is given, a temperature correction is applied.
#'
#' @param p A numerical vector of barometer observations in any unit of length.
#' @param f Conversion factor to mm (e.g., 25.4 for English inches).
#' @param lat Station latitude (degrees North in decimal).
#' @param alt Station altitude (metres). Assumed zero if not given. Ignored if
#' \code{lat} is NA.
#' @param atb A vector of the attached temperature observations in Celsius.
#' 
#' @return
#' A numerical vector of pressure values in hPa.
#' 
#' @author Yuri Brugnara
#' 
#' @references
#' WMO, 2008: Guide to meteorological instruments and methods of observation,
#' WMO-No. 8, World Meteorological Organization, Geneva.
#' 
#' @examples
#' convert_pressure(760) # Gives a standard pressure of 1013.25 hPa
#' 
#' convert_pressure(760, lat=70, alt=100) # Gives a higher pressure because of higher g
#' 
#' convert_pressure(760, lat=70, alt=100, atb=20) # Gives a lower pressure because the
#'                                                # temperature correction is larger than 
#'                                                # the gravity correction.
#' 
#' @export

convert_pressure <- function(p, f = 1, lat = NA, alt = NA, atb = NULL) {
  
  rho <- 13595.1
  gamma <- 0.000182
  lat <- lat * pi / 180
  
  if(is.na(alt)) alt <- 0
  
  g <- ifelse(!is.na(lat),
              9.80620 * (1 - 0.0026442 * cos(2*lat) - 0.0000058 * (cos(2*lat))**2) 
              - 0.000003086 * alt,
              9.80665)
  
  if (!is.null(atb)) {
    if (length(atb) != length(p)) stop("atb must have the same length of p")
    p <- (1 - gamma * atb) * p
  }
  
  p <- p * f * rho * g * 1e-05
  
  return(p)
  
}


###############################################################################


## Convert data to the units used by the QC functions.

check_units <- function(x, v, u) {
  
  ## Define error message for unknown units
  error_msg <- paste0("Unknown units for ", v, ": ", u)
  
  ## Temperature
  if (v %in% c("ta","tb","td","t_air","t_wet","t_dew","Tx","Tn","dep_dew","ibt",
               "atb","Txs","TGs","Tns","TGn","t_snow","Ts","t_water")) {
    if (u %in% c("C","K","F","R","c","k","f","r")) {
      if (u %in% c("K","k")) {
        x <- round(x - 273.15, 1)
      } else if (u %in% c("F","f")) {
        x <- round((x - 32) * 5 / 9, 1)
      } else if (u %in% c("R", "r")) {
        x <- round(x * 1.25, 1)
      }
    } else {
      stop(error_msg)
    }
  }
  
  ## Pressure
  if (v %in% c("p","mslp","pppp")) {
    if (u %in% c("hPa","Pa","mm","mmHg","in","hpa","pa","mmhg","\"")) {
      if (u %in% c("Pa","pa")) {
        x <- round(x / 100, 1)
      } else if (u %in% c("mm","mmHg","mmhg")) {
        x <- round(x * 1013.25 / 760, 1)
      } else if (u %in% c("in", "\"")) {
        x <- round(x * 25.4 * 1013.25 / 760, 1)
      }
    } else {
      stop(error_msg)
    }
  }
  
  ## Precipitation
  if (v %in% c("rr","sw","rrls")) {
    if (u %in% c("mm","in","\"")) {
      if (u %in% c("in","\"")) {
        x <- round(x * 25.4, 1)
      }
    } else {
      stop(error_msg)
    }
  }
  
  ## Snow
  if (v %in% c("sd","fs")) {
    if (u %in% c("cm","mm","m","in","\"","ft")) {
      if (u %in% c("in","\"")) {
        x <- round(x * 2.54, 1)
      } else if (u == "mm") {
        x <- round(x / 10, 1)
      } else if (u == "m") {
        x <- x * 100
      } else if (u == "ft") {
        x <- round(x * 30.48, 1)
      }
    } else {
      stop(error_msg)
    }
  }
  
  ## Wind
  if (v == "w") {
    if (u %in% c("m/s","mps","km/h","kph","mph","kn","kt")) {
      if (u %in% c("km/h","kph")) {
        x <- round(x / 3.6, 1)
      } else if (u == "mph") {
        x <- round(x / 2.2369, 1)
      } else if (u %in% c("kn","kt")) {
        x <- round(x / 1.9438, 1)
      }
    } else {
      stop(error_msg)
    }
  }
  
  ## Humidity
  if (v == "rh") {
    if (!u %in% c("%","perc","percent")) {
      stop(error_msg)
    }
  }
  
  ## Cloud cover
  if (v == "n") {
    if (!u %in% c("%","perc","percent","Okta","okta","Oktas","oktas","okt","Okt")) {
      stop(error_msg)
    }
  }
  
  return(x)
  
}


################################################################################


#' Download a GHCN-Daily data file from the Climate Explorer and convert it
#' into the Station Exchange Format
#'
#' @param url Character string giving the url of the data file.
#' @param outpath Character string giving the path where to save the file.
#'
#' @author Yuri Brugnara
#'
#' @import utils
#' @export


climexp_to_sef <- function(url, outpath) {
  
  filename <- rev(strsplit(url, "/")[[1]])[1]
  download.file(url, paste(file.path(outpath), filename, sep = "/"))
  
  ## Find out variable, units, statistic
  v <- substr(filename, 1, 1)
  if (v == "x") {
    vbl <- "Tx"
    u <- "C"
    st <- "maximum"
  } else if (v == "n") {
    vbl <- "Tn"
    u <- "C"
    st <- "minimum"
  } else if (v == "v") {
    vbl <- "ta"
    u <- "C"
    st <- "mean"
  } else if (v == "p") {
    vbl <- "rr"
    u <- "mm"
    st <- "sum"
  } else if (v == "f") {
    vbl <- "fs"
    u <- "mm"
    st <- "sum"
  } else if (v == "d") {
    vbl <- "sd"
    u <- "mm"
    st <- "point"
  }
  
  ## Read data
  Data <- read.table(paste(file.path(outpath), filename, sep = "/"))
  names(Data) <- c("Year", "Month", "Day", "Value")
  Data$Hour <- ""
  Data$Minute <- ""
  Data <- Data[, c("Year", "Month", "Day", "Hour", "Minute", "Value")]
  
  ## Read metadata
  header <- read.table(paste(file.path(outpath), filename, sep = "/"), 
                       comment.char = "", nrows = 2, sep = "\t")[2, 1]
  lat <- substr(header, 16, 21)
  lon <- substr(header, 25, 31)
  alt <- substr(header, 37, 42)
  
  ## Remove original file
  file.remove(paste(file.path(outpath), filename, sep = "/"))
  
  ## Write SEF file
  write_sef(Data, outpath,
            variable = vbl,
            cod = substr(header, 67, 77),
            nam = substr(header, 79, 108),
            lat = substr(header, 16, 21),
            lon = substr(header, 25, 31),
            alt = substr(header, 37, 42),
            sou = "GHCN-D",
            units = u,
            stat = st,
            metaHead = "Data policy=U.S. Government Work (non-commercial)",
            period = ifelse(vbl == "sd", 0, "day"))
  
}