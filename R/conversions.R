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
        x <- x - 273.15
      } else if (u %in% c("F","f")) {
        x <- (x - 32) * 5 / 9
      } else if (u %in% c("R", "r")) {
        x <- x * 1.25
      }
    } else {
      stop(error_msg)
    }
  }
  
  ## Pressure
  if (v %in% c("p","mslp","pppp")) {
    if (u %in% c("hPa","Pa","mm","mmHg","in","hpa","pa","mmhg","\"")) {
      if (u %in% c("Pa","pa")) {
        x <- x / 100
      } else if (u %in% c("mm","mmHg","mmhg")) {
        x <- x * 1013.25 / 760
      } else if (u %in% c("in", "\"")) {
        x <- x * 25.4 * 1013.25 / 760
      }
    } else {
      stop(error_msg)
    }
  }
  
  ## Precipitation
  if (v %in% c("rr","sw","rrls")) {
    if (u %in% c("mm","in","\"")) {
      if (u %in% c("in","\"")) {
        x <- x * 25.4
      }
    } else {
      stop(error_msg)
    }
  }
  
  ## Snow
  if (v %in% c("sd","fs")) {
    if (u %in% c("cm","mm","m","in","\"","ft")) {
      if (u %in% c("in","\"")) {
        x <- x * 2.54
      } else if (u == "mm") {
        x <- x / 10
      } else if (u == "m") {
        x <- x * 100
      } else if (u == "ft") {
        x <- x * 30.48
      }
    } else {
      stop(error_msg)
    }
  }
  
  ## Wind
  if (v == "w") {
    if (u %in% c("m/s","mps","km/h","kph","mph","kn","kt")) {
      if (u %in% c("km/h","kph")) {
        x <- x / 3.6 
      } else if (u == "mph") {
        x <- x / 2.2369 
      } else if (u %in% c("kn","kt")) {
        x <- x / 1.9438
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