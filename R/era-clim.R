#' WMO Time Consistency Test for Pressure, Temperature and Dew Point.
#' 
#' Applicable to a series of sub-daily air pressure (p, mslp),  
#' air temperature (ta) or dew point temperature (td) observations with at least
#' some time intervals between observations less or equal to twelve hours. 
#' Flags the records where
#' the observations exceed the WMO suggested tolerances for the temperatures and
#' pressure tendency as function of time period between consecutive reports.
#' 
#' @details 
#' \strong{Input:}
#' \itemize{
#' \item A SEF file or a data frame and metadata. The observations data frame
#' must have seven columns: variable code, year (YYYY), month (MM), day
#' (DD), hour (HH), minute (MM), observation.
#' }
#' \strong{The WMO time consistency test:}
#' \itemize{
#' \item WMO suggested tolerances for the temperatures and pressure tendency as
#' function of time period between consecutive reports (WMO, 1993: VI.21):
#' \tabular{llllll}{Parameter | \tab dt = 1 hour | \tab dt = 2 hours | \tab dt
#' = 3 hours | \tab dt = 6 hours | \tab dt = 12 hours \cr ta_tol \tab 4 ºC \tab
#' 7 ºC  \tab 9 ºC \tab 15 ºC  \tab 25 ºC \cr td_tol \tab 4 ºC \tab 6 ºC \tab 8
#' ºC \tab 12 ºC  \tab 20 ºC \cr pp_tol \tab 3 hPa \tab 6 hPa  \tab 9 hPa \tab 18
#' hPa \tab 36 hPa \cr}
#' \item The temperatures tolerance - \emph{ta_tol} and \emph{td_tol} -
#' considered for 1, 2, 3, 6 and 12 hours is given by the table above.
#' \item The pressure tolerance \emph{p_tol} is determined for time intervals
#' belonging to [1, 12] hours, assuming that there is a linear variation of 3
#' hPa per hour, based in the table above.
#' \item Time consistency test (WMO, 1993: VI.21):
#' 
#' \eqn{| obs(t) - obs(t - dt) | > tol => flag_obs(t) = suspect and
#' flag_obs(t-dt) = suspect}
#' 
#' \eqn{obs - observation (ta, td or p), t - hour in decimal, dt - hour
#' difference in decimal, tol - tolerance}
#' 
#' \item The flag, correspondent to suspect values, is always associated with
#' two consecutive observations within twelve hours.
#' }
#' \strong{Output:}
#' \itemize{
#' \item A text file of flagged observations with eight columns: variable code,
#' year, month, day, hour, minute, value, test.
#' The test column has the description "wmo_time_consistency".
#' }
#' @references 
#' WMO, 1993: Chapter 6 - Quality Control Procedures. Guide on the Global
#' Data-processing System, World Meteorological Organization, Geneva, No. 305,
#' VI.1-VI.27, ISBN 92-63-13305-0.
#' 
#' @param series A character string giving the path of the input file,
#' or a 7-column matrix with following columns: variable code, year, month,
#' day, hour, minute, value.
#' @param meta A character vector with 6 elements: station ID, latitude, longitude,
#' altitude, variable code, units. If \code{series} is a path, \code{meta} is
#' ignored.
#' @param outpath Character string giving the path for the QC results.
#' 
#' @author Clara Ventura, Yuri Brugnara
#' 
#' @examples 
#' wmo_time_consistency(series = Bern$p, meta = Meta$p[which(Meta$p$id=="Bern"),],
#'                      outpath = tempdir())
#' 
#' @export
 
wmo_time_consistency <- function(series, meta = NULL, outpath) {
  
  #Read data and metadata
  if (is.null(dim(series))) {
    meta <- read_meta(series, c("id","lat","lon","alt","var","units"))
    series <- read_sef(series)
  }
  if (length(meta) != 6) stop("Incorrect metadata (must be a vector of 6 elements)")
  n <- dim(series)[2]
  if (n != 7) stop("Incorrect dimension of series")
  names(series) <- c("vcod", "year", "month", "day", "hour", "minut", "obs")
  subd <- series[order(series$year,series$month,series$day,series$hour,series$minut), ]
  idsta <- as.character(meta[1])
  vcode <- as.character(meta[5])
  vunit <- as.character(meta[6])
  
  #Define flag name
  flag <- "wmo_time_consistency"
  
  #Check variable
  if (vcode %in% c("ta","p","td","mslp")) {
    
    #Check units
    subd[,7] <- check_units(subd[,7], vcode, vunit)
    
    # Creates output file name
    out_file <- paste0("qc_", idsta, "_", vcode, "_subdaily.txt")
    out_file <- paste(file.path(outpath), out_file, sep="/")
    
    # Create vectors of differences
    times <- ISOdate(subd[,2], subd[,3], subd[,4], subd[,5], subd[,6])
    dt <- as.numeric(difftime(times[2:length(times)], times[1:(length(times)-1)],
                              units = "hours"))
    dsubd <- abs(subd[2:length(times), 7] - subd[1:(length(times)-1), 7])

    # Tolerable pressure variation: 3 hPa per hour
    # Rule applied until 12 hours as considered in the WMO table
    if (vcode == "p" || vcode == "mslp") {
      tol <- 3 * dt

    # Test applied to {1, 2, 3, 6, 12} hours with the tolerances suggested by WMO,
    # respectively {4, 7, 9, 15, 25} Celsius degrees
    } else if (vcode == "ta") {
      tol <- array(dim = length(dt))
      tol[which(dt <= 1)] <- 4
      tol[which(dt > 1 & dt <= 2)] <- 7
      tol[which(dt > 2 & dt <= 3)] <- 9
      tol[which(dt > 3 & dt <= 4)] <- 11
      tol[which(dt > 4 & dt <= 5)] <- 13
      tol[which(dt > 5 & dt <= 6)] <- 15
      tol[which(dt > 6 & dt <= 7)] <- 17
      tol[which(dt > 7 & dt <= 8)] <- 18
      tol[which(dt > 8 & dt <= 9)] <- 20
      tol[which(dt > 9 & dt <= 10)] <- 22
      tol[which(dt > 10 & dt <= 11)] <- 23
      tol[which(dt > 11 & dt <= 12)] <- 25
      
    # Test applied to {1, 2, 3, 6, 12} hours with the tolerances suggested by WMO,
    # respectively {4, 6, 8, 12, 20} Celsius degrees
    } else if (vcode == "td") {
      tol <- array(dim = length(dt))
      tol[which(dt <= 1)] <- 4
      tol[which(dt > 1 & dt <= 2)] <- 6
      tol[which(dt > 2 & dt <= 3)] <- 8
      tol[which(dt > 3 & dt <= 4)] <- 9
      tol[which(dt > 4 & dt <= 5)] <- 11
      tol[which(dt > 5 & dt <= 6)] <- 12
      tol[which(dt > 6 & dt <= 7)] <- 13
      tol[which(dt > 7 & dt <= 8)] <- 15
      tol[which(dt > 8 & dt <= 9)] <- 16
      tol[which(dt > 9 & dt <= 10)] <- 17
      tol[which(dt > 10 & dt <= 11)] <- 19
      tol[which(dt > 11 & dt <= 12)] <- 20
    } 
    
    # Create data frame of flagged observations
    flags <- which(dt <= 12 & dsubd-tol > 0)
    flags <- append(flags, flags+1)
    err_tc <- subd[flags, ]

    # Output
    if (length(flags) != 0) {
      err_tc <- err_tc[!duplicated(err_tc), ]
      err_tc$Test <- flag
      names(err_tc) <- c("Var","Year","Month","Day","Hour","Minute","Value","Test")
      write_intermediate_subdaily(out_file, err_tc)
    }
    outMsg <- "WMO time consistency test completed"
    
  } else {
    outMsg <- "Variable not supported by this test"
  }
  
  message(outMsg)
  
}


################################################################################


#' Gross Errors Test for Cloud Cover and Relative Humidity.
#' 
#' Applicable to a series (daily or sub-daily) of relative humidity (rh) in percent
#' or to a series of cloud cover (n) in percent or oktas.
#' 
#' @details 
#' \strong{Input:}
#' \itemize{
#' \item A SEF file or a data frame and metadata. The observations data frame
#' must have five or seven columns: variable code, year (YYYY), month (MM), day
#' (DD), (hour (HH), minute (MM)), observation.
#' }
#' \strong{Output:}
#' \itemize{
#' \item A text file of flagged observations with six or eight columns: variable
#' code, year (YYYY), month (MM), day (DD), (hour (HH), minute (MM)),
#' observation, test. The test column has the description "gross_errors".
#' \item The flagged observations correspond to values that don't belong to the
#' integer interval (0, 100) if the unit is percent or that don't belong to the
#' integer interval (0, 9) if the unit is oktas.
#' }
#' 
#' @param series A character string giving the path of the SEF file, or a
#'   five or seven-column (daily or subdaily) data frame with the series.
#' @param meta A character vector with 6 elements: station ID, latitude, longitude,
#' altitude, variable code, units. If \code{series} is a path, \code{meta} is
#' ignored.
#' @param outpath Character string giving the path for the QC results.
#' 
#' @author Clara Ventura, Yuri Brugnara
#' 
#' @examples
#' impossible_values(series = Rosario$n, meta = Meta$n, outpath = tempdir())
#' impossible_values(series = Rosario$rh, meta = Meta$rh, outpath = tempdir())
#' 
#' @export
#'

impossible_values <- function(series, meta = NULL, outpath) {
  
  #Read data and metadata
  if (is.null(dim(series))) {
    meta <- read_meta(series, c("id","lat","lon","alt","var","units"))
    series <- read_sef(series)
  }
  if (length(meta) != 6) stop("Incorrect metadata (must be a vector of 6 elements)")
  n <- dim(series)[2]
  if (!n %in% c(5,7)) stop("Incorrect dimension of series")
  namcols <- c("vcod", "year", "month", "day", "hour", "minut", "obs")
  if (n == 7) {
    names(series) <- namcols
  } else {
    names(series) <- namcols[c(1:4,7)]
  }
  drec <- series
  idsta <- as.character(meta[1])
  vcode <- as.character(meta[5])
  vunit <- as.character(meta[6])
  
  #Define flag name
  flag <- "impossible_values"
  
  #Check variable
  if (vcode %in% c("rh","n")) {
    
    #Check units
    drec[,n] <- check_units(drec[,n], vcode, vunit)
    
    # Creates output file name
    res <- ifelse(n == 5, "daily", "subdaily")
    out_file <- paste0("qc_", idsta, "_", vcode, "_", res, ".txt")
    out_file <- paste(file.path(outpath), out_file, sep="/")
    
    # Test for cloud cover in percent and relative humidity
    if (vunit == "%") {
      err_gel <- drec[!(drec$obs %in% 0:100), ]
      # Test for cloud cover in oktas
    } else if (vunit %in% c("Okta","okta","Oktas","oktas","okt","Okt")) {
      err_gel <- drec[!(drec$obs %in% 0:9), ]
    }
    
    # Output
    if (nrow(err_gel) != 0) {
      err_gel$flag <- flag
      if (n == 7) {
        names(err_gel) <- c("Var","Year","Month","Day","Hour","Minute","Value","Test")
        write_intermediate_subdaily(out_file, err_gel)
      } else {
        names(err_gel) <- c("Var","Year","Month","Day","Value","Test")
        write_intermediate_daily(out_file, err_gel)
      }
    }
    outMsg <- "Impossible values test completed"
    
  } else {
    outMsg <- "Variable not supported by this test"
  }
  
  message(outMsg)
  
}


################################################################################


#' WMO Gross Errors Tests for Pressure, Temperature, Dew Point, and Wind Speed.
#' 
#' Applicable to a series (daily or sub-daily) of air pressure, air temperature (ta),
#' dew point temperature (td), wind speed (w).
#' The pressure series can be at mean sea level (mslp) or at station level (p). 
#' Flags the records where the observations values exceed the limit values given 
#' by WMO (1993).
#' 
#' @details 
#' \strong{Input:}
#' \itemize{
#' \item A SEF file or a data frame and metadata. The observations data frame
#' must have five or seven columns: variable code, year (YYYY), month (MM), day
#' (DD), (hour (HH), minute (MM)), observation. The required metadata are the
#' station identifier, the station latitude and variable units.
#' }
#' \strong{The WMO gross error limits for air pressure, air temperature, dew
#' point temperature, and wind speed:}
#' \itemize{
#' \item For station level pressure the gross error limits are latitude and
#' meteorological season independent (WMO, 1993: VI.7). According to the same
#' reference, for mean sea level pressure, temperature, dew point, and wind speed,
#' the WMO establishes the gross error limits as
#' function of the station latitude and the meteorological season in which the
#' observations were collected.
#' \item The tests divide the meteorological seasons in Winter and Summer. So,
#' based on the meteorological calendar for the Northern Hemisphere, which
#' defines seasons as Spring (March, April, May), Summer (June, July, August),
#' Autumn (September, October, November) and Winter (December, January,
#' February), it was here considered:
#' \itemize{
#' \item Northern Hemisphere Winter / Southern Hemisphere Summer - January,
#' February, March, October, November, December;
#' \item Northern Hemisphere Summer / Southern Hemisphere Winter  - April, May,
#' June, July, August, September.
#' }
#' \item The gross error limits for each variable divide the flagged values in
#' suspect and erroneous (WMO, 1993: VI.6 - VI.8).
#' \item \strong{Latitude independent}
#' \itemize{
#' \item \strong{Meteorological season independent}
#' \itemize{
#' \item \strong{Station Level Pressure (p):} 
#' \itemize{
#' \item Suspect: \emph{300 <= p < 400 hPa or 1080 < p <= 1100 hPa}
#' \item Erroneous: \emph{p < 300 or p > 1100 hPa}
#' }
#' }
#' }
#' \item \strong{Latitudes belonging to the interval [-45, +45]}
#' \itemize{
#' \item \strong{Winter}
#' \itemize{
#' \item \strong{Mean Sea Level Pressure (mslp)}
#' \itemize{
#' \item Suspect: \emph{870 <= mslp < 910 hPa or 1080 < mslp <= 1100 hPa}
#' \item Erroneous: \emph{mslp < 870 hPa or mslp > 1100 hPa}
#' }
#' \item \strong{Air Temperature (ta)}
#' \itemize{
#' \item Suspect: \emph{-40 <= ta < -30 ºC or 50 < ta <= 55 ºC}
#' \item Erroneous: \emph{ta < -40 ºC or ta > 55 ºC}
#' }
#' \item \strong{Dew Point Temperature (td)}
#' \itemize{
#' \item Suspect: \emph{-45 <= td < -35 ºC or 35 < td <= 40 ºC}
#' \item Erroneous: \emph{td < -45 ºC or td > 40 ºC}
#' }
#' \item \strong{Wind Speed (w)}
#' \itemize{
#' \item Suspect: \emph{w > 60 m/s and w <= 125 m/s}
#' \item Erroneous: \emph{w > 125 m/s}
#' }
#' }
#' \item \strong{Summer}
#' \itemize{
#' \item \strong{Mean Sea Level Pressure (mslp)}
#' \itemize{
#' \item Suspect: \emph{850 <= mslp < 900 hPa or 1080 < mslp <= 1100 hPa}
#' \item Erroneous: \emph{mslp < 850 hPa or mslp > 1100 hPa}
#' }
#' \item \strong{Air Temperature (ta)}
#' \itemize{
#' \item Suspect: \emph{-30 <= ta < -20 ºC or 50 < ta <= 60 ºC}
#' \item Erroneous: \emph{ta < -30 ºC or ta > 60 ºC}
#' }
#' \item \strong{Dew Point Temperature (td)}
#' \itemize{
#' \item Suspect: \emph{-35 <= td < -25 ºC or 35 < td <= 40 ºC}
#' \item Erroneous: \emph{td < -35 ºC or td > 40 ºC}
#' }
#' \item \strong{Wind Speed (w)}
#' \itemize{
#' \item Suspect: \emph{w > 90 m/s and w <= 150 m/s}
#' \item Erroneous: \emph{w > 150 m/s}
#' }
#' }
#' }
#' \item \strong{Latitudes belonging to the interval [-90, -45[ U ]+45, +90]}
#' \itemize{
#' \item \strong{Winter}
#' \itemize{
#' \item \strong{Mean Sea Level Pressure (mslp)}
#' \itemize{
#' \item Suspect: \emph{910 <= mslp < 940 hPa or 1080 < mslp <= 1100 hPa}
#' \item Erroneous: \emph{mslp < 910 hPa or mslp > 1100 hPa}
#' }
#' \item \strong{Air Temperature (ta)}
#' \itemize{
#' \item Suspect: \emph{-90 <= ta < -80 ºC or 35 < ta <= 40 ºC}
#' \item Erroneous: \emph{ta < -90 ºC or ta > 40 ºC}
#' }
#' \item \strong{Dew Point Temperature (td)}
#' \itemize{
#' \item Suspect: \emph{-99 <= td < -85 ºC or 30 < td <= 35 ºC}
#' \item Erroneous: \emph{td < -99 ºC or td > 35 ºC}
#' }
#' \item \strong{Wind Speed (w)}
#' \itemize{
#' \item Suspect: \emph{w > 50 m/s and w <= 100 m/s}
#' \item Erroneous: \emph{w > 100 m/s}
#' }
#' }
#' \item \strong{Summer}
#' \itemize{
#' \item \strong{Mean Sea Level Pressure (mslp)}
#' \itemize{
#' \item Suspect: \emph{920 <= mslp < 950 hPa or 1080 < mslp <= 1100 hPa}
#' \item Erroneous: \emph{mslp < 920 hPa or mslp > 1100 hPa}
#' }
#' \item \strong{Air Temperature (ta)}
#' \itemize{
#' \item Suspect: \emph{-40 <= ta < -30 ºC or 40 < ta <= 50 ºC}
#' \item Erroneous: \emph{ta < -40 ºC or ta > 50 ºC}
#' }
#' \item \strong{Dew Point Temperature (td)}
#' \itemize{
#' \item Suspect: \emph{-45 <= td < -35 ºC or 35 < td <= 40 ºC}
#' \item Erroneous: \emph{td < -45 ºC or td > 40 ºC}
#' }
#' \item \strong{Wind Speed (w)}
#' \itemize{
#' \item Suspect values: \emph{w > 40 m/s and w <= 75 m/s}
#' \item Erroneous values: \emph{w > 75 m/s}
#' }
#' }
#' }
#' }
#' \strong{Output:}
#' \itemize{
#' \item A text file of flagged observations with six or eight columns: variable
#' code, year (YYYY), month (MM), day (DD), (hour (HH), minute (MM)),
#' observation, test. The test column has the description "gross_errors".
#' }
#' 
#' @references 
#' WMO, 1993: Chapter 6 - Quality Control Procedures. Guide on the Global
#' Data-processing System, World Meteorological Organization, Geneva, No. 305,
#' VI.1-VI.27, ISBN 92-63-13305-0.
#' 
#' @param series A character string giving the path of the SEF file, or a
#' five or seven-column (daily or subdaily) data frame with the series.
#' @param meta A character vector with 6 elements: station ID, latitude, longitude,
#' altitude, variable code, units. If \code{series} is a path, \code{meta} is
#' ignored.
#' @param outpath Character string giving the path for the QC results.
#' 
#' @author Clara Ventura, Yuri Brugnara
#' 
#' @examples
#' wmo_gross_errors(series = Rosario$p, meta = Meta$p[which(Meta$p$id=="Rosario"),],
#'                  outpath = tempdir())
#' wmo_gross_errors(series = Rosario$ta, meta = Meta$ta[which(Meta$p$id=="Rosario"),],
#'                  outpath = tempdir())
#' wmo_gross_errors(series = Rosario$td, meta = Meta$td, outpath = tempdir())
#' 
#' @export

wmo_gross_errors <- function(series, meta = NULL, outpath) {
  
  #Read data and metadata
  if (is.null(dim(series))) {
    meta <- read_meta(series, c("id","lat","lon","alt","var","units"))
    series <- read_sef(series)
  }
  if (length(meta) != 6) stop("Incorrect metadata (must be a vector of 6 elements)")
  n <- dim(series)[2]
  if (!n %in% c(5,7)) stop("Incorrect dimension of series")
  namcols <- c("vcod", "year", "month", "day", "hour", "minut", "obs")
  if (n == 7) {
    names(series) <- namcols
  } else {
    names(series) <- namcols[c(1:4,7)]
  }
  drec <- series
  idsta <- as.character(meta[1])
  lat <- as.numeric(meta[2])
  vcode <- as.character(meta[5])
  vunit <- as.character(meta[6])
  
  #Define flag name
  flag <- "wmo_gross_errors"
  
  #Check variable
  if (vcode %in% c("mslp","p","ta","td","w")) {
    
    #Check units
    drec[,n] <- check_units(drec[,n], vcode, vunit)
    
    # Creates output file name
    res <- ifelse(n == 5, "daily", "subdaily")
    out_file <- paste0("qc_", idsta, "_", vcode, "_", res, ".txt")
    out_file <- paste(file.path(outpath), out_file, sep="/")
    
    # Separates data collected in the Winter from data collected in the Summer
    # Data collected in Northern Hemisphere Winter / Southern Hemisphere Summer
    new_ses <- drec[drec$month == 1 | drec$month == 2 | drec$month == 3 |
                      drec$month == 10 | drec$month == 11 | drec$month == 12, ]
    # Data collected in Northern Hemisphere Summer / Southern Hemisphere Winter
    nes_sew <- drec[drec$month == 4 | drec$month == 5 | drec$month == 6 |
                      drec$month == 7 | drec$month == 8 | drec$month == 9, ]
    # Gross errors tests
    if (vcode == "p" || vcode == "mslp") {
      # err_gel <- gel_pressure(drec, new_ses, nes_sew, vcode)
      if (vcode == "mslp") {
        # Latitudes belonging to the interval [-45, +45]
        # Southern Hemisphere
        if (lat >= -45 && lat < 0) {
          # Southern Hemisphere Winter
          # Suspect values if: 870 <= p < 910 hPa or 1080 < p <= 1100 hPa
          wsusp <- nes_sew[(nes_sew$obs >= 870 & nes_sew$obs < 910) | 
                             (nes_sew$obs > 1080 & nes_sew$obs <= 1100), ]
          # Erroneous if: p < 870 hPa or p > 1100 hPa
          werro <- nes_sew[nes_sew$obs < 870 | nes_sew$obs > 1100, ]
          # Southern Hemisphere Summer
          # Suspect values if: 850 <= p < 900 hPa or 1080 < p <= 1100 hPa
          ssusp <- new_ses[(new_ses$obs >= 850 & new_ses$obs < 900) | 
                             (new_ses$obs > 1080 & new_ses$obs <= 1100), ]
          # Erroneous if: p < 850 hPa or p > 1100 hPa
          serro <- new_ses[new_ses$obs < 850 | new_ses$obs > 1100, ]
          # Northern Hemisphere
        } else if (lat >= 0 && lat <= 45) {
          # Northern Hemisphere Winter
          # Suspect values if: 870 <= p < 910 hPa or 1080 < p <= 1100 hPa
          wsusp <- new_ses[(new_ses$obs >= 870 & new_ses$obs < 910) | 
                             (new_ses$obs > 1080 & new_ses$obs <= 1100), ]
          # Erroneous if: p < 870 hPa or p > 1100 hPa
          werro <- new_ses[new_ses$obs < 870 | new_ses$obs > 1100, ]
          # Northern Hemisphere Summer
          # Suspect values if: 850 <= p < 900 hPa or 1080 < p <= 1100 hPa
          ssusp <- nes_sew[(nes_sew$obs >= 850 & nes_sew$obs < 900) | 
                             (nes_sew$obs > 1080 & nes_sew$obs <= 1100), ]
          # Erroneous if: p < 850 hPa or p > 1100 hPa
          serro <- nes_sew[nes_sew$obs < 850 | nes_sew$obs > 1100, ]
          # Latitudes belonging to the interval [-90, -45[ U ]+45, +90]
          # Southern Hemisphere
        } else if (lat >= -90 && lat < -45) {
          # Southern Hemisphere Winter
          # Suspect values if: 910 <= p < 940 hPa or 1080 < p <= 1100 hPa
          wsusp <- nes_sew[(nes_sew$obs >= 910 & nes_sew$obs < 940) | 
                             (nes_sew$obs > 1080 & nes_sew$obs <= 1100), ]
          # Erroneous if: p < 910 hPa or p > 1100 hPa
          werro <- nes_sew[nes_sew$obs < 910 | nes_sew$obs > 1100, ]
          # Southern Hemisphere Summer
          # Suspect values if: 920 <= p < 950 hPa or 1080 < p <= 1100 hPa
          ssusp <- new_ses[(new_ses$obs >= 920 & new_ses$obs < 950) | 
                             (new_ses$obs > 1080 & new_ses$obs <= 1100), ]
          # Erroneous if: p < 920 hPa or p > 1100 hPa
          serro <- new_ses[new_ses$obs < 920 | new_ses$obs > 1100, ]
          # Northern Hemisphere  
        } else if (lat > 45 && lat <= 90) {
          # Northern Hemisphere Winter
          # Suspect values if: 910 <= p < 940 hPa or 1080 < p <= 1100 hPa
          wsusp <- new_ses[(new_ses$obs >= 910 & new_ses$obs < 940) | 
                             (new_ses$obs > 1080 & new_ses$obs <= 1100), ]
          # Erroneous if: p < 910 hPa or p > 1100 hPa
          werro <- new_ses[new_ses$obs < 910 | new_ses$obs > 1100, ]
          # Northern Hemisphere Summer
          # Suspect values if: 920 <= p < 950 hPa or 1080 < p <= 1100 hPa
          ssusp <- nes_sew[(nes_sew$obs >= 920 & nes_sew$obs < 950) | 
                             (nes_sew$obs > 1080 & nes_sew$obs <= 1100), ]
          # Erroneous if: p < 920 hPa or p > 1100 hPa
          serro <- nes_sew[nes_sew$obs < 920 | nes_sew$obs > 1100, ]
        }
        err_gel <- na.omit(rbind(wsusp, werro, ssusp, serro))
      } else if (vcode == "p") {
        # Suspect values if: 300 <= p < 400 hPa or 1080 < p <= 1100 hPa
        susp <- drec[(drec$obs >= 300 & drec$obs < 400) | 
                       (drec$obs > 1080 & drec$obs <= 1100), ]
        # Erroneous if: p < 300 or p > 1100 hPa
        erro <- drec[drec$obs < 300 | drec$obs > 1100, ]
        err_gel <- na.omit(rbind(susp, erro))
      }
    } else if (vcode == "ta") {
      # err_gel <- gel_airtemp(new_ses, nes_sew)
      # Latitudes belonging to the interval [-45, +45]
      # Southern Hemisphere
      if (lat >= -45 && lat < 0) {
        # Southern Hemisphere Winter
        # Suspect values if: -40 <= ta < -30 ºC or 50 < ta <= 55 ºC
        wsusp <- nes_sew[(nes_sew$obs >= -40 & nes_sew$obs < -30) | 
                           (nes_sew$obs > 50 & nes_sew$obs <= 55), ]
        # Erroneous if: ta < -40 ºC or ta > 55 ºC
        werro <- nes_sew[nes_sew$obs < -40 | nes_sew$obs > 55, ]
        # Southern Hemisphere Summer
        # Suspect values if: -30 <= ta < -20 ºC or 50 < ta <= 60 ºC
        ssusp <- new_ses[(new_ses$obs >= -30 & new_ses$obs < -20) | 
                           (new_ses$obs > 50 & new_ses$obs <= 60), ]
        # Erroneous if: ta < -30 ºC or ta > 60 ºC
        serro <- new_ses[new_ses$obs < -30 | new_ses$obs > 60, ]
        # Northern Hemisphere
      } else if (lat >= 0 && lat <= 45) {
        # Northern Hemisphere Winter
        # Suspect values if: -40 <= ta < -30 ºC or 50 < ta <= 55 ºC
        wsusp <- new_ses[(new_ses$obs >= -40 & new_ses$obs < -30) | 
                           (new_ses$obs > 50 & new_ses$obs <= 55), ]
        # Erroneous if: ta < -40 ºC or ta > 55 ºC
        werro <- new_ses[new_ses$obs < -40 | new_ses$obs > 55, ]
        # Northern Hemisphere Summer
        # Suspect values if: -30 <= ta < -20 ºC or 50 < ta <= 60 ºC
        ssusp <- nes_sew[(nes_sew$obs >= -30 & nes_sew$obs < -20) | 
                           (nes_sew$obs > 50 & nes_sew$obs <= 60), ]
        # Erroneous if: ta < -30 ºC or ta > 60 ºC
        serro <- nes_sew[nes_sew$obs < -30 | nes_sew$obs > 60, ]
        # Latitudes belonging to the interval [-90, -45[ U ]+45, +90]
        # Southern Hemisphere
      } else if (lat >= -90 && lat < -45) {
        # Southern Hemisphere Winter
        # Suspect values if: -90 <= ta < -80 ºC or 35 < ta <= 40 ºC
        wsusp <- nes_sew[(nes_sew$obs >= -90 & nes_sew$obs < -80) | 
                           (nes_sew$obs > 35 & nes_sew$obs <= 40), ]
        # Erroneous if: ta < -90 ºC or ta > 40 ºC
        werro <- nes_sew[nes_sew$obs < -90 | nes_sew$obs > 40, ]
        # Southern Hemisphere Summer
        # Suspect values if: -40 <= ta < -30 ºC or 40 < ta <= 50 ºC
        ssusp <- new_ses[(new_ses$obs >= -40 & new_ses$obs < -30) | 
                           (new_ses$obs > 40 & new_ses$obs <= 50), ]
        # Erroneous if: ta < -40 ºC or ta > 50 ºC
        serro <- new_ses[new_ses$obs < -40 | new_ses$obs > 50, ]
        # Northern Hemisphere  
      } else if (lat > 45 && lat <= 90) {
        # Northern Hemisphere Winter
        # Suspect values if: -90 <= ta < -80 ºC or 35 < ta <= 40 ºC
        wsusp <- new_ses[(new_ses$obs >= -90 & new_ses$obs < -80) | 
                           (new_ses$obs > 35 & new_ses$obs <= 40), ]
        # Erroneous if: ta < -90 ºC or ta > 40 ºC
        werro <- new_ses[new_ses$obs < -90 | new_ses$obs > 40, ]
        # Northern Hemisphere Summer
        # Suspect values if: -40 <= ta < -30 ºC or 40 < ta <= 50 ºC
        ssusp <- nes_sew[(nes_sew$obs >= -40 & nes_sew$obs < -30) | 
                           (nes_sew$obs > 40 & nes_sew$obs <= 50), ]
        # Erroneous if: -40 ºC or ta > 50 ºC
        serro <- nes_sew[nes_sew$obs < -40 | nes_sew$obs > 50, ]
      }
      err_gel <- na.omit(rbind(wsusp, werro, ssusp, serro))
    } else if (vcode == "td") {
      # err_gel <- gel_dewpoint(new_ses, nes_sew)
      if (lat >= -45 && lat < 0) {
        # Southern Hemisphere Winter
        # Suspect values if: -45 <= td < -35 ºC or 35 < td <= 40 ºC
        wsusp <- nes_sew[(nes_sew$obs >= -45 & nes_sew$obs < -35) | 
                           (nes_sew$obs > 35 & nes_sew$obs <= 40), ]
        # Erroneous if: td < -45 ºC or td > 40 ºC
        werro <- nes_sew[nes_sew$obs < -45 | nes_sew$obs > 40, ]
        # Southern Hemisphere Summer
        # Suspect values if: -35 <= td < -25 ºC or 35 < td <= 40 ºC
        ssusp <- new_ses[(new_ses$obs >= -35 & new_ses$obs < -25) | 
                           (new_ses$obs > 35 & new_ses$obs <= 40), ]
        # Erroneous if: td < -35 ºC or td > 40 ºC
        serro <- new_ses[new_ses$obs < -35 | new_ses$obs > 40, ]
        # Northern Hemisphere
      } else if (lat >= 0 && lat <= 45) {
        # Northern Hemisphere Winter
        # Suspect values if: -45 <= td < -35 ºC or 35 < td <= 40 ºC
        wsusp <- new_ses[(new_ses$obs >= -45 & new_ses$obs < -35) | 
                           (new_ses$obs > 35 & new_ses$obs <= 40), ]
        # Erroneous if: td < -45 ºC or td > 40 ºC
        werro <- new_ses[new_ses$obs < -45 | new_ses$obs > 40, ]
        # Northern Hemisphere Summer
        # Suspect values if: -35 <= td < -25 ºC or 35 < td <= 40 ºC
        ssusp <- nes_sew[(nes_sew$obs >= -35 & nes_sew$obs < -25) | 
                           (nes_sew$obs > 35 & nes_sew$obs <= 40), ]
        # Erroneous if: td < -35 ºC or td > 40 ºC
        serro <- nes_sew[nes_sew$obs < -35 | nes_sew$obs > 40, ]
        # Latitudes belonging to the interval [-90, -45[ U ]+45, +90]
        # Southern Hemisphere
      } else if (lat >= -90 && lat < -45) {
        # Southern Hemisphere Winter
        # Suspect values if: -99 <= td < -85 ºC or 30 < td <= 35 ºC
        wsusp <- nes_sew[(nes_sew$obs >= -99 & nes_sew$obs < -85) | 
                           (nes_sew$obs > 30 & nes_sew$obs <= 35), ]
        # Erroneous if: td < -99 ºC or td > 35 ºC
        werro <- nes_sew[nes_sew$obs < -99 | nes_sew$obs > 35, ]
        # Southern Hemisphere Summer
        # Suspect values if: -45 <= td < -35 ºC or 35 < td <= 40 ºC
        ssusp <- new_ses[(new_ses$obs >= -45 & new_ses$obs < -35) | 
                           (new_ses$obs > 35 & new_ses$obs <= 40), ]
        # Erroneous if: td < -45 ºC or td > 40 ºC
        serro <- new_ses[new_ses$obs < -45 | new_ses$obs > 40, ]
        # Northern Hemisphere  
      } else if (lat > 45 && lat <= 90) {
        # Northern Hemisphere Winter
        # Suspect values if: -99 <= td < -85 ºC or 30 < td <= 35 ºC
        wsusp <- new_ses[(new_ses$obs >= -99 & new_ses$obs < -85) | 
                           (new_ses$obs > 30 & new_ses$obs <= 35), ]
        # Erroneous if: td < -99 ºC or td > 35 ºC
        werro <- new_ses[new_ses$obs < -99 | new_ses$obs > 35, ]
        # Northern Hemisphere Summer
        # Suspect values if: -45 <= td < -35 ºC or 35 < td <= 40 ºC
        ssusp <- nes_sew[(nes_sew$obs >= -45 & nes_sew$obs < -35) | 
                           (nes_sew$obs > 35 & nes_sew$obs <= 40), ]
        # Erroneous if: td < -45 ºC or td > 40 ºC
        serro <- nes_sew[nes_sew$obs < -45 | nes_sew$obs > 40, ]
      }
      err_gel <- na.omit(rbind(wsusp, werro, ssusp, serro))
    } else if (vcode == "w") {
      # err_gel <- gel_windspeed(new_ses, nes_sew)
      # Latitudes belonging to the interval [-45, +45]
      # Southern Hemisphere
      if (lat >= -45 && lat < 0) {
        # Southern Hemisphere Winter
        # Suspect values if: w > 60 m/s and w <= 125 m/s
        wsusp <- nes_sew[(nes_sew$obs > 60 & nes_sew$obs <= 125), ]
        # Erroneous if: w > 125 m/s
        werro <- nes_sew[nes_sew$obs > 125, ]
        # Southern Hemisphere Summer
        # Suspect values if: w > 90 m/s and w <= 150 m/s
        ssusp <- new_ses[(new_ses$obs > 90 & new_ses$obs <= 150), ]
        # Erroneous if: w > 150 m/s
        serro <- new_ses[new_ses$obs > 150, ]
        # Northern Hemisphere
      } else if (lat >= 0 && lat <= 45) {
        # Northern Hemisphere Winter
        # Suspect values if: w > 60 m/s and w <= 125 m/s
        wsusp <- new_ses[(new_ses$obs > 60 & new_ses$obs <= 125), ]
        # Erroneous if: w > 125 m/s
        werro <- new_ses[new_ses$obs > 125, ]
        # Northern Hemisphere Summer
        # Suspect values if: w > 90 m/s and w <= 150 m/s
        ssusp <- nes_sew[(nes_sew$obs > 90 & nes_sew$obs <= 150), ]
        # Erroneous if: w > 150 m/s
        serro <- nes_sew[nes_sew$obs > 150, ]
        # Latitudes belonging to the interval [-90, -45[ U ]+45, +90]
        # Southern Hemisphere
      } else if (lat >= -90 && lat < -45) {
        # Southern Hemisphere Winter
        # Suspect values if: w > 50 m/s and w <= 100 m/s
        wsusp <- nes_sew[(nes_sew$obs > 50 & nes_sew$obs <= 100), ]
        # Erroneous if: w > 100 m/s
        werro <- nes_sew[nes_sew$obs > 100, ]
        # Southern Hemisphere Summer
        # Suspect values if: w > 40 m/s and w <= 75 m/s
        ssusp <- new_ses[(new_ses$obs > 40 & new_ses$obs <= 75), ]
        # Erroneous if: w > 75 m/s
        serro <- new_ses[new_ses$obs > 75, ]
        # Northern Hemisphere  
      } else if (lat > 45 && lat <= 90) {
        # Northern Hemisphere Winter
        # Suspect values if: w > 50 m/s and w <= 100 m/s
        wsusp <- new_ses[(new_ses$obs > 50 & new_ses$obs <= 100), ]
        # Erroneous if: w > 100 m/s
        werro <- new_ses[new_ses$obs > 100, ]
        # Northern Hemisphere Summer
        # Suspect values if: w > 40 m/s and w <= 75 m/s
        ssusp <- nes_sew[(nes_sew$obs > 40 & nes_sew$obs <= 75), ]
        # Erroneous if: w > 75 m/s
        serro <- nes_sew[nes_sew$obs > 75, ]
      }
      err_gel <- na.omit(rbind(wsusp, werro, ssusp, serro))
    }
    
    # Output
    if (nrow(err_gel) != 0) {
#      if (!is.null(err_gel$obs_ou)) {
#        err_gel$obs <- err_gel$obs_ou
#        err_gel$obs_ou <- NULL
#      }
      err_gel$flag <- flag
      if (n == 7) {
        names(err_gel) <- c("Var","Year","Month","Day","Hour","Minute","Value","Test")
        write_intermediate_subdaily(out_file, err_gel)
      } else {
        names(err_gel) <- c("Var","Year","Month","Day","Value","Test")
        write_intermediate_daily(out_file, err_gel)
      }
    }
    outMsg <- "WMO gross errors test completed"
    
  } else {
    outMsg <- "Variable not supported by this test"
  }
  
  message(outMsg)
  
}