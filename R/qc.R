#' Apply all tests
#' 
#' Perform all quality tests at once on multiple stations and multiple
#' variables.
#' 
#' @param Data Either a character vector of paths to SEF files, a data frame
#' or a list of data frames with 7 columns (one data frame for each station):
#' variable code, year, month, day, hour, minute, value.
#' Each data frame can contain more than one variable code.
#' @param Metadata A data frame with 7 columns: station ID, latitude, longitude,
#' altitude, variable code, units, resolution. If \code{Data} is a list, the station 
#' IDs must appear in the same order of the respective elements in the list.
#' 'resolution' can be either 'daily' (or 'd') or 'subdaily' (or 's').
#' If \code{Data} is a vector, \code{Metadata} is ignored and all the
#' required information is read from the SEF files.
#' @param outpath Character string giving the path where the output is saved.
#' @param time_offset Numerical vector (of length 1 or equal to the number of 
#' analysed stations) of the number of hours to add to the time 
#' to obtain local time. This is used for tests on day and night temperature.
#' Recycled for all stations if only one value is given. Data not stored
#' in SEF files (i.e., not in UTC) are typically already expressed in local time:
#' in this case the offset is zero (the default).
#' 
#' @details
#' This is a wrapper of all functions that can be applied to the variables given
#' in \code{Data} (except the plotting functions).
#' 
#' \code{Data} can include any supported variable (see \link{Variables}) from 
#' different stations. The algorithm
#' will select the tests that can be applied to each variable. Note that some
#' tests require more than one variable from the same station.
#' 
#' This function produces flag files (one for each variable at each
#' station). The filenames follow the standard 
#' 'qc_<stationID>_<varcode>_<resolution>.txt'.
#' Each files contains a table of flagged values, with the last column
#' indicating the tests failed by each flagged observation.
#' 
#' The flag files can be edited by hand to remove or add flags. The flags can 
#' then be added to the 'Meta' column of SEF files by using the function 
#' \link{write_flags}.
#' 
#' @note
#' The tests will use their default parameters (e.g. thresholds). To use
#' custum parameters run the tests one by one.
#' 
#' @author Yuri Brugnara
#' 
#' @examples
#' # Testing all variables for Rosario de Santa Fe
#' 
#' # Create a data frame with all data from list Rosario
#' # For daily data we need to add the hour and minute columns (NAs)
#' Ros <- Rosario
#' Ros$Tx[,c("Hour","Minute")] <- NA
#' Ros$Tn[,c("Hour","Minute")] <- NA
#' Ros$rr[,c("Hour","Minute")] <- NA
#' Ros <- do.call("rbind", Ros)
#' Ros <- Ros[, c("Var","Year","Month","Day","Hour","Minute","Value")]
#' 
#' # Create a data frame with metadata including data resolution
#' df_meta <- do.call("rbind", Meta)
#' df_meta <- df_meta[which(df_meta$id=="Rosario"), ]
#' df_meta$res <- c("s", "s", "d", "d", "s", "s", "d", rep("s",4))
#' 
#' # Run all qc tests at once
#' # Time for Rosario is in UTC, therefore an offset is needed to get local time 
#' qc(Ros, df_meta, outpath = tempdir(), time_offset=-4.28)
#' 
#' 
#' # Testing one variable at one station
#' qc(Bern$ta, cbind(Meta$ta[which(Meta$ta$id=="Bern"),],"s"), 
#'    outpath = tempdir(), time_offset=0)
#'    
#' @export

qc <- function(Data, Metadata = NULL, outpath, time_offset = 0) {
  
  ## Check metadata and read SEF files
  if (class(Data) == "list") {
    if (is.null(dim(Metadata))) stop("Metadata must be a data frame or a matrix")
    if (dim(Metadata)[2] != 7) stop("Metadata must have 7 columns")
    if (length(Data) != length(unique(Metadata[,1]))) {
      stop(paste0("Number of unique station IDs in Metadata must correspond to the 
           number of elements of Data (", length(Data), ")"))
    }
    names(Data) <- unique(Metadata[,1])
  } else if (class(Data) %in% c("data.frame", "matrix")) {
    if (is.null(dim(Metadata))) stop("Metadata must be a data frame or a matrix")
    if (dim(Metadata)[2] != 7) stop("Metadata must have 7 columns")
    if (length(unique(Metadata[,1])) != 1) {
      stop("Number of unique station IDs in Metadata must correspond to the 
           number of elements of Data (1)")      
    }
    Data <- list(Data)
    names(Data) <- unique(Metadata[,1])
  } else {
    paths <- Data
    Data <- list()
    Metadata <- matrix(nrow=length(paths), ncol=7)
    for (i in 1:length(paths)) {
      id <- read_meta(paths[i], "id")
      v <- read_meta(paths[i], "var")
      if (id %in% names(Data)) {
        Data[[id]] <- rbind(Data[[id]], read_sef(paths[i], all=TRUE)[, c(1:6,8)])
      } else {
        Data[[id]] <- read_sef(paths[i], all=TRUE)[, c(1:6,8)]
      }
      stat <- read_meta(paths[i], "stat")
      if (v %in% c("Tx","Tn", "sd", "sc")) {
        datares <- "d"
      } else if (stat == "point") { 
        datares <- "s"
      } else if (sum(!is.na(Data[[id]][,5])) == 0) {
        datares <- "d"
      } else if (length(grep("h",Data[[id]][,7])) > 0 | 
                 length(unique(Data[[id]][,7])) > 1 |
                 unique(Data[[id]][,7])[1] == "p") {
        datares <- "s"
      } else {
        datares <- "d"
      }
      Metadata[i, ] <- c(read_meta(paths[i], c("id","lat","lon","alt","var","units")),
                         datares)
    }
  }
  Metadata <- as.data.frame(Metadata, stringsAsFactors = FALSE)
  names(Metadata) <- c("id", "lat", "lon", "alt", "var", "units", "res")

  
  ## Initialize data frames for flags and tests
  Flags <- data.frame(var = character(), year = numeric(), month = numeric(),
                      day = numeric(), hour = numeric(), minute = numeric(),
                      value = numeric(), test = character())
  
  ## Recycle offset
  if (length(time_offset) == 1) time_offset <- rep(time_offset, length(Data))
  if (length(time_offset) != length(Data)) stop("Length of time_offset must correspond
                                                 to the number of stations")
  
  ## Loop on stations
  for (st in names(Data)) {
    message("")
    message("----------------")
    message(st)
    message("----------------")
    names(Data[[st]]) <- c("var", "year", "month", "day", "hour", "minute", "value")
    Data$var <- as.character(Data$var)
    t_offset <- time_offset[which(names(Data)==st)]
    
    ## Look for possible tests for this station
    vars <- unique(Data[[st]]$var)
    st_tests <- Tests[which(Tests$Variable%in%c(vars,"any")), ]
#    st_tests <- st_tests[!duplicated(st_tests$Test), ]

    ## Loop on variables and tests
    if (nrow(st_tests) > 0) {
      for (vbl in vars) {
        
        ## Read metadata
        Meta <- subset(Metadata, Metadata$id==st & Metadata$var==vbl)
        
        ## Loop on resolutions
        for (res in Meta$res) {
          
          meta <- subset(Meta, Meta$res==res)
          if (nrow(meta) > 1) stop("Duplicate in Metadata (same variable and resolution)")
          
          ## Define columns to be used by the tests
          if (substr(res,1,1) == "s") {
            j <- 1:7
            r <- "subdaily"
          } else if (substr(res,1,1) == "d") {
            j <- c(1:4, 7)
            r <- "daily"
          } else {
            stop (paste("Unrecognized resolution for", vbl))
          }
          message("")
          message(paste0("[ ", vbl, " (", r, ") ]"))
          
          for (i_test in 1:nrow(st_tests)) {
            
            ## Check test requirements
            if (st_tests$Resolution[i_test] %in% c(r,"any") & 
                st_tests$Variable[i_test] %in% c(vbl,"any")) {
              #          if (is.na(st_tests$Requires2[i_test])) {
              if (st_tests$Requires[i_test] == "") {
                if (st_tests$Test[i_test] == "subdaily_out_of_range") {
                  get(st_tests$Test[i_test])(subset(Data[[st]][,j],Data[[st]][,1]==vbl),
                                             meta[1:6], outpath = outpath, 
                                             time_offset = t_offset)
                } else {
                  get(st_tests$Test[i_test])(subset(Data[[st]][,j],Data[[st]][,1]==vbl),
                                             meta[1:6], outpath = outpath)
                }
              } else if (st_tests$Requires[i_test] %in% vars) {
                meta2 <- rbind(meta, Metadata[which(Metadata$id==st & 
                                                      Metadata$var==st_tests$Requires[i_test] &
                                                      Metadata$res==res),])
                get(st_tests$Test[i_test])(subset(Data[[st]][,j], Data[[st]][,1]%in%
                                                    c(vbl,st_tests$Requires[i_test])),
                                           meta2[,1:6], outpath = outpath)
              }
              #          } else if (st_tests$Requires1[i_test] %in% vars & 
              #                     st_tests$Requires2[i_test] %in% vars) {
              #            meta2 <- rbind(meta, Metadata[which(Metadata$id==st & 
              #                                                 Metadata$var==st_tests$Requires1[i_test])])
              #            meta2 <- rbind(meta2, Metadata[which(Metadata$id==st & 
              #                                                 Metadata$var==st_tests$Requires2[i_test])])              
              #            get(st_tests$Test[i_test])(subset(Data[[st]][,j],
              #                                              Data[[st]][,1]%in%c(vbl,st_tests$Requires2[i_test])),
              #                                       meta2, outpath = outpath)
              #          }
            }
            
          }
          
        }
        
      }
    } else {
      warning(paste("No supported variables in station", st))
    }
    
  }
  
}