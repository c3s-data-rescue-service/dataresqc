#' Duplicate columns test
#' 
#' Looks for data that have been digitized twice by mistake. For sud-daily data,
#' this is done by looking for series of zero differences between adjacent
#' observation times. For daily data, by looking for series of zero differences 
#' between the same days of adjacent months.
#' 
#' @param Data A character string giving the path of the input file,
#' or a matrix with 5 (7) columns for daily (sub-daily) data: variable code, year, 
#' month, day, (hour), (minute), value.
#' @param meta A character vector with 6 elements: station ID, latitude, longitude,
#' altitude, variable code, units. If \code{Data} is a path, \code{meta} is
#' ignored.
#' @param outpath Character string giving the path for the QC results.
#' @param ndays Number of consecutive days with zero difference required to
#' flag the data. The default is 5.
#' 
#' @details The input file must follow the Copernicus Station Exchange 
#' Format (SEF). This function works with any numerical variable.
#' 
#' Zeroes are automatically excluded in bounded variables such as precipitation.
#' 
#' @author Yuri Brugnara
#'
#' @examples 
#' climatic_outliers(Rosario$Tn, Meta$Tn, outpath = tempdir(), ndays = 3) 
#'
#' @export
 
duplicate_columns <- function(Data, meta = NULL, outpath, ndays = 5) {
  
  #Read data and metadata
  if (is.null(dim(Data))) {
    meta <- read_meta(Data, c("id","lat","lon","alt","var","units"))
    Data <- read_sef(Data)
  }
  if (length(meta) != 6) stop("Incorrect metadata (must be a vector of 6 elements)")
  n <- dim(Data)[2]
  if (!n %in% c(5,7)) stop("Incorrect dimension of Data")
  meta <- as.character(meta)
  
  #Define flag
  flag <- "duplicate_columns"
  
  #Check units
  Data[,n] <- check_units(Data[,n], meta[5], meta[6])
  
  #Remove zeroes and missing values
  if (meta[5] %in% c("rr","sd","fs","sc","sw")) {
    Data <- Data[which(Data[,n] != 0), ]
  }
  Data <- Data[which(!is.na(Data[,n])), ]
  
  #Flag daily data
  if (n == 5) {
    out <- data.frame(Var=character(), Year=numeric(), Month=numeric(),
                      Day=numeric(), Value=numeric())
    Data <- Data[order(Data[,2], Data[,3], Data[,4]), ]
    Data$ym <- paste(Data[,2], Data[,3])
    ym <- unique(Data$ym)
    for (i in 2:length(ym)) {
      difference <- suppressWarnings(Data[which(Data$ym == ym[i-1]), 5] - 
                                       Data[which(Data$ym == ym[i]), 5])
      counts <- rle(difference)
      j <- which(counts$values == 0 & counts$lengths >= ndays)
      if (length(j) > 0) {
        for (k in j) {
          dupl_subset <- Data[which(Data$ym == ym[i]), ]
          dupl_subset <- dupl_subset[(max(cumsum(counts$lengths[1:k])) - counts$lengths[k] + 1):
                                       max(cumsum(counts$lengths[1:k])), 1:5]
          out <- rbind(out, dupl_subset)
        }
      }
    }
    
  #Flag subdaily data
  } else {
    out <- data.frame(Var=character(), Year=numeric(), Month=numeric(),
                      Day=numeric(), Hour=numeric(), Minute=numeric(),
                      Value=numeric())
    Data <- Data[order(Data[,2], Data[,3], Data[,4], Data[,5], Data[,6]), ]
    Data$ymd <- paste(Data[,2], Data[,3], Data[,4])
    nobs <- aggregate(Data[,7], list(Data[,2], Data[,3], Data[,4]), length)
    nobs <- nobs$x[order(nobs[,1], nobs[,2], nobs[,3])]
    Data$obs_number <- unlist(lapply(nobs, function(x) 1:x))
    nobs <- max(Data$obs_number)
    if (nobs > 1) {
      Data_by_obstime <- list()
      for (i in 1:nobs) {
        Data_by_obstime[[i]] <- aggregate(Data[,7], list(Data[,2], Data[,3], Data[,4]), 
                                          function(x) x[i])
        Data_by_obstime[[i]] <- Data_by_obstime[[i]][order(Data_by_obstime[[i]][,1],
                                                           Data_by_obstime[[i]][,2],
                                                           Data_by_obstime[[i]][,3]), ]
        Data_by_obstime[[i]]$ymd <- paste(Data_by_obstime[[i]][,1],
                                          Data_by_obstime[[i]][,2],
                                          Data_by_obstime[[i]][,3])
      }
      for (i in 2:nobs) {
        difference <- Data_by_obstime[[i]]$x - Data_by_obstime[[i-1]]$x
        counts <- rle(difference)
        j <- which(counts$values == 0 & counts$lengths >= ndays)
        if (length(j) > 0) {
          for (k in j) {
            dupl_subset <- 
              Data[which(Data$ymd %in% Data_by_obstime[[i]]$ymd[(max(cumsum(counts$lengths[1:k])) - 
                                                             counts$lengths[k] + 1):
                                                            max(cumsum(counts$lengths[1:k]))] &
                           Data$obs_number %in% c(i-1,i)), 1:7]
            out <- rbind(out, dupl_subset)
          }
        }
      }
    }
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
  outMsg <- "Duplicate columns test completed"
    
  message(outMsg)
  
}