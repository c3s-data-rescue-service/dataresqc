#' Write data in Station Exchange Format version 1.0.0
#'
#' @param Data A data frame with 6 variables in this order:
#' year, month, day, hour, minute, value.
#' @param outpath Character string giving the output path (note that the
#' filename is generated from the source identifier, station code, start
#' and end dates, and variable code).
#' @param variable Variable code. This is a required field.
#' @param cod Station code. This is a required field.
#' @param nam Station name.
#' @param lat Station latitude (degrees North in decimal).
#' @param lon Station longitude (degreees East in decimal).
#' @param alt Station altitude (metres).
#' @param sou Character string giving the source identifier.
#' @param link Character string giving an url for metadata (e.g., link to the
#' C3S Data Rescue registry).
#' @param stat Character string giving the statistic code. This is a required
#' field.
#' @param units Character string giving the units. This is a required field.
#' @param metaHead Character string giving metadata entries for the header
#' (pipe separated).
#' @param meta Character vector with length equal to the number of rows
#' of \code{Data}, giving metadata entries for the single observations (pipe
#' separated).
#' @param period Observation time period code. Must be a character vector with
#' length equal to the number of rows of \code{Data} unless all observations
#' have the same period code.
#' @param time_offset Numerical vector of offsets from UTC in hours.
#' This value will be subtracted from the observation times to obtain UTC times,
#' so for instance the offset of Central European Time is +1 hour.
#' Recycled for all observations if only one value is given.
#' @param note Character string to be added to the end of the standard output
#' filename. It will be separated from the rest of the name by an underscore.
#' Blanks will be also replaced by underscores.
#' @param keep_na If FALSE (the default), lines where observations are NA are
#' removed.
#' @param outfile Output filename. If specified, ignores \code{note}.
#'
#' @note
#' Times in SEF files must be expressed in UTC.
#'
#' If \code{outfile} is not specified, the output filename is generated
#' automatically as \code{sou}_\code{cod}_startdate_enddate_\code{variable}.tsv
#'
#' @author Yuri Brugnara
#'
#' @examples
#' # Create a basic SEF file for air temperature in Bern
#' # (assuming the observation times are in local solar time)
#' # The file will be written in the working directory
#' meta_bern <- Meta$ta[which(Meta$ta$id == "Bern"), ]
#' write_sef(Bern$ta[, 2:7], outpath = tempdir(), variable = "ta", 
#'           cod = meta_bern$id, nam = "Bern", lat = meta_bern$lat, 
#'           lon = meta_bern$lon, alt = meta_bern$alt, 
#'           units = meta_bern$units, stat = "point", period = "0", 
#'           time_offset = meta_bern$lon * 24 / 360)
#'
#' @import utils
#' @export

write_sef <- function(Data, outpath, variable, cod, nam = "", lat = "",
                      lon = "", alt = "", sou = "", link = "", units,
                      stat, metaHead = "", meta = "", period = "",
                      time_offset = 0, note = "", keep_na = FALSE, outfile = NA) {
  
  ## Build filename
  if (substr(outpath, nchar(outpath), nchar(outpath)) != "/") {
    outpath <- paste0(outpath, "/")
  }
  if (is.na(outfile)) {
    datemin <- paste(formatC(unlist(Data[1, 1:3]), width=2, flag=0),
                     collapse = "")
    datemax <- paste(formatC(unlist(Data[dim(Data)[1], 1:3]), width=2, flag=0),
                     collapse = "")
    if (substr(datemin,5,6) == "NA") datemin <- substr(datemin, 1, 4)
    if (substr(datemax,5,6) == "NA") datemax <- substr(datemax, 1, 4)
    if (substr(datemin,3,4) == "NA") datemin <- substr(datemin, 1, 2)
    if (substr(datemax,3,4) == "NA") datemax <- substr(datemax, 1, 2)
    dates <- paste(datemin, datemax, sep = "-")
    filename <- paste(sou, cod, dates, variable, sep = "_")
    if (sou %in% c(NA,"")) filename <- sub("_", "", filename)
    if (note != "") {
      note <- paste0("_", gsub(" ", "_", note))
    }
    filename <- paste0(outpath, filename, note, ".tsv")
  } else {
    filename <- paste0(outpath, outfile)
    if (substr(filename, nchar(filename)-3, nchar(filename)) != ".tsv") {
      filename <- paste0(filename, ".tsv")
    }
  }
  
  ## Build header
  header <- array(dim = c(12, 2), data = "")
  header[1, ] <- c("SEF", "1.0.0")
  header[2, ] <- c("ID", trimws(as.character(cod)))
  header[3, ] <- c("Name", trimws(as.character(nam)))
  header[4, ] <- c("Lat", trimws(as.character(lat)))
  header[5, ] <- c("Lon", trimws(as.character(lon)))
  header[6, ] <- c("Alt", trimws(as.character(alt)))
  header[7, ] <- c("Source", trimws(as.character(sou)))
  header[8, ] <- c("Link", trimws(as.character(link)))
  header[9, ] <- c("Vbl", trimws(as.character(variable)))
  header[10, ] <- c("Stat", trimws(as.character(stat)))
  header[11, ] <- c("Units", trimws(as.character(units)))
  header[12, ] <- c("Meta", trimws(as.character(metaHead)))
  
  ## For instantaneous observations the period must be 0
  if (stat == "point" & !all(as.character(period) == "0")) {
    period <- "0"
    warning("Period forced to 0 because of 'stat'")
  }
  
  ## Convert times to UTC
  if (!all(time_offset == 0) & !all(is.na(as.integer(Data[,4])+as.integer(Data[,5])))) {
    times <- ISOdate(Data[,1], Data[,2], Data[,3], Data[,4], Data[,5])
    times <- times - time_offset * 3600
    Data[which(!is.na(times)), 1] <- as.integer(substr(times[which(!is.na(times))],1,4))
    Data[which(!is.na(times)), 2] <- as.integer(substr(times[which(!is.na(times))],6,7))
    Data[which(!is.na(times)), 3] <- as.integer(substr(times[which(!is.na(times))],9,10))
    Data[which(!is.na(times)), 4] <- as.integer(substr(times[which(!is.na(times))],12,13))
    Data[which(!is.na(times)), 5] <- as.integer(substr(times[which(!is.na(times))],15,16))
  }
  
  ## Build data frame with SEF structure
  DataNew <- data.frame(Year = as.integer(Data[, 1]),
                        Month = as.integer(Data[, 2]),
                        Day = as.integer(Data[, 3]),
                        Hour = Data[, 4],
                        Minute = Data[, 5],
                        Period = as.character(period),
                        Value = as.character(Data[, 6]),
                        Meta = as.character(meta),
                        stringsAsFactors = FALSE)
  
  ## Remove lines with missing data
  if (!keep_na) DataNew <- DataNew[which(!is.na(DataNew$Value)), ]
  
  ## Write header to file
  write.table(header, file = filename, quote = FALSE, row.names = FALSE,
              col.names = FALSE, sep = "\t", dec = ".", fileEncoding = "UTF-8")
  
  ## Write column names to file
  write.table(t(names(DataNew)), file = filename, quote = FALSE, row.names = FALSE,
              col.names = FALSE, sep = "\t", fileEncoding = "UTF-8",
              append = TRUE)
  
  ## Write data to file
  write.table(DataNew, file = filename, quote = FALSE, row.names = FALSE,
              col.names = FALSE, sep = "\t", dec = ".", fileEncoding = "UTF-8",
              append = TRUE)
  
  message(paste("Data written to file", filename))
  
}


###############################################################################


#' Add quality flags to a data file in Station Exchange Format version 1.0.0
#'
#' @param infile Character string giving the path of the SEF file.
#' @param qcfile Character string giving the path of the file with
#' the quality flags as produced with the QC tests. 
#' This file must have 6 (8) tab-separated columns
#' for daily (sub-daily) data: variable code, year, month, day, (hour),
#' (minute), value, semicolon(';')-separated failed tests. 
#' @param outpath Character string giving the output path.
#' @param note Character string to be added to the end of the name of the
#' input file to form the output filename. 
#' It will be separated from the rest of the name by an underscore.
#' Blanks will be also replaced by underscores.
#' If not specified, input and output filenames will be identical.
#' 
#' @author Yuri Brugnara
#' 
#' @note
#' The data will be converted to the standard units adopted by the qc.
#' An exception is made for cloud cover (oktas will not be converted).
#' 
#' @import utils
#' @export

write_flags <- function(infile, qcfile, outpath, note = "") {
  
  ## Read SEF file
  Data <- read_sef(infile, all = TRUE)
  header <- read.table(file = infile, quote = "", comment.char = "", sep = "\t",
                       nrows = 12, stringsAsFactors = FALSE, fill = TRUE)
  header[which(is.na(header[,2])), 2] <- ""
  
  ## Check units
  vbl <- read_meta(infile,"var")
  uts <- read_meta(infile,"units")
  Data$Value <- check_units(Data$Value, vbl, uts)
  if (vbl %in% c("ta","tb","td","t_air","t_wet","t_dew","Tx","Tn","dep_dew","ibt",
               "atb","Txs","TGs","Tns","TGn","t_snow","Ts","t_water")) {
    uts <- "C"
  } else if (vbl %in% c("p","mslp","pppp")) {
    uts <- "hPa"
  } else if (vbl %in% c("rr","sw","rrls")) {
    uts <- "mm"
  } else if (vbl %in% c("sd","fs")) {
    uts <- "cm"
  } else if (vbl == "w") {
    uts <- "m/s"
  } else if (vbl == "rh") {
    uts <- "%"
  }
  
  ## Read flags
  flags <- read.table(qcfile, stringsAsFactors = FALSE, header = TRUE, sep = "\t")
  if (flags[1,1] != Data[1,1]) {
    stop(paste("Variable mismatch:", flags[1,1], Data[1,1]))
  }
  if (dim(flags)[2] == 6) { # daily resolution
    dates <- paste(flags[,2], flags[,3], flags[,4], flags[,5])
    i <- which(paste(Data$Year, Data$Month, Data$Day, Data$Value) %in% dates)
  } else if (dim(flags)[2] == 8) { # subdaily resolution
    dates <- paste(flags[,2], flags[,3], flags[,4], flags[,5], flags[,6], flags[,7])
    i <- which(paste(Data$Year, Data$Month, Data$Day, Data$Hour, Data$Minute, Data$Value) 
               %in% dates)   
  } else {
    stop("Wrong number of columns in flags file")
  }
  
  if (length(i) > 0) {
    
    ## Stop if there are duplicates
    if (length(i) > nrow(flags)) stop("Duplicates need one flag for each instance")
    
    ## Check if all flagged times are present in the SEF file
    if (length(i) < nrow(flags)) {
      k <- which(!dates %in% paste(Data$Year, Data$Month, Data$Day, Data$Hour, Data$Minute,
                                   Data$Value))
      warning(paste("The SEF file does not contain all flagged observations.",
                    "Flags for the following observations could not be written:\n",
                    paste(dates[k], collapse="\n")))
      flags <- flags[-k, ]
    }
    
    ## Write flags
    Data$Meta[i] <- paste0(Data$Meta[i], "|qc=", flags[, dim(flags)[2]])
    Data$Meta[which(substr(Data$Meta, 1, 1) == "|")] <- 
      sub("[|]", "", Data$Meta[which(substr(Data$Meta, 1, 1) == "|")])
    
    ## Add name of QC package to the header
    meta_string <- paste0("QC software=dataresqc v", packageVersion("dataresqc"))
    if (header[12, 2] == "") {
      header[12, 2] <- meta_string
    } else {
      header[12, 2] <- paste(header[12, 2], meta_string, sep = "|")
    }
    
    ## Write SEF file with flags
    filename <- strsplit(infile, "/")[[1]]
    filename <- filename[length(filename)]
    if (note != "") filename <- paste0(substr(filename, 1, nchar(filename)-4),
                                       "_", gsub(" ", "_" , note), ".tsv")
    write_sef(Data = Data[, c(2:6,8)],
              outpath = outpath,
              variable = vbl,
              cod = header[2, 2],
              nam = header[3, 2],
              lat = header[4, 2],
              lon = header[5, 2],
              alt = header[6, 2],
              sou = header[7, 2], 
              link = header[8, 2],
              stat = header[10, 2],
              units = uts, 
              metaHead = header[12, 2], 
              meta = Data[, 9], 
              period = Data[, 7], 
              outfile = filename)
    
  } else warning("No matches found: possibly incorrect input files")
  
}


###############################################################################


## Write daily quality flags in the intermediate format

write_intermediate_daily <- function(namefile, out) {
  
  cols <- c("Var", "Year", "Month", "Day", "Value")
  names(out) <- c(cols, "Test")
  out <- out[order(out$Year,out$Month,out$Day), ]
  if (file.exists(namefile)) {
    qcfile <- read.table(namefile, header = TRUE)
    names(qcfile) <- c(cols, "Test")
    new_data <- merge(qcfile, out, by = cols, all = TRUE)
    
    ## Remove duplicates
    new_data$Test.y[grep(out$Test[1], new_data$Test.x)] <- NA
    
    ## Paste flags
    new_data$Test <- paste(new_data$Test.x, new_data$Test.y, sep = ";")
    
    ## Remove NAs from Test column
    new_data$Test[which(is.na(new_data$Test.x))] <- 
      as.character(new_data$Test.y[which(is.na(new_data$Test.x))])
    new_data$Test[which(is.na(new_data$Test.y))] <- 
      as.character(new_data$Test.x[which(is.na(new_data$Test.y))])
    
    ## Write to file
    new_data <- as.data.frame(new_data[, c(cols, "Test")])
    new_data <- new_data[order(new_data$Year,new_data$Month,new_data$Day), ]
    write.table(new_data, file = namefile, quote = FALSE, sep = "\t", row.names = FALSE)
  } else {
    write.table(out, file = namefile, quote = FALSE, sep = "\t", row.names = FALSE)
  }
  
}


###############################################################################


## Write sub-daily quality flags in the intermediate format

write_intermediate_subdaily <- function(namefile, out) {
  
  cols <- c("Var", "Year", "Month", "Day", "Hour", "Minute", "Value")
  names(out) <- c(cols, "Test")
  out <- out[order(out$Year,out$Month,out$Day,out$Hour,out$Minute), ]
  if (file.exists(namefile)) {
    qcfile <- read.table(namefile, header = TRUE)
    names(qcfile) <- c(cols, "Test")
    new_data <- merge(qcfile, out, by = cols, all = TRUE)
    
    ## Remove duplicates
    new_data$Test.y[grep(out$Test[1], new_data$Test.x)] <- NA
    
    ## Paste flags
    new_data$Test <- paste(new_data$Test.x, new_data$Test.y, sep = ";")
    
    ## Remove NAs from Test column
    new_data$Test[which(is.na(new_data$Test.x))] <- 
      as.character(new_data$Test.y[which(is.na(new_data$Test.x))])
    new_data$Test[which(is.na(new_data$Test.y))] <- 
      as.character(new_data$Test.x[which(is.na(new_data$Test.y))])
    
    ## Write to file
    new_data <- as.data.frame(new_data[, c(cols, "Test")])
    new_data <- new_data[order(new_data$Year,new_data$Month,new_data$Day,
                               new_data$Hour,new_data$Minute), ]
    write.table(new_data, file = namefile, quote = FALSE, sep = "\t", row.names = FALSE)
  } else {
    write.table(out, file = namefile, quote = FALSE, sep = "\t", row.names = FALSE)
  }
  
}
