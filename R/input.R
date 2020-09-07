#' Read data files in Station Exchange Format version 1.0.0
#'
#' @param file Character string giving the path of the SEF file.
#' @param all If FALSE (the default), omit the columns 'Period' and 'Meta'
#' (also 'Hour' and 'Minute' for non-instantaneous data)
#'
#' @return A data frame with up to 9 variables, depending on whether
#' \code{all} is set to TRUE.
#' The variables are: variable code, year, month, day, hour, minute,
#' value, period, metadata.
#'
#' @author Yuri Brugnara
#'
#' @import utils
#' @export

read_sef <- function(file = file.choose(), all = FALSE) {
  
  ## Read from the header
  meta <- read_meta(file)
  varcode <- meta["var"]
  timeres <- meta["stat"]
  
  ## Read the data
  Data <- read.table(file = file, skip = 12, header = TRUE, fill = TRUE,
                     sep = "\t", stringsAsFactors = FALSE, quote = "")
  
  ## Select columns
  if (!all) {
    if (timeres == "point" & !varcode %in% c("sd","sc")) {
      Data <- Data[, c(1:5,7)]
      colnames(Data) <- c("Year", "Month", "Day", "Hour", "Minute", "Value")
    } else {
      Data <- Data[, c(1:3,7)]
      colnames(Data) <- c("Year", "Month", "Day", "Value")
    }
  } else {
    Data$Meta[which(is.na(Data$Meta))] <- ""
  }
  
  Data <- cbind(rep(varcode, dim(Data)[1]), Data, stringsAsFactors = FALSE)
  colnames(Data)[1] <- "Var"
  
  return(Data)
  
}


#' Read metadata from the Station Exchange Format version 1.0.0
#'
#' @param file Character string giving the path of the data file.
#' @param parameter Character vector of required parameters. Accepted
#' values are \code{"version"}, \code{"id"}, \code{"name"}, \code{"lat"},
#' \code{"lon"}, \code{"alt"}, \code{"source"}, \code{"link"},
#' \code{"var"}, \code{"stat"}, \code{"units"}, \code{"meta"}.
#' By default all parameters are read at once.
#'
#' @return A character vector with the required parameters.
#'
#' @author Yuri Brugnara
#'
#' @import utils
#' @export

read_meta <- function(file = file.choose(), parameter = NULL) {
  
  ## Read header
  header <- read.table(file = file, quote = "", comment.char = "", sep = "\t",
                       nrows = 12, stringsAsFactors = FALSE, fill = TRUE)
  pars <-c("version", "id", "name", "lat", "lon", "alt", "source", "link",
           "var", "stat", "units", "meta")
  
  ## Check format
  if (header[1,1] != "SEF") stop("This is not a SEF file")
  if (!header[1,2] %in% c("0.2.0","1.0.0")) {
    stop(paste("This function is not compatible with SEF version", header[1,2]))
  }
  
  ## Extract metadata
  if (is.null(parameter)) {
    out <- header[, 2]
    names(out) <- pars
  } else {
    out <- header[match(parameter, pars), 2]
    names(out) <- parameter
  }
  
  return(out)
  
}


#' Check compliance with SEF guidelines
#'
#' @param file Character string giving the path of the SEF file.
#'
#' @return TRUE if no errors are found, FALSE otherwise.
#'
#' @author Yuri Brugnara
#'
#' @note
#' For more information on error/warning messages produced by this
#' function see the SEF documentation.
#'
#' @import utils
#' @export

check_sef <- function(file = file.choose()) {
  
  meta <- read_meta(file)
  header <- read.table(file = file, nrows = 12, fill = TRUE,
                       sep = "\t", stringsAsFactors = FALSE)
  x <- read.table(file = file, skip = 12, header = TRUE, fill = TRUE,
                  sep = "\t", stringsAsFactors = FALSE)
  e <- 0
  w <- 0
  
  
  ## Check header
  if (!all(header[, 1] == c("SEF", "ID", "Name", "Lat", "Lon", "Alt",
                            "Source","Link", "Vbl", "Stat", "Units", "Meta"))) {
    message("ERROR: One or more labels in the header were not recognized, 
            or they are in the wrong order")
    e <- e + 1
  }
  
  if (is.na(meta["id"]) | meta["id"] == "") {
    message("ERROR: Missing Station ID")
    e <- e + 1
  } else {
    if (any(is.na(utf8ToInt(meta["id"])))) {
      message("ERROR: Special characters are not allowed in the Station ID")
      e <- e + 1
    }
    if (grepl(" ", trimws(meta["id"]))) {
      message("ERROR: Blanks are not allowed in the Station ID")
      e <- e + 1
    }
  }
    
  if (is.na(meta["name"])) {
    message("ERROR: Missing Station Name")
    e <- e + 1
  }
  
  if (is.na(meta["lat"])) {
    message("ERROR: Missing Latitude")
    e <- e + 1
  }
  
  if (is.na(meta["lon"])) {
    message("ERROR: Missing Longitude")
    e <- e + 1
  }
  
  if (is.na(meta["alt"])) {
    message("Warning: Missing Altitude")
    w <- w + 1
  }
  
  if (is.na(meta["source"])) {
    message("Warning: Missing Data Source")
    w <- w + 1
  }
  
  if (!meta["var"] %in% Variables$abbr) {
    message("Warning: Non-standard variable code")
    w <- w + 1
  }
  
  if (is.na(meta["stat"])) {
    message("ERROR: Missing statistic")
    e <- e + 1
  }
  
  if (!meta["stat"] %in% c("point","maximum","minimum","mean","median",
                           "mid_range","mode","sum","variance","standard_deviation")) {
    message("Warning: stat entry not recognized")
    w <- w + 1
  }
  
  if (is.na(meta["units"])) {
    message("ERROR: Missing units")
    e <- e + 1
  }
  
  if (!is.na(meta["meta"]) & meta["meta"] != "") {
    n1 <- length(strsplit(meta["meta"], "|", fixed = TRUE)[[1]])
    n2 <- length(strsplit(meta["meta"], "=", fixed = TRUE)[[1]])
    if (n2 != (n1+1)) {
      message("Warning: Uncorrect metadata format in Meta in the header")
      w <- w + 1
    }
  }
  
  
  ## Check data
  if (!all(names(x) == c("Year", "Month", "Day", "Hour", "Minute", "Period",
                         "Value", "Meta"))) {
    message("ERROR: One or more data headers were not recognized")
    e <- e + 1
  }
  
  for (lab in c("Year", "Month", "Day", "Hour", "Minute")) {
    x[[lab]] <- suppressWarnings(as.numeric(x[[lab]]))
    i <- which(as.integer(x[[lab]]) != x[[lab]])
    if (length(i) > 0) {
      message(paste("ERROR: Non-integer", lab))
      e <- e + 1
    }
  }
  
  if (sum(is.na(x$Year)) > 0) {
    message("Warning: There are missing or non-numeric values in the Year column")
    w <- w + 1
  }
  
  current_year <- as.integer(substr(date(), nchar(date())-4, nchar(date())))
  i <- which(x$Year < 1600 | x$Year > current_year)
  if (length(i) > 0) {
    message ("Warning: Year outside 1600-present range")
    w <- w + 1
  }
  
  i <- which(x$Month < 1 | x$Month > 12)
  if (length(i) > 0) {
    message ("ERROR: Month outside 1-12 range")
    e <- e + 1
  }
  
  n <- sum(!is.na(x$Month) & grepl("year", x$Period, ignore.case = TRUE))
  if (n > 0) {
    message("Warning: If the annual values refer to the calendar year the Month column should contain missing values")
    w <- w + 1
  }
  
  i <- which(x$Day < 1 | x$Day > 31)
  if (length(i) > 0) {
    message ("ERROR: Day outside 1-31 range")
    e <- e + 1
  }
  
  n <- sum(!is.na(x$Day) & grepl("month", x$Period, ignore.case = TRUE))
  if (n > 0) {
    message("Warning: If the monthly values refer to the calendar month the Day column should contain missing values")
    w <- w + 1
  }
  
  j <- grep("year", x$Period, ignore.case = TRUE)
  if (sum(is.na(x$Month)) > 0) {
    i <- which(is.na(x$Month))
    if(!all(i %in% j)) {
      message("Warning: There are missing or non-numeric values in the Month column where data are not annual")
      w <- w + 1
    }
  }
  
  i <- which(x$Hour < 0 | x$Hour > 24)
  if (length(i) > 0) {
    message("ERROR: Hour outside 0-24 range")
    e <- e + 1
  }
  
  i <- which(x$Minute < 0 | x$Minute > 59)
  if (length(i) > 0) {
    message("ERROR: Minute outside 0-59 range")
    e <- e + 1
  }
  
  i <- which(x$Hour == 0 & x$Minute == 0)
  if (length(i) > 0 & meta["stat"] != "point") {
    message("Warning: 24 is recommended instead of 0 for the Hour column")
    w <- w + 1
  }
  
  k <- grep("month", x$Period, ignore.case = TRUE)
  for (lab in c("Day", "Hour", "Minute")) {
    if (sum(is.na(x[[lab]])) > 0) {
      i <- which(is.na(x[[lab]]))
      if(!all(i %in% c(j, k))) {
        message(paste("Warning: There are missing or non-numeric values in the",
                      lab, "column but data are neither monthly nor annual"))
        w <- w + 1
      }
    }
  }
  
  n <- sum(is.na(x$Period))
  if (n > 0) {
    message("ERROR: There are missing values in the Period column")
    e <- e + 1
  }
  
  n <- length(grep(" ", trimws(x$Period)))
  if (n > 0) {
    message("ERROR: Blanks are not allowed in the Period column")
    e <- e + 1
  }
  
  if (length(c(which(trimws(x$Period) == "0"),
               which(trimws(x$Period) == "p"),
               grep("second", x$Period, ignore.case = TRUE),
               grep("minute", x$Period, ignore.case = TRUE),
               grep("hour", x$Period, ignore.case = TRUE),
               grep("day", x$Period, ignore.case = TRUE),
               grep("month", x$Period, ignore.case = TRUE),
               grep("year", x$Period, ignore.case = TRUE))) != (nrow(x) - n)) {
    message("Warning: There are unrecognized values in the Period column")
    e <- e + 1
  }
  
  if (length(grep("p", x$Period)) > 0) {
    chronol <- order(x$Year, x$Month, x$Day)
    if (!all(chronol == 1:nrow(x))) {
      message("Warning: Data should be in cronological order when using p in the Period column")
      w <- w + 1
    }
  }
  
  n <- sum(is.na(x$Value))
  x$Value <- suppressWarnings(as.numeric(x$Value))
  if (sum(is.na(x$Value)) > n) {
    message("Warning: There are non-numeric values in the Value column")
    w <- w + 1
  }
  
  if (sum(!is.na(x$Meta)) > 0) {
    n1 <- sapply(x$Meta, function(x) length(strsplit(x, "|", fixed = TRUE)[[1]]))
    n2 <- sapply(x$Meta, function(x) length(strsplit(x, "=", fixed = TRUE)[[1]]))
    if (any(n2 != (n1+1))) {
      message("Warning: Uncorrect format detected in the Meta column")
      w <- w + 1
    }
  }
  
  if (meta["var"] == "p" & !grepl("PGC", meta["meta"])) {
    message("Warning: PGC not specified")
    w <- w + 1
  }
  
  if (meta["var"] == "p" & !grepl("PTC", meta["meta"])) {
    message("Warning: PTC not specified")
    w <- w + 1
  }
  
  
  ## Check consistency between stat and period
  if (meta["stat"] == "point" & !all(trimws(x$Period) == "0")) {
    i <- which(trimws(x$Period) != "0")
    j <- grep("stat", x$Meta, ignore.case = TRUE)
    if (!all(i %in% j)) {
      message("ERROR: Period must be 0 when stat is point")
      e <- e + 1
    }
  }
  
  if (meta["stat"] != "point" & any(trimws(x$Period) == "0")) {
    i <- which(trimws(x$Period) == "0")
    j <- grep("stat", x$Meta, ignore.case = TRUE)
    if (!all(i %in% j)) {
      message("ERROR: Period cannot be 0 when stat is not point")
      e <- e + 1
    }
  }
  
  
  ## Print recap and return a boolean
  message(paste(e, "errors and", w, "warnings"))
  if (e == 0) {
    return(TRUE)
  } else {
    return(FALSE)
  }
  
}
