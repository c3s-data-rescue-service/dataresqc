#' Read data files in Station Exchange Format version 0.2.0
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
  v <- meta["version"]
  
  ## Read the data
  Data <- read.table(file = file, skip = 12, header = TRUE, fill = TRUE,
                     sep = "\t", stringsAsFactors = FALSE)
  
  ## Select columns
  if (!all) {
    if (timeres == "point") {
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


#' Read metadata from the Station Exchange Format version 0.2.0
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
  if (header[1,2] != "0.2.0") stop(paste("This function is not compatible with
                                         SEF version", header[1,2]))
  
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