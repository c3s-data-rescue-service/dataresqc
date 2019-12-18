# Run in the package main directory to import example data and metadata,
# plus a table of available tests

library(devtools)

# Load and format Bern data from CHIMES project
load("data-raw/Bern_Studer.RData")
n <- dim(df)[1]
df <- df[order(df$year,df$month,df$day,df$hour,df$minutes),]
Bern <- list()
Bunits <- c("C","hPa")
names(Bunits) <- c("ta","p")
Bern$ta <- data.frame(Var = rep("ta", n),
                      Year = df$year,
                      Month = df$month,
                      Day = df$day,
                      Hour = df$hour,
                      Minutes = df$minutes,
                      Obs = df$air_temperature,
                      stringsAsFactors = FALSE)
Bern$ta <- Bern$ta[which(!is.na(Bern$ta$Obs)), ]
Bern$p <- data.frame(Var = rep("p", n),
                     Year = df$year,
                     Month = df$month,
                     Day = df$day,
                     Hour = df$hour,
                     Minutes = df$minutes,
                     Obs = df$air_pressure * 1013.25 / 760,
                     stringsAsFactors = FALSE)
Bern$p <- Bern$p[which(!is.na(Bern$p$Obs)), ]


# Load data from Rosario de Santa Fe (SEF files)
Rosario <- list()
Runits <- c()
for (variable in c("ta_maximum","ta_minimum","dd","n","p","rr","ta","tb","wind_force")) {
  Rosario[[variable]] <- read_sef(list.files("data-raw", 
                                             pattern=paste0("_",variable,".tsv"), 
                                             full.names = TRUE))
  if (variable == "ta_minimum") {
    Rosario[[variable]] <- cbind(variable, aggregate(Rosario[[variable]]$Value,
                                                     list(Rosario[[variable]]$Year,
                                                          Rosario[[variable]]$Month,
                                                          Rosario[[variable]]$Day),
                                                     FUN = min), 
                                 stringsAsFactors = FALSE)
    names(Rosario[[variable]]) <- c("Var", "Year", "Month", "Day", "Value")
    Rosario[[variable]] <- Rosario[[variable]][order(Rosario[[variable]]$Year,
                                                     Rosario[[variable]]$Month,
                                                     Rosario[[variable]]$Day), ]
    Rosario[[variable]] <- Rosario[[variable]][, 1:which(names(Rosario[[variable]])=="Value")]
    Rosario$Tn <- Rosario[[variable]]
    Rosario$Tn[, 1] <- "Tn"
    Rosario[[variable]] <- NULL
  }
  if (variable == "ta_maximum") {
    Rosario[[variable]] <- cbind(variable, aggregate(Rosario[[variable]]$Value,
                                                     list(Rosario[[variable]]$Year,
                                                          Rosario[[variable]]$Month,
                                                          Rosario[[variable]]$Day),
                                                     FUN = max), 
                                 stringsAsFactors = FALSE)
    names(Rosario[[variable]]) <- c("Var", "Year", "Month", "Day", "Value")
    Rosario[[variable]] <- Rosario[[variable]][order(Rosario[[variable]]$Year,
                                                     Rosario[[variable]]$Month,
                                                     Rosario[[variable]]$Day), ]
    Rosario[[variable]] <- Rosario[[variable]][, 1:which(names(Rosario[[variable]])=="Value")]
    Rosario$Tx <- Rosario[[variable]]
    Rosario$Tx[, 1] <- "Tx"
    Rosario[[variable]] <- NULL
  }
  if (variable == "rr") {
    Rosario[[variable]] <- cbind(variable, aggregate(Rosario[[variable]]$Value,
                                                     list(Rosario[[variable]]$Year,
                                                          Rosario[[variable]]$Month,
                                                          Rosario[[variable]]$Day),
                                                     FUN = sum), 
                                 stringsAsFactors = FALSE)
    names(Rosario[[variable]]) <- c("Var", "Year", "Month", "Day", "Value")
    Rosario[[variable]] <- Rosario[[variable]][order(Rosario[[variable]]$Year,
                                                     Rosario[[variable]]$Month,
                                                     Rosario[[variable]]$Day), ]    
  }
  if (variable == "wind_force") {
    Rosario[[variable]] <- subset(Rosario[[variable]], Rosario[[variable]]$Year %in% 1894:1900)
  }
  if (!variable %in% c("ta_minimum","ta_maximum")) {
    Rosario[[variable]] <- Rosario[[variable]][, 1:which(names(Rosario[[variable]])=="Value")]
  }
  Runits <- append(Runits, read_meta(list.files("data-raw", 
                                                pattern=paste0("_",variable,".tsv"), 
                                                full.names = TRUE), "units"))
}
# Convert to hPa
Rosario$p$Value <- Rosario$p$Value / 100
# Make ta and tb have same dimensions
j <- match(apply(Rosario$tb[,2:6], 1, paste, collapse="-"), 
           apply(Rosario$ta[,2:6], 1, paste, collapse="-"))
Rosario$ta <- Rosario$ta[j, ]
# Calculate relative humidity
Rosario$rh <- Rosario$tb
Rosario$rh$Var <- "rh"
ed <- 6.112 * exp(17.502*Rosario$ta$Value/(240.97+Rosario$ta$Value))
ew <- 6.112 * exp(17.502*Rosario$tb$Value/(240.97+Rosario$tb$Value))
Rosario$rh$Value <- 100 * (ew-0.6687451584*(1+0.00115*Rosario$tb$Value)*
                           (Rosario$ta$Value-Rosario$tb$Value)) / ed
Rosario$rh$Value <- round(Rosario$rh$Value, 0)
Runits <- append(Runits, "%")
# Calculate dew point (August-Roche-Magnus approximation)
Rosario$td <- Rosario$tb
Rosario$td$Var <- "td"
Rosario$td$Value <- 243.04 * (log(Rosario$rh$Value/100)+((17.625*Rosario$ta$Value)/
                                                     (243.04+Rosario$ta$Value))) /
  (17.625-log(Rosario$rh$Value/100)-((17.625*Rosario$ta$Value)/(243.04+Rosario$ta$Value)))
Rosario$td$Value <- round(Rosario$td$Value, 1)
Runits <- append(Runits, "C")
names(Runits) <- c("Tn","Tx","dd","n","p","rr","ta","tb","w","rh","td")
Runits["p"] <- "hPa"
# Transforme wind force to wind speed
Rosario$w <- Rosario$wind_force
Rosario$w$Var <- "w"
Rosario$wind_force <- NULL
beaufort <- c(0, 1, 2, 4, 7, 9, 12, 16, 19, 23, 26, 31, 35)
Rosario$w$Value <- beaufort[Rosario$w$Value+1]
Runits["w"] <- "m/s"

# Create metadata table
Meta <- list()
ids <- c("Bern", "Rosario")
lats <- c(46.94812, -32.945)
lons <- c(7.45196, -60.333)
alts <- c(534, 36)
for (variable in unique(c(names(Bern),names(Rosario)))) {
  if (variable %in% names(Bern) & variable %in% names(Rosario)) {
    Meta[[variable]] <- data.frame(id = ids, lat = lats, lon = lons, alt = alts,
                                   var = rep(variable, 2),
                                   units = c(Bunits[variable],Runits[variable]),
                                   stringsAsFactors = FALSE)
  } else if (variable %in% names(Bern)) {
    Meta[[variable]] <- data.frame(id = ids[1], lat = lats[1], lon = lons[1], alt = alts[1],
                                   var = variable, units = Bunits[variable],
                                   stringsAsFactors = FALSE)    
  } else {
    Meta[[variable]] <- data.frame(id = ids[2], lat = lats[2], lon = lons[2], alt = alts[2],
                                   var = variable, units = Runits[variable],
                                   stringsAsFactors = FALSE)     
  }
}

# Create tests table
Tests <- read.csv("data-raw/tests.csv", stringsAsFactors = FALSE)
names(Tests)[1] <- "Test"

# Create variables table
Variables <- read.table("data-raw/variables.txt", sep = "\t", stringsAsFactors = FALSE)
names(Variables) <- c("abbr", "full_name")

# Import using devtools
use_data(Bern, Rosario, Meta, Tests, Variables, overwrite = TRUE)