# US ICE Wholesale Electricity Data processing step 2
# 1. use 5 hubs started from 2001
# 2. select variables
# 3. tidy up
# 4. inflation adjustment (base year = 2022) using US CPI

# @author: Konan Hara

# @date: November 18, 2024

# @software: R version 4.4.1

# clear environment
rm(list = ls())

# package loading
library(here)
library(data.table)

# path to data
pathtodata <- readLines(file(here("data.txt")), n = 1, warn = FALSE)

# load US CPI data
cpi_data <- fread(file.path(pathtodata, "processed", "cpi_2001_2022.csv"))

# load and clean US ICE Wholesale Electricity Data

# format from 2001--2013
# list of files for selected five hubs
ice_electric_to2013 <- list()
# Mass, PJM West, Mid-C, Palo Verde, and SP-15
hub_filelist_to2013 <- c(
	  "Mass"
	, "PJMWest"
	, "MidC"
	, "PaloVerde"
	, "SP15_1"
	, "SP15_2")
var_2001 <- c(
	  "Price Hub"
	, "Trade Date"
	, "Wtd Avg Price $/MWh"
	, "Daily Volume MWh"
	)
var <- c(
	  "hub"
	, "date"
	, "price"
	, "volume"
	)
for (k in 1:length(hub_filelist_to2013)) {
	h <- hub_filelist_to2013[k]
	# load US ICE Wholesale Electricity Data
	ice_electric <- fread(file.path(pathtodata
		, "processed", "step1", "ICE"
		, paste0("ice_electric_2013before_", h, ".csv")))
	# check and change variable names
	setnames(ice_electric, var_2001, var)
	# select variables
	ice_electric <- ice_electric[, var, with = FALSE]
	ice_electric_to2013[[k]] <- ice_electric
}
ice_electric_to2013 <- rbindlist(ice_electric_to2013)

# rename hub names
unique(ice_electric_to2013$hub)
hub_Mass <- c(
	 "NEPOOL"
	,"NEPOOL Mass Hub"
	,"NEPOOL MH DA LMP"
	,"Nepool MH Da Lmp "
	,"Nepool MH Da Lmp"
	,"Nepool MH Da Lmp Peak"
	,"Nepool MH DA LMP"
	,"Nepool MH Da LMP Peak"
	)
ice_electric_to2013[hub %in% hub_Mass, hub := "Mass"]
hub_PJMWest <- c(
	 "PJM West"
	,"PJM Wh Real Time Peak"
	,"PJM-West Real Time Peak"
	,"PJM-Wh Real Time Peak"
	)
ice_electric_to2013[hub %in% hub_PJMWest, hub := "PJMWest"]
ice_electric_to2013[hub == "Mid Columbia Peak", hub := "MidC"]
ice_electric_to2013[hub == "Palo Verde", hub := "PaloVerde"]
hub_SP15 <- c(
	 "SP 15"
	,"SP-15 Peak"
	,"SP-15 Gen DA LMP Peak"
	)
ice_electric_to2013[hub %in% hub_SP15, hub := "SP15"]

# format date
ice_electric_to2013[, date := as.Date(date, format = "%Y/%m/%d")]

# format from 2014--2022
ice_electric_from2014 <- list()
var_2014 <- c(
	  "Price hub"
	, "Trade date"
	, "Wtd avg price $/MWh"
	, "Daily volume MWh"
	)
year_list <- 2014:2022
for (k in 1:length(year_list)) {
	y <- year_list[k]
	# load US ICE Wholesale Electricity Data
	ice_electric <- fread(file.path(pathtodata
		, "processed", "step1", "ICE"
		, paste0("ice_electric_", y, ".csv")), fill=TRUE)
	# check and change variable names
	setnames(ice_electric, var_2014, var)
	# select variables
	ice_electric <- ice_electric[, var, with = FALSE]
	# select the five hubs (Mass, PJM West, Mid-C, Palo Verde, and SP-15) and rename
	# # check hub names
	# print(unique(ice_electric$hub))
	ice_electric[hub %in% c("Nepool MH Da LMP Peak","Nepool MH DA LMP Peak"), hub := "Mass"]
	ice_electric[hub %in% c("PJM WH Real Time Peak","PJM-Wh Real Time Peak"), hub := "PJMWest"]
	ice_electric[hub %in% c("Mid C Peak","Mid Columbia Peak"), hub := "MidC"]
	ice_electric[hub %in% c("Palo Verde","Palo Verde Peak"), hub := "PaloVerde"]
	ice_electric[hub %in% c("SP15 EZ Gen DA LMP Peak","SP-15 Gen DA LMP Peak"), hub := "SP15"]
	ice_electric <- copy(ice_electric[hub %in% c("Mass","PJMWest","MidC","PaloVerde","SP15")])
	# format date
	if (y <= 2016) {
		ice_electric[, date := as.Date(date, format = "%Y/%m/%d")]
	} else {
		ice_electric[, date := as.Date(date, format = "%m-%d-%y")]
	}
	ice_electric_from2014[[k]] <- ice_electric
}
ice_electric_from2014 <- rbindlist(ice_electric_from2014)

# bind two formats
ice_electric_data <- rbind(ice_electric_to2013, ice_electric_from2014)

# remove duplicate entry
ice_electric_data <- unique.data.frame(ice_electric_data)

# fix two or more entries for same (hub,date)
ice_electric_data <- copy(ice_electric_data[, .(price = sum(price*volume)/sum(volume)
	, volume = sum(volume)), by = .(hub,date)])

# record year
ice_electric_data[, year := as.numeric(format(date, "%Y"))]

# inflation adjustment, base year = 2022
setnames(cpi_data, c("Year","Annual"), c("year","cpi"))
base_value <- cpi_data[year == 2022]$cpi
cpi_data[, adjust_2022 := base_value/cpi]
ice_electric_data <- merge(ice_electric_data, cpi_data[,.(year,adjust_2022)], by = "year")
ice_electric_data[, price := price*adjust_2022]
ice_electric_data[, adjust_2022 := NULL]

# save
fwrite(ice_electric_data, file.path(pathtodata, "processed", "ice_electric_data.csv"))
