# data for Brazil LEE auction
# add technology designation information from `CCEE_auction_description_LEE/`
# @author: Konan Hara

# @date: November 23, 2024

# @software: R version 4.4.1

# clear environment
rm(list = ls())

# package loading
library(here)
library(data.table)

# path to data
pathtodata <- readLines(file(here("data.txt")), n = 1, warn = FALSE)

# load Brazil auction data
generation_auction <- fread(file.path(pathtodata, "processed", "generation_auction.csv"))

# existing energy auctions (LEE)
generation_auction_lee <- copy(generation_auction[auction_type == "LEE"])

# auction list
unique.data.frame(generation_auction_lee[,.(auction_name, auction_category, source_type)])
# insert technology from `CCEE_auction_description_LEE/`
generation_auction_lee[auction_name == "01ºLEE", technology := "neutral"]
generation_auction_lee[auction_name == "02ºLEE", technology := "neutral"]
generation_auction_lee[auction_name == "03ºLEE", technology := "neutral"]
generation_auction_lee[auction_name == "04ºLEE", technology := "neutral"]
generation_auction_lee[auction_name == "05ºLEE", technology := "neutral"]
generation_auction_lee[auction_name == "08ºLEE", technology := "neutral"]
generation_auction_lee[auction_name == "08ºLEE", technology := "neutral"]
generation_auction_lee[auction_name == "09ºLEE" &
		auction_category == "2011-Q03"
	, technology := "hydro"]
generation_auction_lee[auction_name == "09ºLEE" &
		auction_category == "2011-D03"
	, technology := "thermo"]
generation_auction_lee[auction_name == "10ºLEE", technology := "neutral"]
generation_auction_lee[auction_name == "12ºLEE", technology := "neutral"]
generation_auction_lee[auction_name == "13ºLEE", technology := "neutral"]
generation_auction_lee[auction_name == "14ºLEE" &
		auction_category == "DIS-03ANOS-2015"
	, technology := "thermo"]
generation_auction_lee[auction_name == "14ºLEE" &
		auction_category == "QTD2-03ANOS-2015"
	, technology := "neutral"]
generation_auction_lee[auction_name == "15ºLEE" &
		auction_category == "DIS2-03ANOS"
	, technology := "thermo"]
generation_auction_lee[auction_name == "15ºLEE" &
		auction_category == "DIS3-01ANO"
	, technology := "thermo"]
generation_auction_lee[auction_name == "15ºLEE" &
		auction_category == "QTD-03ANOS"
	, technology := "neutral"]
generation_auction_lee[auction_name == "16ºLEE", technology := "neutral"]
generation_auction_lee[auction_name == "17ºLEE", technology := "neutral"]
generation_auction_lee[auction_name == "18ºLEE", technology := "neutral"]
generation_auction_lee[auction_name == "19ºLEE", technology := "neutral"]
generation_auction_lee[auction_name == "20ºLEE" &
		auction_category == "DIS-2020-2A"
	, technology := "thermo"]
generation_auction_lee[auction_name == "20ºLEE" &
		auction_category == "QTD-2020-2A"
	, technology := "neutral"]
generation_auction_lee[auction_name == "21ºLEE", technology := "neutral"]
generation_auction_lee[auction_name == "22ºLEE" &
		auction_category == "DIS-2021-2A"
	, technology := "thermo"]
generation_auction_lee[auction_name == "22ºLEE" &
		auction_category == "QTD-2021-2A"
	, technology := "neutral"]
generation_auction_lee[auction_name == "23ºLEE", technology := "thermo"]
generation_auction_lee[auction_name == "24ºLEE", technology := "thermo"]
generation_auction_lee[auction_name == "25ºLEE", technology := "neutral"]
generation_auction_lee[auction_name == "26ºLEE", technology := "neutral"]
generation_auction_lee[auction_name == "27ºLEE", technology := "neutral"]

# save
fwrite(generation_auction_lee, file.path(pathtodata, "processed"
	, "generation_auction_lee.csv"))
