# Brazil new energy auction data processing
# structural model data preparation

# @author: Konan Hara

# @date: September 13, 2026

# @software: R version 4.5.1
# clear environment
rm(list = ls())

# package loading
library(here)
library(data.table)
library(dplyr)
library(readxl)
library(lubridate)
library(zoo)

# path to data
pathtodata <- readLines(file(here("data.txt")), n = 1, warn = FALSE)
# configuration parameters
# dollar scale for Brazilian Real to USD conversion 
dscale <- 5 # by USD 1
# parent_category:
# load generation auction data
generation_auction <- fread(file.path(pathtodata
	, "processed", "generation_auction.csv"))
generation_auction_history <- copy(generation_auction)

# load generation auction participation data
generation_auction_participation <- fread(file.path(pathtodata
	, "processed", "generation_auction_participation.csv"))
# registered numbers for entry cost calculation in counterfactual analysis
generation_auction_registration <- copy(generation_auction_participation[stage == "registered"])
# use qualified numbers for number of participants
generation_auction_participation <- generation_auction_participation[stage == "qualified"]

# load generation auction firm data (for long-term (20 years) wind energy auctions)
generation_auction_represent <- fread(file.path(pathtodata
	, "processed", "generation_auction_represent.csv"))

# load spot market price data
spot_price_data <- fread(file.path(pathtodata
	, "processed", "spot_price_annual.csv"))

# load marginal cost of operation forecast data
mc_operation_forecast <- as.data.table(read_excel(file.path(pathtodata
	, "original", "Brazil", "marginal_cost_operation_forecast.xlsx")))

# use long-term (20 years) wind energy auctions
generation_auction <- generation_auction[
	source_type == "Wind Energy" & 
	contract_length > 15]
# round contract length
generation_auction[, contract_length := round(contract_length, 0)]
table(generation_auction$contract_length)
# -> only 20 years

# auction type
unique(generation_auction$auction_type)
# -> "LER" "LFA" "LEN"
# contract type
unique(generation_auction$contract_type)
# -> "QTD" "DIS"
# DIS (availability contract); QTD (quantity contract)
# may need sensitivity analysis for the contract types since descending-clock
# auctions are mostly QTD

# there are some duplicate plants
generation_auction[, dup := .N, by = .(plant_id)]
generation_dup <- copy(generation_auction[dup > 1])
generation_dup[order(plant_id)][,.(auction_name,auction_date,plant_id,contract_start_date,physical_assurance,contract_quantity_year,bid_price)]
# it seemed to be possible to participate in several auctions before starting operation
# only happens in pairs (27ºLEN,28ºLEN) and (33ºLEN,34ºLEN)
# subtract the first auction's physical assurance from the second auction's physical assurance
# since the firms knew that that portion would not be available for the second auction
generation_dup[, contract_quantity_year_previous := shift(contract_quantity_year, type = "lag"), by = .(plant_id)]
generation_dup[auction_name %in% c("28ºLEN","34ºLEN"), physical_assurance := physical_assurance - contract_quantity_year_previous]
generation_dup[, contract_quantity_year_previous := NULL]
# the other cases are pair (19ºLEN,26ºLEN) which seem to be a fresh re-participation
generation_auction <- rbind(generation_auction[dup == 1], generation_dup)
generation_auction[, dup := NULL]

# there could be projects with agreements since aneel_id == "004/2017" according to "BIDDING DOCUMENT FOR AUCTION"
# check "AUCTION RESULTS" to see the reference price
# -> all of them have the reference price of projects without agreements

# merge auction registration and participation information
cols <- c("number","total_capacity")
# registration
generation_auction_registration <- copy(generation_auction_registration[
	, lapply(.SD, sum), .SDcols = cols, by = .(auction_name,auction_category,source_category)])
generation_auction_registration <- copy(generation_auction_registration[source_category == "wind"])
generation_auction_registration[, source_category := NULL]
generation_auction <- merge(generation_auction, generation_auction_registration
	, by = c("auction_name","auction_category"))
setnames(generation_auction, "number", "number_registered")
setnames(generation_auction, "total_capacity", "total_capacity_registered")
# participation
generation_auction_participation <- copy(generation_auction_participation[
	, lapply(.SD, sum), .SDcols = cols, by = .(auction_name,auction_category,source_category)])
generation_auction_participation <- copy(generation_auction_participation[source_category == "wind"])
generation_auction_participation[, source_category := NULL]
generation_auction <- merge(generation_auction, generation_auction_participation
	, by = c("auction_name","auction_category"))
setnames(generation_auction, "number", "number_participant")

# merge parent firm
generation_auction <- merge(generation_auction, generation_auction_represent
	, by = c("aneel_id","plant_name","bidder_name"))

# Publicly observed wind-award history at the time of each auction. Count each
# project once at its first award and include it once contract delivery starts.
bidder_keys <- c("aneel_id", "plant_name", "bidder_name")
wind_history <- merge(
	generation_auction_history[source_type == "Wind Energy"]
	, generation_auction_represent
	, by = bidder_keys
	, allow.cartesian = TRUE
)
wind_history[, project_key := fifelse(
	!is.na(plant_id) & plant_id != ""
	, plant_id
	, paste(provider_id, plant_name, sep = "|")
)]
wind_history <- wind_history[
	order(auction_date), .SD[1], by = .(parent, project_key)
]
generation_auction[, prior_started_wind_capacity := vapply(
	seq_len(.N)
	, function(i) wind_history[
		parent == generation_auction$parent[i] &
			auction_date < generation_auction$auction_date[i] &
			contract_start_date <= generation_auction$auction_date[i]
		, sum(physical_assurance, na.rm = TRUE)
	]
	, numeric(1)
)]

# auction year
generation_auction[, year := as.numeric(format(auction_date, "%Y"))]
# decimal year from 2010
generation_auction[, year_decimal := decimal_date(auction_date) - 2010]
# descending-clock auction indicator (year >= 2017)
generation_auction[, descending_clock := as.numeric(year >= 2017)]

# auction indicator
# order by auction date
generation_auction[, auction_id :=
	as.numeric(as.factor(paste0(auction_date,auction_name,auction_category)))]

# exclude the two 2021 auctions conducted during the COVID-period market disruption
generation_auction <- generation_auction[!auction_id %in% c(23, 24)]

# regard yearly contract quantity as the bid quantity
setnames(generation_auction, "contract_quantity_year", "bid_quantity")
# regard physically assured quantity as the capacity
setnames(generation_auction, "capacity", "nameplate_capacity")
setnames(generation_auction, "physical_assurance", "capacity")

# Standardize the final CEG component from, for example, ".01" to ".1".
generation_auction[, plant_id := gsub("0([1-9])$", "\\1", plant_id)]

# number of winners
generation_auction[, number_winner := .N, by = .(auction_id)]

# procurement quantity
# sum of bid quantity of winners
generation_auction[, procurement_quantity := sum(bid_quantity), by = .(auction_id)]

# bid share
generation_auction[, bid_share := bid_quantity/capacity]

# bid price rank
# if tie in bid price, lower bid share ranks lower
setorder(generation_auction, "auction_id", "bid_price", -"bid_share", "plant_id")
generation_auction[, bid_price_rank := 1:.N, by = .(auction_id)]

# clearing price
generation_auction[, clearing_price := max(bid_price), by = .(auction_id)]

# bid share upper bound is less than 1 when quantity_minimum_unit == 1
# bid quantity upper bound
# sqrt(.Machine$double.eps) is for numerical precision problem in floor()
generation_auction[, bid_quantity_ub := floor(capacity/quantity_minimum_unit + sqrt(.Machine$double.eps))*quantity_minimum_unit]
# bid share upper bound
generation_auction[, bid_share_ub := bid_quantity_ub/capacity]

generation_auction[quantity_minimum_unit == 1]
generation_auction[year == 2009]

# bid share lower bound from ANEEL Energy Auction Notice ("BIDDING DOCUMENT FOR AUCTION")
# 0 for LER
generation_auction[auction_type == "LER", bid_share_lb_auction := 0.0]
# 0 until 2011
generation_auction[year <= 2011, bid_share_lb_auction := 0.0]
# 0.7 from 2012--2015 other than LER
generation_auction[year >= 2012 & year <= 2015 & auction_type != "LER", bid_share_lb_auction := 0.7]
# 0.3 from 2017
generation_auction[year >= 2017, bid_share_lb_auction := 0.3]
# bid share lower bound is not exact because of minimum quantity unit
# bid quantity lower bound
# sqrt(.Machine$double.eps) is for numerical precision problem in ceiling()
generation_auction[, bid_quantity_lb := ceiling(capacity*bid_share_lb_auction/quantity_minimum_unit - sqrt(.Machine$double.eps))*quantity_minimum_unit]
# bid share lower bound
generation_auction[, bid_share_lb := bid_quantity_lb/capacity]

# spot market price data
# consider the spot market prices are reflected to the bidding behavior at the end of the year
spot_price_data[, year_decimal := year - 2010 + 1]
auction_year_decimal <- data.table(year_decimal = unique(generation_auction$year_decimal))
spot_price_data <- rbind(spot_price_data[,.(year_decimal,spot_price)], auction_year_decimal, fill = TRUE)
spot_price_data[, spot_price := na.approx(spot_price, x = year_decimal, rule = 2)]
generation_auction <- merge(generation_auction, spot_price_data[,.(year_decimal,spot_price)], by = "year_decimal")

# marginal cost of operation forecast data
# average four regions
mc_operation_forecast[is.na(average), average := (SECO+S+NE+N)/4]
# use the forecasts of the latter five years (out of ten years)
mc_operation_forecast <- mc_operation_forecast[year >= year_forecast + 5]
# average forecasts of the latter five years
mc_operation_forecast <- mc_operation_forecast[
	, .(mc_operation = mean(average)), by = .(year_forecast)]
# linearly interpolated annual marginal cost of operation forecasts
# consider the forecasts are reflected to the bidding behavior at the end of the year
mc_operation_forecast[, year_decimal := year_forecast - 2010 + 1]
auction_year_decimal <- data.table(year_decimal = unique(generation_auction$year_decimal))
mc_operation_forecast <- rbind(mc_operation_forecast[,.(year_decimal,mc_operation)], auction_year_decimal, fill = TRUE)
mc_operation_forecast[, mc_operation := na.approx(mc_operation, x = year_decimal, rule = 2)]
generation_auction <- merge(generation_auction, mc_operation_forecast, by = "year_decimal")
# inflation factor
generation_auction[, inflation_factor := bid_price/bid_price_raw]
# adjust for rounding error to make the inflation factor unique to each auction
generation_auction[, inflation_factor := max(inflation_factor), by = .(auction_id)]
# inflation adjustment
generation_auction[, mc_operation := mc_operation*inflation_factor]

# marginal winner who may got partial bid quantity (only from aneel_id == "003/2018")
generation_auction[auction_date > "2018-08-30"
	, partial_bid_quantity := ifelse(bid_price == max(bid_price), 1, 0)
	, by = .(auction_id)]
generation_auction[is.na(partial_bid_quantity), partial_bid_quantity := 0]

# scale price
generation_auction[, bid_price := bid_price/dscale]
generation_auction[, clearing_price := clearing_price/dscale]
generation_auction[, mc_operation := mc_operation/dscale]

# select variables
var <- c(
	  "auction_id"
	, "auction_type"
	, "contract_type"
	, "descending_clock"
	, "year"
	, "year_decimal"
	, "lead_time"
	, "number_registered"
	, "number_participant"
	, "number_winner"
	, "total_capacity_registered"
	, "total_capacity"
	, "procurement_quantity"
	, "bid_share_lb_auction"
	, "bid_share_lb"
	, "bid_share_ub"
	, "quantity_minimum_unit"
	, "spot_price"
	, "mc_operation"
	, "plant_id"
	, "nameplate_capacity"
	, "capacity"
	, "bid_quantity"
	, "bid_share"
	, "bid_price"
	, "bid_price_rank"
	, "clearing_price"
	, "partial_bid_quantity"
	, "parent"
	, "prior_started_wind_capacity"
)
generation_auction <- copy(generation_auction[, var, with = FALSE])

# save bid data
fwrite(generation_auction, file.path(pathtodata, "processed", "generation_auction_wind20.csv"))
