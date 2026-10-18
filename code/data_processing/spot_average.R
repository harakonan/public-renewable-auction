# Brazil new energy auction data processing
# calculate average spot prices

# @author: Konan Hara

# @date: April 20, 2025

# @software: R version 4.4.1

# clear environment
rm(list = ls())

# package loading
library(here)
library(data.table)
library(dplyr)
library(readxl)

# path to data
pathtodata <- readLines(file(here("data.txt")), n = 1, warn = FALSE)

# configuration parameters
# dollar scale
dscale <- 5 # by USD 1

# load spot market price data
spot_price_monthly <- as.data.table(read_excel(file.path(pathtodata
	, "original", "Brazil", "monthly_average_price_200101_202212.xls")))

# average price across regions
spot_price_monthly[, spot_price := (SUDESTE+SUL+NORDESTE+NORTE)/4]
# use annual average price
setnames(spot_price_monthly, "MES", "month")
spot_price_monthly[, year := as.numeric(format(month,"%Y"))]
spot_price_annual <- spot_price_monthly[, .(spot_price = mean(spot_price)), by = .(year)]

# spot market monthly data
spot_price_monthly <- copy(spot_price_monthly[,.(month,spot_price)])

# scale price
spot_price_annual[, spot_price := spot_price/dscale]
spot_price_monthly[, spot_price := spot_price/dscale]

# save
# annual data is available up to 2022
spot_price_annual <- copy(spot_price_annual[year <= 2022])
fwrite(spot_price_annual, file.path(pathtodata, "processed", "spot_price_annual.csv"))
fwrite(spot_price_monthly, file.path(pathtodata, "processed", "spot_price_monthly.csv"))
