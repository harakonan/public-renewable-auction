# US ICE Wholesale Electricity Data processing
# monthly/annual spot price calculation

# @author: Konan Hara

# @date: November 18, 2024

# @software: R version 4.4.1

# clear environment
rm(list = ls())

# package loading
library(here)
library(data.table)
library(zoo)

# path to data
pathtodata <- readLines(file(here("data.txt")), n = 1, warn = FALSE)

# load US ICE Wholesale Electricity Data
ice_electric_data <- fread(file.path(pathtodata, "processed", "ice_electric_data.csv"))

# monthly average for each hub
ice_electric_data[, month := as.yearmon(date)]
ice_electric_data[, month := as.Date(month)]
spot_price_monthly_us_hub <- ice_electric_data[, .(spot_price = mean(price)), by = .(hub,month)]
# average across hubs
spot_price_monthly_us <- spot_price_monthly_us_hub[, .(spot_price = mean(spot_price)), by = .(month)]
# bind them in one
spot_price_monthly_us[, hub := "Average"]
setcolorder(spot_price_monthly_us, names(spot_price_monthly_us_hub))
spot_price_monthly_us <- rbind(spot_price_monthly_us, spot_price_monthly_us_hub)

# annual average for each hub
spot_price_annual_us_hub <- ice_electric_data[, .(spot_price = mean(price)), by = .(hub,year)]
# average across hubs
spot_price_annual_us <- spot_price_annual_us_hub[, .(spot_price = mean(spot_price)), by = .(year)]
# bind them in one
spot_price_annual_us[, hub := "Average"]
setcolorder(spot_price_annual_us, names(spot_price_annual_us_hub))
spot_price_annual_us <- rbind(spot_price_annual_us, spot_price_annual_us_hub)

# save
fwrite(spot_price_monthly_us, file.path(pathtodata, "processed", "spot_price_monthly_us.csv"))
fwrite(spot_price_annual_us, file.path(pathtodata, "processed", "spot_price_annual_us.csv"))
