# Brazil central bank target interest rate data processing
# 1. select variables
# 2. tidy up

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

# load Brazil central bank target interest rate data
base_interest_rate <- as.data.table(read_excel(
	  file.path(pathtodata, "original", "Brazil", "central_bank", "base_interest_rates_history.xlsx")
	, sheet = "base_interest_rates_history"
	, col_names = FALSE
	, skip = 852))

# trim after the table
base_interest_rate[`...1` == "1º", which = TRUE]
# -> 262
base_interest_rate <- copy(base_interest_rate[1:262])

# select variables
# V4: period
# V5: target interest rate
base_interest_rate <- copy(base_interest_rate[,.(`...4`,`...5`)])
setnames(base_interest_rate, c("date","rate"))

# collect start date
base_interest_rate[, date := substr(date, 1, 10)]

# tidy up
base_interest_rate[, date := as.Date(date, format = "%d/%m/%Y")]
base_interest_rate[, rate := as.numeric(gsub(",", ".",rate))]

# save
fwrite(base_interest_rate, file.path(pathtodata
	, "processed", "base_interest_rate.csv"))
