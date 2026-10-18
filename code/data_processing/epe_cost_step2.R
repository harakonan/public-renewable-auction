# EPE power generation cost processing step 2

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

# load cost plot data
# median points
CAPEX_median <- fread(file.path(pathtodata, "original", "Brazil", "epe_cost/CAPEX_median.csv"))
OM_median <- fread(file.path(pathtodata, "original", "Brazil", "epe_cost/OM_median.csv"))
# violin plots
CAPEX_violin <- fread(file.path(pathtodata, "original", "Brazil", "epe_cost/CAPEX_violin.csv"))
OM_violin <- fread(file.path(pathtodata, "original", "Brazil", "epe_cost/OM_violin.csv"))

# configuration parameters
# dollar scale
dscale <- 5 # by USD 1
# observation period
obs_year <- 2010:2020

# scale prices (USD/kW)
CAPEX_median[, V2 := V2/dscale]
OM_median[, V2 := V2/dscale]
CAPEX_violin[, V2 := V2/dscale]
OM_violin[, V2 := V2/dscale]

# rename columns
CAPEX_median[, year := obs_year]
OM_median[, year := obs_year]
CAPEX_median[, V1 := NULL]
OM_median[, V1 := NULL]
setnames(CAPEX_median, "V2", "CAPEX")
setnames(OM_median, "V2", "OM")
setnames(CAPEX_violin, c("density","CAPEX"))
setnames(OM_violin, c("density","OM"))

# process median data
# combine CAPEX and OM
cost_median <- merge(CAPEX_median, OM_median, by = "year")

# process violin plot data
# convenient functions
# calculate mean from density f
find_mean <- function(f, ..., from = -Inf, to = Inf, subdivisions=100) {
	integrate(function(x) x * f(x, ...), from, to, subdivisions=subdivisions)
}
# calculate variance from density f
find_variance <- function(f, ..., from = -Inf, to = Inf, subdivisions=100) {
	mu <- find_mean(f, ..., from = from, to = to)$value
	integrate(function(x) (x - mu)^2 * f(x, ...), from, to, subdivisions=subdivisions)
}
# calculate mean and variance from violin plot
cost_data <- data.table(year = obs_year)
for (t in obs_year){
	# CAPEX
	# convert violin plot to density function
	data <- CAPEX_violin[density > (t-0.5) & density < t]
	data[, density := t - density]
	# scale density
	range <- c(min(data$CAPEX), max(data$CAPEX))
	f <- Vectorize(approxfun(data$CAPEX, data$density))
	density_scale <- integrate(f, range[1], range[2], subdivisions=100)$value
	data[, density := density/density_scale]
	# calculate mean and variance
	f <- Vectorize(approxfun(data$CAPEX, data$density))
	cost_data[year == t, CAPEX_mean := find_mean(f, from = range[1], to = range[2])$value]
	cost_data[year == t, CAPEX_variance := find_variance(f, from = range[1], to = range[2])$value]

	# OM
	# convert violin plot to density function
	data <- OM_violin[density > (t-0.5) & density < t]
	data[, density := t - density]
	# scale density
	range <- c(min(data$OM), max(data$OM))
	f <- Vectorize(approxfun(data$OM, data$density))
	density_scale <- integrate(f, range[1], range[2], subdivisions=100)$value
	data[, density := density/density_scale]
	# calculate mean and variance
	f <- Vectorize(approxfun(data$OM, data$density))
	cost_data[year == t, OM_mean := find_mean(f, from = range[1], to = range[2])$value]
	cost_data[year == t, OM_variance := find_variance(f, from = range[1], to = range[2])$value]
}

# save cost data
fwrite(cost_data, file.path(pathtodata, "processed", "epe_cost_data.csv"))
