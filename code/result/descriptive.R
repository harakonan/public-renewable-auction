# Brazil new energy auction descriptive statistics

# @author: Konan Hara

# @date: September 13, 2026

# @software: R version 4.5.3
# clear environment
rm(list = ls())

# package loading
library(here)
library(data.table)
library(dplyr)
library(kableExtra)
library(ggplot2)
library(lubridate)
library(fixest)

# path to data
pathtodata <- readLines(file(here("data.txt")), n = 1, warn = FALSE)

# output path
fdir <- here("output", "figure")
tdir <- here("output", "table")
sdir <- here("output", "snt")

# create target output directories if they do not exist
dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
dir.create(tdir, showWarnings = FALSE, recursive = TRUE)
dir.create(sdir, showWarnings = FALSE, recursive = TRUE)

# Load the preferred Step 1 and descriptive sample of 20-year PPA wind
# auctions. The two 2021 COVID-period auctions are excluded upstream. The two
# sealed-bid auctions with missing participant counts remain in this file
# because participant counts are not required for Step 1.
auction_data <- fread(file.path(pathtodata, "processed", "generation_auction_wind20.csv"))

# load EPE wind turbine cost data
epe_cost_data <- fread(file.path(pathtodata, "processed", "epe_cost_data.csv"))

# load central bank base interest rate data
base_interest_rate_data <- fread(file.path(pathtodata, "processed", "base_interest_rate.csv"))

# load spot market price data
spot_price_annual <- fread(file.path(pathtodata, "processed", "spot_price_annual.csv"))
spot_price_monthly <- fread(file.path(pathtodata, "processed", "spot_price_monthly.csv"))
# use IDate
spot_price_monthly[, month := as.IDate(month)]

# load US spot market price data
spot_price_annual_us <- fread(file.path(pathtodata, "processed", "spot_price_annual_us.csv"))
spot_price_monthly_us <- fread(file.path(pathtodata, "processed", "spot_price_monthly_us.csv"))
# use average across the hubs
spot_price_annual_us <- copy(spot_price_annual_us[hub == "Average"])
spot_price_annual_us[, hub := NULL]
spot_price_monthly_us <- copy(spot_price_monthly_us[hub == "Average"])
spot_price_monthly_us[, hub := NULL]

# load Brazil short-term purchase agreement auction (LEE) data
auction_lee <- fread(file.path(pathtodata, "processed", "generation_auction_lee.csv"))

# configuration parameters
# purchase agreement length
T <- 20
# discount factor
delta <- 0.95
# delta tilde at lead time = 0
delta_tilde_0 <- sum(delta^(0:(T - 1))) / T

# function to format
# value
value_format <- function(value, d) {
	paste0("$", sprintf(paste0("%.", d, "f"), value), "$")
}
# estimates with SEs
est_format <- function(coef, se, d) {
	as.vector(rbind(
		  value_format(coef, d)
		, paste0("$(", sprintf(paste0("%.", d, "f"), se), ")$")
	))
}
# snt
value_format_snt <- function(x, d, s) {
	paste0(sprintf(paste0("%.", d, "f"), x * s), "%")
}
value_format_snt_pct <- function(x, d, s) {
	paste0(sprintf(paste0("%.", d, "f"), x * s), "\\%%")
}
write_snt <- function(x, filename, d = 2, s = 1, pct = FALSE) {
	val <- if (pct) value_format_snt_pct(x, d, s) else value_format_snt(x, d, s)
	cat(val, file = file.path(sdir, filename))
}
write_snt_raw <- function(x, filename) {
	cat(paste0(x, "%"), file = file.path(sdir, filename))
}
fmt_num <- function(x, digits) {
	if (length(x) == 0 || is.na(x)) return("$NA$")
	sprintf(paste0("$%0.", digits, "f$"), x)
}
fmt_range <- function(lo, hi, digits) {
	if (is.na(lo) || is.na(hi)) return("$[NA, NA]$")
	paste0(
		  "$["
		, sprintf(paste0("%0.", digits, "f"), lo)
		, ", "
		, sprintf(paste0("%0.", digits, "f"), hi)
		, "]$"
	)
}

################################################################################
# overall summary statistics
################################################################################

# maximum contract length in technology neutral short-term purchase agreement auctions since 2011 (years)
value <- max(auction_lee[auction_date >= "2011-01-01" & technology == "neutral"]$contract_length)
sprintf(paste0("maximum contract length in technology neutral short-term purchase agreement auctions since 2011 (years): "
	,"$%.", 1, "f$"), value)

# ANEEL auction-SIGA linkage for construction outcomes. These objects are
# descriptive only and do not enter estimation. Collapse SIGA's ownership rows
# to one row per physical wind plant.
generation_data <- fread(file.path(
	pathtodata, "processed", "generation_data.csv"
))[type == "wind"]
generation_data <- unique(generation_data[, .(
	  plant_id
	, phase
	, operation_start
	, nameplate_capacity_siga = capacity_kw/1e3
)])
stopifnot(generation_data[, max(uniqueN(phase), na.rm = TRUE), by = plant_id]$V1 == 1)

# quantity and availability contracts
# number of auctions by contract type
data <- distinct(auction_data[, .(auction_id, contract_type)])
auction_N_qtd <- data[contract_type == "QTD", .N]
auction_N_dis <- data[contract_type == "DIS", .N]
print(paste0("number of auctions: quantity contracts, ", auction_N_qtd))
print(paste0("number of auctions: availability contracts, ", auction_N_dis))
write_snt_raw(auction_N_qtd, "auction_N_qtd.tex")
write_snt_raw(auction_N_dis, "auction_N_dis.tex")
# time trend?
distinct(auction_data[, .(auction_id, contract_type, year)])

# number of bids
val <- auction_data[, .N]
write_snt(val, "bidder_N_all.tex", 0)
# remove observations with the possibility of partial bid quantity
# number of such cases
val <- auction_data[partial_bid_quantity == 1, .N]
write_snt(val, "bidder_N_partial_contract.tex", 0)
auction_data <- auction_data[partial_bid_quantity == 0]
stopifnot(
	  auction_data[, all(is.finite(bid_share))]
	, auction_data[, all(bid_share >= 0)]
	, auction_data[, all(bid_share <= bid_share_ub + sqrt(.Machine$double.eps))]
)

# number of bids
# total
write_snt_raw(auction_data[, .N], "bidder_N.tex")
# sealed-bid auctions
write_snt_raw(auction_data[descending_clock == 0, .N], "bidder_N_sb.tex")
# sealed-bid auctions with non-missing number of participants
write_snt_raw(auction_data[descending_clock == 0 & !is.na(number_participant), .N]
	, "bidder_N_sb_nomissing.tex")
# descending-clock auctions
write_snt_raw(auction_data[descending_clock == 1, .N], "bidder_N_dc.tex")

# total capacity (GW)
val <- auction_data[, sum(capacity)]/1e3
write_snt(val, "bidder_total_capacity.tex", 1)

################################################################################
# auction-SIGA wind construction linkage
################################################################################

# Treat each physical plant once. For projects awarded in multiple estimation-
# sample auctions, assign the project to its first retained auction cohort.
winner_projects <- copy(auction_data[order(auction_id), .SD[1], by = plant_id])
winner_projects <- merge(
	winner_projects
	, generation_data
	, by = "plant_id"
	, all.x = TRUE
)
winner_projects[, observed_operating :=
	!is.na(phase) & phase == "operation" & year(operation_start) <= 2025]

# Construction outcomes by first auction cohort.
winner_construction_cohort <- winner_projects[, .(
	  share = mean(observed_operating)
	, projects = .N
), by = .(cohort = year)]

# Conservative round-number lower bound for all cohorts outside 2011--2013.
write_snt(
	floor(winner_construction_cohort[
		!cohort %in% 2011:2013, min(share)
	] * 10)/10
	, "winner_project_operating_prop_other_years_floor_10.tex"
	, 0
	, 100
	, TRUE
)

# Average construction rate for the lower-construction 2011--2013 cohorts.
write_snt(
	winner_projects[year %in% 2011:2013, mean(observed_operating)]
	, "winner_project_operating_prop_2011_2013.tex"
	, 0
	, 100
	, TRUE
)

# Among wind plants operating in SIGA, identify matches to winners in the
# 2009--2019 analysis sample. Use 2011 as the first plausible commissioning
# year for projects in the 2009 auction and end the sample in 2023.
estimation_winner_ids <- unique(auction_data$plant_id)
siga_operating <- copy(generation_data[
	phase == "operation" & year(operation_start) >= 2011 &
		year(operation_start) <= 2023
])
siga_operating[, `:=`(
	  cohort = year(operation_start)
	, estimation_winner = plant_id %in% estimation_winner_ids
)]

# Main-text statistic for the central commissioning period.
write_snt(
	siga_operating[cohort <= 2020, mean(estimation_winner)]
	, "operating_project_auction_winner_prop_2011_2020.tex"
	, 0
	, 100
	, TRUE
)

operating_origin_cohort <- siga_operating[, .(
	share = mean(estimation_winner)
), by = cohort]

# Round the post-2020 cohort range to the nearest 10 percentage points for the
# concise discussion in the main text.
post_2020_linkage_range <- operating_origin_cohort[
	cohort %in% 2021:2023, range(share)
]
write_snt(
	round(post_2020_linkage_range[1] * 10)/10
	, "operating_project_auction_winner_prop_2021_2023_min.tex"
	, 0
	, 100
	, TRUE
)
write_snt(
	round(post_2020_linkage_range[2] * 10)/10
	, "operating_project_auction_winner_prop_2021_2023_max.tex"
	, 0
	, 100
	, TRUE
)

################################################################################
# auction-level summary statistics
################################################################################

variable <- c(
  "lead_time",
  "number_participant",
  "total_capacity",
  "number_winner",
  "number_winning_firm",
  "procurement_quantity"
)

variable_int <- c(
  "number_participant",
  "number_winner",
  "number_winning_firm"
)

variable_name <- c(
  "Lead Time (years)",
  "Number of Participants",
  "Participants Total Capacity (GW)",
  "Number of Winners",
  "Number of Winning Firms",
  "Procurement Capacity (MW)"
)

# rescale total capacity from MW to GW
auction_data[, total_capacity := total_capacity/1e3]
# count representative parent companies among sample winners in each auction
auction_data[, number_winning_firm := uniqueN(parent), by = auction_id]

d <- 1

auction_data_a <- copy(distinct(auction_data[
  , c("auction_id","descending_clock","year",variable), with = FALSE
]))

# summary table by auction type
compute_stats <- function(dt, ua_label) {
  rows <- vector("list", length(variable))
  for (k in seq_along(variable)) {
    v <- variable[k]
    vname <- variable_name[k]
    digits_here <- if (v %in% variable_int) 0 else d
    x <- dt[[v]]

    n_nonmiss <- sum(!is.na(x))
    m  <- if (n_nonmiss > 0) mean(x, na.rm = TRUE) else NA_real_
    s  <- if (n_nonmiss > 0) sd(x, na.rm = TRUE) else NA_real_
    lo <- if (n_nonmiss > 0) suppressWarnings(min(x, na.rm = TRUE)) else NA_real_
    hi <- if (n_nonmiss > 0) suppressWarnings(max(x, na.rm = TRUE)) else NA_real_

    rows[[k]] <- data.table(
      auction_type = ua_label,
      var_order = k,
      variable = vname,
      N = n_nonmiss,
      Mean = fmt_num(m, digits_here),
      SD = fmt_num(s, digits_here),
      Range = fmt_range(lo, hi, digits_here)
    )
  }
  rbindlist(rows)
}

# split data once to keep auction blocks contiguous
sb_dt <- auction_data_a[descending_clock == 0]
dc_dt <- auction_data_a[descending_clock == 1]

sb_res <- compute_stats(sb_dt, "Sealed-Bid")[order(var_order)]
dc_res <- compute_stats(dc_dt, "Descending-Clock")[order(var_order)]

res <- rbindlist(list(sb_res, dc_res), use.names = TRUE)
res[, var_order := NULL]

# overall auction counts and years
auction_N_sb <- sb_dt[, .N]
auction_year_min_sb <- sb_dt[, min(year)]
auction_year_max_sb <- sb_dt[, max(year)]
auction_N_dc <- dc_dt[, .N]
auction_year_min_dc <- dc_dt[, min(year)]
auction_year_max_dc <- dc_dt[, max(year)]

write_snt_raw(auction_data_a[, .N], "auction_N.tex")
write_snt_raw(auction_data_a[, min(year)], "auction_year_min.tex")
write_snt_raw(auction_data_a[, max(year)], "auction_year_max.tex")
write_snt_raw(auction_N_sb, "auction_N_sb.tex")
write_snt_raw(sb_dt[!is.na(number_participant), .N], "auction_N_sb_nomissing.tex")
write_snt_raw(auction_year_min_sb, "auction_year_min_sb.tex")
write_snt_raw(auction_year_max_sb, "auction_year_max_sb.tex")
write_snt_raw(auction_N_dc, "auction_N_dc.tex")
write_snt_raw(auction_year_min_dc, "auction_year_min_dc.tex")
write_snt_raw(auction_year_max_dc, "auction_year_max_dc.tex")

# create latex file
n_sb <- nrow(sb_res)
n_dc <- nrow(dc_res)
kbl(
  res[, .(variable, N, Mean, SD, Range)],
  format = "latex",
  col.names = c("", "N", "Mean", "SD", "Range"),
  align = "lcccc",
  escape = FALSE,
  booktabs = TRUE,
  position = "htb!"
) %>%
  pack_rows(
    group_label = paste0(auction_N_sb
    	, " Sealed-Bid Auctions (", auction_year_min_sb, "--", auction_year_max_sb, ")"),
    start_row = 1,
    end_row = n_sb,
    indent = FALSE
  ) %>%
  pack_rows(
    group_label = paste0(auction_N_dc
    	, " Descending-Clock Auctions (", auction_year_min_dc, "--", auction_year_max_dc, ")"),
    start_row = n_sb + 1,
    end_row = n_sb + n_dc,
    indent = FALSE
  ) %>%
  save_kable(file.path(tdir, "summary_auction.tex"))


################################################################################
# bidder-level summary statistics
################################################################################

# proportion of bid share partial allocation (percent)
auction_data[, partial := 1]
auction_data[bid_share > (bid_share_ub-sqrt(.Machine$double.eps)) |
		(bid_share > (bid_share_lb-sqrt(.Machine$double.eps)) & bid_share <= (bid_share_lb+sqrt(.Machine$double.eps)))
	, partial := 0]
write_snt(auction_data[, sum(partial == 1)/.N], "bidder_prop_partial.tex", 1, 100, TRUE)
# # by contract type
# auction_data[, (sum(partial == 1)/.N)*100, by = .(contract_type)]

# bidder-level variables and labels
b_variables <- c("bid_share", "bid_price", "capacity")
b_variable_name <- c("Bid Share", "Purchase Agreement Price (\\$/MWh)", "Capacity (MW)")

# display precision
digits_map <- c("bid_share" = 2, "bid_price" = 1, "capacity" = 1)

# quantile ranges for bidder-level summaries
fmt_quantile <- function(x, probs, digits) {
	if (length(x) == 0 || all(is.na(x))) return("$[NA, NA]$")
	q <- suppressWarnings(quantile(x, probs = probs, na.rm = TRUE, type = 7))
	if (any(is.na(q))) return("$[NA, NA]$")
	fmt_range(q[[1]], q[[2]], digits)
}

# summary table by auction type
compute_bidder_stats <- function(dt, ua_label) {
  rows <- vector("list", length(b_variables))
  for (k in seq_along(b_variables)) {
    v <- b_variables[k]
    vname <- b_variable_name[k]
    digits_here <- digits_map[[v]]

    x <- dt[[v]]
    n_nonmiss <- sum(!is.na(x))
    m  <- if (n_nonmiss > 0) mean(x, na.rm = TRUE) else NA_real_
    s  <- if (n_nonmiss > 0) sd(x, na.rm = TRUE) else NA_real_

    rows[[k]] <- data.table(
      auction_type = ua_label,
      var_order = k,                     # keep the specified variable order
      variable = vname,
      N = n_nonmiss,
      Mean = fmt_num(m, digits_here),
      SD = fmt_num(s, digits_here),
      IQR = fmt_quantile(x, c(0.25, 0.75), digits_here),
      P5P95 = fmt_quantile(x, c(0.05, 0.95), digits_here)
    )
  }
  rbindlist(rows)
}

# use all winning bids
sb_dt <- auction_data[descending_clock == 0]
dc_dt <- auction_data[descending_clock == 1]

sb_res <- compute_bidder_stats(sb_dt, "Sealed-Bid")[order(var_order)]
dc_res <- compute_bidder_stats(dc_dt, "Descending-Clock")[order(var_order)]

res <- rbindlist(list(sb_res, dc_res), use.names = TRUE)
res[, var_order := NULL]

# number of auctions and year span
auction_data_a <- distinct(auction_data[, .(auction_id, descending_clock, year)])
auction_N_sb <- auction_data_a[descending_clock == 0, .N]
auction_year_min_sb <- auction_data_a[descending_clock == 0, min(year)]
auction_year_max_sb <- auction_data_a[descending_clock == 0, max(year)]
auction_N_dc <- auction_data_a[descending_clock == 1, .N]
auction_year_min_dc <- auction_data_a[descending_clock == 1, min(year)]
auction_year_max_dc <- auction_data_a[descending_clock == 1, max(year)]

# create latex file
n_sb <- nrow(sb_res)
n_dc <- nrow(dc_res)
kbl(
  res[, .(variable, N, Mean, SD, IQR, P5P95)],
  format = "latex",
  col.names = c("", "N", "Mean", "SD", "IQR", "[P5, P95]"),
  align = "lccccc",
  escape = FALSE,
  booktabs = TRUE,
  position = "htb!"
) %>%
  pack_rows(
    group_label = paste0(auction_N_sb, " Sealed-Bid Auctions (", auction_year_min_sb, "--", auction_year_max_sb, ")"),
    start_row = 1,
    end_row = n_sb,
    indent = FALSE
  ) %>%
  pack_rows(
    group_label = paste0(auction_N_dc, " Descending-Clock Auctions (", auction_year_min_dc, "--", auction_year_max_dc, ")"),
    start_row = n_sb + 1,
    end_row = n_sb + n_dc,
    indent = FALSE
  ) %>%
  save_kable(file.path(tdir, "summary_bidder.tex"))

res_snt <- res[auction_type == "Sealed-Bid" & variable == "Bid Share"]
write_snt_raw(res_snt[, Mean], "bid_share_mean_sb.tex")
write_snt_raw(res_snt[, IQR], "bid_share_iqr_sb.tex")
write_snt_raw(res_snt[, P5P95], "bid_share_p5p95_sb.tex")
res_snt <- res[auction_type == "Descending-Clock" & variable == "Bid Share"]
write_snt_raw(res_snt[, Mean], "bid_share_mean_dc.tex")
write_snt_raw(res_snt[, IQR], "bid_share_iqr_dc.tex")
write_snt(mean(auction_data$capacity), "bidder_capacity_mean.tex", 0)

# bid share and price bin scatter
data <- copy(auction_data)

# create decile bins
data[, price_decile := cut(bid_price
	, breaks = quantile(bid_price, probs = seq(0, 1, 0.1), na.rm = TRUE)
	, include.lowest = TRUE)]

# bid price and bid share summaries by decile
data_bid <- data[!is.na(price_decile), .(
		  bid_price_mean = mean(bid_price, na.rm = TRUE)
		, bid_share_mean = mean(bid_share, na.rm = TRUE)
		, bid_share_q1 = quantile(bid_share, 0.25, na.rm = TRUE)
		, bid_share_q3 = quantile(bid_share, 0.75, na.rm = TRUE)
	), by = .(price_decile)]

# reshape for plotting
data_bid <- melt(data_bid
	, id.vars = c("price_decile", "bid_price_mean")
	, measure.vars = paste0("bid_share_", c("mean", "q1", "q3")))

# label name
data_bid[, variable := factor(variable
	, levels = paste0("bid_share_", c("mean","q1","q3"))
	, labels = c("Mean","Interquartile Range","Q3"))]

# figure
fig <- ggplot() +
	geom_line(aes(x=bid_price_mean, y=value, linetype=variable)
		, data_bid[variable %in% c("Mean","Interquartile Range")]
		, linewidth = 0.6) +
	geom_point(aes(x=bid_price_mean, y=value)
		, data_bid[variable == "Mean"]
		, size = 1.5) +
	geom_line(aes(x=bid_price_mean, y=value)
		, data_bid[variable == "Q3"]
		, linetype=3
		, linewidth = 0.6) +
	scale_linetype_manual(values = c(1,3)) +
	xlab("Purchase Agreement Price ($/MWh)") +
	ylab("Bid Share") +
	ylim(0,1) +
	theme_bw() +
	theme(legend.title=element_blank()
		, legend.position = c(0.98, 0.02)
		, legend.justification = c(1, 0)
		, legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
		, legend.key = element_blank()
	)
ggsave(file.path(fdir, "bidshare_price_bin.png")
	, plot = fig, width = 8, height = 5, units = "cm")

################################################################################
# preliminary regressions
################################################################################

# digits to the right of the decimal point to be rounded
d <- 4

# Prior-started wind capacity is expressed in units of 100 MW in the reduced
# form. The structural model applies its log transformation separately.
auction_data[, prior_started_wind_capacity_100 :=
	prior_started_wind_capacity/100]

# bid share ~ price
est <- list()
m <- feols(
	bid_share ~ bid_price
	, auction_data
	, cluster = ~auction_id + parent
)
est[[1]] <- c(est_format(coef(m)[names(coef(m)) != "(Intercept)"]
		, se(m)[names(se(m)) != "(Intercept)"], d)
	, rep("", 4)
	, m$nobs
	, length(unique(auction_data$auction_id))
	, length(unique(auction_data$parent))
	, value_format(r2(m)["r2"], 2))

# + expected marginal cost control
m <- feols(
	bid_share ~ bid_price + mc_operation
	, auction_data
	, cluster = ~auction_id + parent
)
est[[2]] <- c(est_format(coef(m)[names(coef(m)) != "(Intercept)"]
		, se(m)[names(se(m)) != "(Intercept)"], d)
	, rep("", 2)
	, m$nobs
	, length(unique(auction_data$auction_id))
	, length(unique(auction_data$parent))
	, value_format(r2(m)["r2"], 2))

# + prior started wind capacity
m <- feols(
	bid_share ~ bid_price + mc_operation + prior_started_wind_capacity_100
	, auction_data
	, cluster = ~auction_id + parent
)
est[[3]] <- c(est_format(coef(m)[names(coef(m)) != "(Intercept)"]
		, se(m)[names(se(m)) != "(Intercept)"], d)
	, m$nobs
	, length(unique(auction_data$auction_id))
	, length(unique(auction_data$parent))
	, value_format(r2(m)["r2"], 2))
write_snt(se(m)["prior_started_wind_capacity_100"]
	, "prior_capacity_col3_se.tex", 4)

# Background robustness regression for the contract-type discussion in the
# institutional-details footnote. Table 3 retains the three specifications
# above; this regression adds a quantity-contract indicator to Column (3).
m_contract <- feols(
	bid_share ~ bid_price + mc_operation + prior_started_wind_capacity_100 + contract_type
	, auction_data
	, cluster = ~auction_id + parent
)
contract_term <- "contract_typeQTD"
write_snt(coef(m_contract)[contract_term], "contract_type_coef.tex", 3)
write_snt(se(m_contract)[contract_term], "contract_type_se.tex", 3)

# Background diagnostic for the capacity-invariant portfolio rule. This
# regression adds the project's own capacity to Column (3), while retaining
# the same controls and two-way clustered standard errors.
m_capacity <- feols(
	bid_share ~ bid_price + mc_operation + prior_started_wind_capacity_100 + capacity
	, auction_data
	, cluster = ~auction_id + parent
)
write_snt(coef(m_capacity)["capacity"], "capacity_bidshare_coef.tex", 4)
write_snt(se(m_capacity)["capacity"], "capacity_bidshare_se.tex", 4)
write_snt(coef(m_capacity)["prior_started_wind_capacity_100"]
	, "prior_capacity_with_capacity_coef.tex", 4)
write_snt(se(m_capacity)["prior_started_wind_capacity_100"]
	, "prior_capacity_with_capacity_se.tex", 4)

# table
table <- data.table(
	  rowname = c(
	  	  "Purchase Agreement Price (\\$/MWh)",""
	  	, "Forecasted Marginal Cost of Electricity (\\$/MWh)",""
		, "Prior Started Wind Capacity (100 MW)",""
	  	, "Number of observations"
	  	, "Number of auctions"
		, "Number of firms"
	  	, "$R^2$"
	)
	, do.call(cbind, est)
)

kbl(table
	, format = "latex"
	, col.names = c("", "(1)", "(2)", "(3)")
	, align = "lcccc"
	, escape = FALSE
	, booktabs = TRUE
	, position = "htb!"
	, linesep = c(rep("", 5), "\\addlinespace")
) %>%
save_kable(file.path(tdir, "preliminary_regression.tex"))

################################################################################
# spot price
################################################################################

# digits to the right of the decimal point to be rounded
d <- 2

val <- min(spot_price_annual$year)
write_snt(val, "spot_price_year_min.tex", 0)
val <- max(spot_price_annual$year)
write_snt(val, "spot_price_year_max.tex", 0)
# average of all periods
val <- mean(spot_price_annual$spot_price)
write_snt(val, "spot_price_mean_all.tex")
# average during auction period
auction_year_min <- auction_data[, min(year)]
val <- mean(spot_price_annual[year >= auction_year_min]$spot_price)
write_snt(val, "spot_price_mean_during_auction.tex")

# Brazil and US spot price figure
linewidth <- c(1.2, 0.6)
color <- c("#000000", "#0072B2")
# annual
data <- copy(spot_price_annual)
data[, variable := "Brazil"]
data_tmp <- copy(spot_price_annual_us)
data_tmp[, variable := "US"]
data <- rbind(data, data_tmp)
data[, variable := factor(variable, levels = c("Brazil", "US"))]
fig <- ggplot(data, aes(x=year, y=spot_price, color=variable, linewidth=variable)) +
	geom_line() +
	ylab("Spot Price ($/MWh)") +
	scale_y_continuous(
		  limits = c(0, 300)
		, breaks = seq(0, 300, by = 100)
		, labels = c("0", "100", "200", "300")
	) +
	scale_color_manual(values=color) +
	scale_linewidth_manual(values=linewidth) +
    guides(
		  color = guide_legend(nrow = 1)
		, linewidth = guide_legend(nrow = 1)
    ) +
	theme_bw() +
	theme(legend.title=element_blank()
		, legend.position = c(0.1, 0.98)
		, legend.justification = c(0, 1)
		, legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
		, axis.title.x = element_blank()
	)
ggsave(file.path(fdir, "spotprice_annual.png")
	 , plot = fig, width = 7.5, height = 5, units = "cm")
# monthly
data <- copy(spot_price_monthly[month <= "2022-12-01"])
data[, variable := "Brazil"]
data_tmp <- copy(spot_price_monthly_us)
data_tmp[, variable := "US"]
data <- rbind(data, data_tmp)
data[, variable := factor(variable, levels = c("Brazil", "US"))]
fig <- ggplot(data, aes(x=month, y=spot_price, color=variable, linewidth=variable)) +
	geom_line() +
	ylab("Spot Price ($/MWh)") +
	scale_y_continuous(
		  limits = c(0, 300)
		, breaks = seq(0, 300, by = 100)
		, labels = c("0", "100", "200", "300")
	) +
	scale_color_manual(values=color) +
	scale_linewidth_manual(values=linewidth) +
    guides(
		  color = guide_legend(nrow = 1)
		, linewidth = guide_legend(nrow = 1)
    ) +
	theme_bw() +
	theme(legend.title=element_blank()
		, legend.position = c(0.1, 0.98)
		, legend.justification = c(0, 1)
		, legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
		, axis.title.x = element_blank()
	)
ggsave(file.path(fdir, "spotprice_monthly.png")
	 , plot = fig, width = 7.5, height = 5, units = "cm")

################################################################################
# engineering estimates of cost
################################################################################

# capacity factor (extrapolate from the winners)
auction_data[, capacity_factor := capacity/nameplate_capacity]
auction_data[, year := floor(2010 + year_decimal)]
total_capacity_a <- auction_data[, .(
	  capacity_factor = weighted.mean(capacity_factor, w = capacity)
	, lead_time = weighted.mean(lead_time, w = capacity)
	, total_capacity = sum(total_capacity, na.rm = TRUE)), by = .(year)]
epe_cost_data <- merge(epe_cost_data, total_capacity_a, by = "year")

# interest rate
base_interest_rate_data[, year := as.numeric(format(date, format = "%Y"))]
rate_annual <- base_interest_rate_data[, .(rate = mean(rate)/100), by = .(year)]
epe_cost_data <- merge(epe_cost_data, rate_annual, by = "year")

# borrow money and pay back over T years
epe_cost_data[, rate_factor := T*((delta^(T-1))*(1+rate)^T)*(1-delta*(1+rate))/(1-(delta^T)*((1+rate)^T))]
# CAPEX costs 1 year before OM starts (construction takes 1 year)
epe_cost_data[, cost20yr_mean := (CAPEX_mean*rate_factor + OM_mean*T*delta*delta_tilde_0)*(1/capacity_factor)]
epe_cost_data[, cost20yr_MWh_mean := cost20yr_mean*(10^3)/(24*365.25*T)]
epe_cost_data[, cost20yr_variance := (CAPEX_variance*(rate_factor^2) + OM_variance*((T*delta*delta_tilde_0)^2)
	 + 2*sqrt(CAPEX_variance)*sqrt(OM_variance)*rate_factor*T*delta*delta_tilde_0)*((1/capacity_factor)^2)]
epe_cost_data[, cost20yr_MWh_variance := cost20yr_variance*(((10^3)/(24*365.25*T))^2)]
epe_cost_data[, cost20yr_MWh_sd := sqrt(cost20yr_MWh_variance)]
# lead time discount inclusive
epe_cost_data[, cost20yr_MWh_mean_lead := (delta^lead_time)*cost20yr_MWh_mean]
epe_cost_data[, cost20yr_MWh_sd_lead := (delta^lead_time)*cost20yr_MWh_sd]
cols <- c("cost20yr_MWh_mean_lead","cost20yr_MWh_sd_lead")
# start year
val <- min(epe_cost_data$year)
write_snt(val, "engineering_cost_year_min.tex", 0)
# summary until 2015 (same as sealed-bid auctions)
epe_cost_summary <- data.table(variable = cols
	, value = epe_cost_data[year <= 2015, sapply(.SD, weighted.mean, w = total_capacity), .SDcols = cols])

# distribution of 20-year discounted sum of engineering cost ($/MWh)
val <- epe_cost_summary[variable == "cost20yr_MWh_mean_lead"]$value
write_snt(val, "engineering_cost_mean.tex")
val <- epe_cost_summary[variable == "cost20yr_MWh_sd_lead"]$value
write_snt(val, "engineering_cost_sd.tex")
