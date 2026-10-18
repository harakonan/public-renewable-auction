# Brazil new energy auction structural model estimation
# estimation results

# @author: Konan Hara

# @date: September 14, 2026

# @software: R version 4.5.3
# clear environment
rm(list = ls())

# package loading
library(here)
library(data.table)
library(dplyr)
library(msm)
library(kableExtra)
library(ggplot2)

# path to data
pathtodata <- readLines(file(here("data.txt")), n = 1, warn = FALSE)

# output path
rdir <- here("output", "result")
fdir <- here("output", "figure")
tdir <- here("output", "table")
sdir <- here("output", "snt")

# create target output directories if they do not exist
dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
dir.create(tdir, showWarnings = FALSE, recursive = TRUE)
dir.create(sdir, showWarnings = FALSE, recursive = TRUE)

# load auction data
auction_data <- fread(file.path(pathtodata, "processed", "generation_auction_wind20.csv"))
# remove observations with the possibility of partial bid quantity (partial_bid_quantity)
auction_data <- auction_data[partial_bid_quantity == 0]

# load spot market parameters
spot_result <- fread(file.path(rdir, "spot_result.csv"))

# load specified parameters
params_specified <- fread(file.path(rdir, "params_specified.csv"))

# load results
# first stage
params1_result <- fread(file.path(rdir, "params1_result.csv"))
# second stage
params2_dist_result <- fread(file.path(rdir, "params2_dist_result.csv"))
winprob_grid <- fread(file.path(pathtodata, "intermediate", "winprob_grid.csv"))
gamma_c_data <- fread(file.path(pathtodata, "intermediate", "gamma_c_data.csv"))
gamma_c_dist <- fread(file.path(rdir, "gamma_c_dist.csv"))
c_dist_lp_result <- fread(file.path(rdir, "c_dist_lp_result.csv"))
gamma_c_dist_result <- fread(file.path(rdir, "gamma_c_dist_result.csv"))
gamma_c_winner_result <- fread(file.path(rdir, "gamma_c_winner_result.csv"))

# goodness-of-fit. Use the expected observed bid share, integrating over the
# censored-normal bid-share shock, rather than the clamped latent optimum.
bid_share_covariates <- fread(file.path(
	pathtodata, "intermediate", "bid_share_predicted.csv"
))
bid_share_diagnostics <- fread(file.path(
	pathtodata, "intermediate", "bid_share_gmm_diagnostics.csv"
))[model == "hetero"]
stopifnot(
	nrow(bid_share_covariates) == nrow(bid_share_diagnostics),
	all(bid_share_covariates$auction_id == bid_share_diagnostics$auction_id),
	all(abs(bid_share_covariates$bid_share - bid_share_diagnostics$bid_share) < 1e-10)
)
bid_share_predicted <- bid_share_covariates
bid_share_predicted[, `:=`(
	  bid_share = bid_share_diagnostics$bid_share
	, bid_share_predicted = bid_share_diagnostics$predicted_mean
)]
bid_share_predictive_cdf <- fread(file.path(
	here("notes", "gmm", "step1_gmm_predictive_cdf.csv")
))[model == "hetero" & group_type == "all" & group_value == "all"]

# configuration parameters
# purchase agreement length
T <- 20
# discount factor
delta <- 0.95
# delta tilde at lead time = 0
delta_tilde_0 <- sum(delta^(0:(T - 1))) / T
# sealed-bid auction IDs
auction_data_a <- distinct(auction_data[
	, .(auction_id, descending_clock, number_participant)
])
A_sb <- sort(auction_data_a[descending_clock == 0]$auction_id)
A_sb_nomissing <- sort(auction_data_a[
	descending_clock == 0 & !is.na(number_participant)
]$auction_id)
# number of bootstrap (same as prep.jl)
bs_num <- 200
# number of sensitivity analysis (same as step1_sa.jl)
sa_num <- 18
# grid interval (same as prep.jl)
grid_interval <- 0.2

# function to format
# value
value_format <- function(value, d) {
	paste0("$", sprintf(paste0("%.", d, "f"), value), "$")
}
# estimates with SEs
est_format <- function(coef, se, d) {
	est <- sprintf(paste0("$%.", d, "f$ $(%.", d, "f)$"), coef, se)
	est[is.na(coef)] <- "-"
	est
}
# estimates with central 95% bootstrap percentile intervals
est_ci_format <- function(coef, lower, upper, d) {
	est <- sprintf(paste0("$%.", d, "f$ $[%.", d, "f, %.", d, "f]$"),
		coef, lower, upper)
	est[is.na(coef)] <- "-"
	est
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

################################################################################
# configuration parameters
################################################################################

write_snt(bs_num, "bs_num.tex", 0)
write_snt(grid_interval, "grid_interval.tex")

################################################################################
# mean reverting process of annual spot prices
################################################################################

# parameter table
# row names in table
variable_name <- c(
	  "Intercept, $A$"
	, "AR(1) Coefficient, $\\rho$"
	, "Standard Deviation of Residual, $\\sigma_\\xi$"
	)
table <- data.table(
	  variable = variable_name
	, est = c(
		  est_format(spot_result[variable == "A"]$est, spot_result[variable == "A"]$se, 1)
		, est_format(spot_result[variable == "rho"]$est, spot_result[variable == "rho"]$se, 2)
		, est_format(
			sqrt(spot_result[variable == "var_xi"]$est),
			spot_result[variable == "var_xi"]$se /
				(2*sqrt(spot_result[variable == "var_xi"]$est)),
			1
		)
		)
)
# create latex file
kbl(table
	, format = "latex"
	, col.names = c("Parameter","Estimate")
	, align = "lc"
	, escape = FALSE
	, booktabs = TRUE
	, position = "htb!"
	) %>%
	save_kable(file.path(tdir, "spotestimates.tex"))

# wholesale market SD ($/MWh)
# range
sigma2_r <- spot_result[grep("var_r", variable)]$est
write_snt(min(sqrt(sigma2_r)), "sigma2_r_min.tex")
write_snt(max(sqrt(sigma2_r)), "sigma2_r_max.tex")
# sealed-bid auction capacity-weighted average
capacity_a <- auction_data[descending_clock == 0, .(capacity = sum(capacity)), by = .(auction_id)]
value <- weighted.mean(sqrt(sigma2_r[A_sb]), w = capacity_a$capacity)
sprintf(paste0("wholesale market SD ($/MWh): sealed-bid auction capacity-weighted average, "
	,"$%.", 2, "f$"), value)

# stationary expected wholesale price implied by the spot price mean reverting process ($/MWh)
A <- spot_result[variable == "A"]$est
rho <- spot_result[variable == "rho"]$est
write_snt(A/(1-rho), "stationary_wholesale_price.tex")

################################################################################
# prespecified auction distribution
################################################################################

# parameter table
# digits to the right of the decimal point to be rounded
d <- 2
# procurement capacity distribution
variable_name <- c(
	  "Intercept, $\\beta_{D 0} $"
	, "Operation Start (year), $\\beta_{D 1} $"
	, "Forecasted Marginal Cost (\\$/MWh), $\\beta_{D 2} $"
	, "100 Participants, $\\beta_{D 3} $"
	, "Standard Deviation of Residual, $\\sigma_D $"
	)
table_D <- data.table(
	  variable = variable_name
	, sb = value_format(c(
			  params_specified[variable %in% paste0(paste0("beta_D", 0:2), "_sb")]$value
			, 100*params_specified[variable == "beta_D3_sb"]$value
			, sqrt(params_specified[variable == "sigma2_D_sb"]$value)
		), d)
	, dc = value_format(c(
			  params_specified[variable %in% paste0(paste0("beta_D", 0:2), "_dc")]$value
			, 100*params_specified[variable == "beta_D3_dc"]$value
			, sqrt(params_specified[variable == "sigma2_D_dc"]$value)
		), d)
)
# equilibrium clearing price distribution
variable_name <- c(
	  "Intercept, $\\beta_{pD 0} $"
	, "Log Procurement Capacity (MW), $\\beta_{pD 1} $"
	, "Operation Start (year), $\\beta_{pD 2} $"
	, "Forecasted Marginal Cost (\\$/MWh), $\\beta_{pD 3} $"
	, "100 Participants, $\\beta_{pD 4} $"
	, "Standard Deviation of Residual, $\\sigma_{pD} $"
	)
table_pD <- data.table(
	  variable = variable_name
	, sb = "-"
	, dc = value_format(c(
			  params_specified[variable == "beta_pD0"]$value
			, params_specified[variable == "beta_pD1"]$value
			, params_specified[variable == "beta_pD2"]$value
			, params_specified[variable == "beta_pD3"]$value
			, 100*params_specified[variable == "beta_pD4"]$value
			, sqrt(params_specified[variable == "sigma2_pD"]$value)
		), d)
)
table <- rbind(table_D, table_pD, fill = TRUE)
# create latex file
kbl(table
	, format = "latex"
	, col.names = c("Parameter","Sealed-Bid","Descending-Clock")
	, align = "lcc"
	, escape = FALSE
	, booktabs = TRUE
	, position = "htb!"
	) %>%
	pack_rows(
	  group_label = "A. Log Procurement Capacity (MW) Distribution"
	, start_row = 1
	, end_row = table_D[,.N]
	, indent = FALSE
	) %>%
	pack_rows(
	  group_label = "B. Equilibrium Clearing Price ($/MWh) Distribution"
	, start_row = table_D[,.N]+1
	, end_row = table_D[,.N]+table_pD[,.N]
	, indent = FALSE
	) %>%
	save_kable(file.path(tdir, "prespecdist.tex"))

# expected procurement capacity range (MW)
procurement_capacity_a <- copy(distinct(auction_data[
	,.(auction_id,descending_clock,year_decimal,lead_time,mc_operation,number_participant)]))
beta_D_sb <- params_specified[variable %in% paste0(paste0("beta_D", 0:3), "_sb")]$value
beta_D_dc <- params_specified[variable %in% paste0(paste0("beta_D", 0:3), "_dc")]$value
procurement_capacity_a[descending_clock == 0
	, mu_D := exp(beta_D_sb[1] + 
			beta_D_sb[2]*(year_decimal+lead_time) + 
			beta_D_sb[3]*mc_operation + 
			beta_D_sb[4]*number_participant)
]
procurement_capacity_a[descending_clock == 1
	, mu_D := exp(beta_D_dc[1] + 
			beta_D_dc[2]*(year_decimal+lead_time) + 
			beta_D_dc[3]*mc_operation + 
			beta_D_dc[4]*number_participant)
]
# sealed-bid auctions
val <- min(procurement_capacity_a[descending_clock == 0]$mu_D, na.rm = TRUE)
write_snt(val, "exp_proccapacity_sb_min.tex", 1)
val <- max(procurement_capacity_a[descending_clock == 0]$mu_D, na.rm = TRUE)
write_snt(val, "exp_proccapacity_sb_max.tex", 1)
# descending-clock auctions
val <- min(procurement_capacity_a[descending_clock == 1]$mu_D)
write_snt(val, "exp_proccapacity_dc_min.tex", 1)
val <- max(procurement_capacity_a[descending_clock == 1]$mu_D)
write_snt(val, "exp_proccapacity_dc_max.tex", 1)
# expected clearing price range ($/MWh)
val <- min(params_specified[grep("mu_p_", variable)]$value)
write_snt(val, "exp_clearprice_min.tex")
val <- max(params_specified[grep("mu_p_", variable)]$value)
write_snt(val, "exp_clearprice_max.tex")
# clearing price SD ($/MWh)
val <- sqrt(params_specified[variable == "sigma2_p"]$value)
write_snt(val, "clearprice_sd.tex")

################################################################################
# first stage
################################################################################

# row names in table
variable_name <- c(
	  "Intercept, $\\alpha_{\\gamma 0}$"
	, "Log [1 + Prior Started Wind Capacity (100 MW)], $\\alpha_{\\gamma 1}$"
	, "Intercept, $\\alpha_{r 0}$"
	, "Forecasted Marginal Cost of Electricity (\\$/MWh), $\\alpha_{r 1}$"
	, "Standard Deviation of Bid Share Shock, $\\sigma_\\eta$"
	)
table <- data.table(
	  variable = variable_name
	, est = c(
		  est_format(params1_result[grep("alpha_gamma", variable)]$hetero
			, params1_result[grep("alpha_gamma", variable)]$hetero_se
			, 3)
		, est_format(params1_result[variable == "alpha_r0"]$hetero
			, params1_result[variable == "alpha_r0"]$hetero_se
			, 2)
		, est_format(params1_result[variable == "alpha_r1"]$hetero
			, params1_result[variable == "alpha_r1"]$hetero_se
			, 3)
		# adjustment for log transformation in the estimation of alpha_eta
		, est_format(exp(params1_result[variable == "alpha_eta"]$hetero)
			, deltamethod(~ exp(x1)
				, params1_result[variable == "alpha_eta"]$hetero
				, params1_result[variable == "alpha_eta"]$hetero_se^2)
			, 3)
		)
)
# create latex file
kbl(table
	, format = "latex"
	, col.names = c("Parameter","Estimate")
	, align = "lc"
	, escape = FALSE
	, booktabs = TRUE
	, position = "htb!"
	, linesep = ""
	) %>%
	pack_rows(
	  group_label = "Log Risk Aversion Coefficient, $\\log(\\\\gamma_i)$"
	, start_row = 1
	, end_row = params1_result[grep("alpha_gamma", variable), .N]
	, escape = FALSE
	, bold = FALSE
	) %>%
	pack_rows(
	  group_label = "Discounted Expected Wholesale Price (\\$/MWh), $\\\\mu_r$"
	, start_row = params1_result[grep("alpha_gamma", variable), .N] + 1
	, end_row = params1_result[grep("alpha_gamma", variable), .N] + params1_result[grep("alpha_r", variable), .N]
	, escape = FALSE
	, bold = FALSE
	) %>%
	save_kable(file.path(tdir, "first_stage_estimates.tex"))

# risk aversion coefficient among winners in the cost-estimation sample
gamma <- gamma_c_winner_result[variable == "gamma"]$estimate
write_snt(gamma, "gamma_mean_winner.tex")
# capacity-weighted average risk aversion coefficient among all winners used in
# the first-stage estimation
alpha_gamma <- c(
	  params1_result[variable == "alpha_gamma0"]$hetero
	, params1_result[variable == "alpha_gammaK"]$hetero
)
gamma_all <- exp(
	alpha_gamma[1] +
	alpha_gamma[2] * log1p(auction_data$prior_started_wind_capacity / 100)
)
gamma_all <- weighted.mean(gamma_all, auction_data$capacity)
write_snt(gamma_all, "gamma_mean_winner_all.tex")
# certainty equivalent of 50-50 chance of l-h $/MWh ($/MWh)
l <- 30
h <- 40
write_snt(l, "lottery_low.tex", 0)
write_snt(h, "lottery_high.tex", 0)
write_snt((l+h)/2, "lottery_mean.tex", 0)
# utility function
u <- function(profit) (1-exp(-gamma_all*profit))
# inverse utility function
uinv <- function(utility) (-(1/gamma_all)*(log(1-utility)))
val <- uinv(0.5*u(l)+0.5*u(h))
write_snt(val, "gamma_mean_winner_all_lotteryce.tex")
# prior started wind capacity coefficient
val <- params1_result[variable == "alpha_gammaK"]$hetero
# change in gamma from moving from zero to 100 MW of prior started capacity
val <- exp(val * log(2)) - 1
write_snt(abs(val), "gamma_prior_capacity_100mw_decline_pct.tex", 1, 100, TRUE)

# expected wholesale price belief
alpha_r <- c(
	  params1_result[variable == "alpha_r0"]$hetero
	, params1_result[variable == "alpha_r1"]$hetero
)
mc_a <- unique(auction_data$mc_operation)
val <- min(alpha_r[1] + alpha_r[2]*mc_a)/delta_tilde_0
write_snt(val, "exp_wholesale_price_min.tex")
val <- max(alpha_r[1] + alpha_r[2]*mc_a)/delta_tilde_0
write_snt(val, "exp_wholesale_price_max.tex")
################################################################################
# equilibrium winning probability function
################################################################################

# parameter table
# digits to the right of the decimal point to be rounded
d <- 2
# row names in table
variable_name <- c(
	  "Bid Price (\\$10/MWh)"
	, "Log Capacity (MW)"
	, paste0("Auction ", A_sb)
	, "Standard Deviation of Residual"
		)
# auction fixed effects
auction_fe <- params2_dist_result[grep("beta_C_a", variable)]
auction_fe[, auction_id := as.integer(gsub("beta_C_a", "", variable))]
setnames(auction_fe, c("estimate","se"), c("est_capacity","se_capacity"))
auction_fe[, variable := NULL]
tmp <- params2_dist_result[grep("beta_b_a", variable)]
tmp[, auction_id := as.integer(gsub("beta_b_a", "", variable))]
setnames(tmp, c("estimate","se"), c("est_bidprice","se_bidprice"))
tmp[, variable := NULL]
auction_fe <- merge(auction_fe, tmp, by = "auction_id", all = TRUE)
# prior started wind capacity auction fixed effects
tmp <- params2_dist_result[grep("beta_K_a", variable)]
tmp[, auction_id := as.integer(gsub("beta_K_a", "", variable))]
setnames(tmp, c("estimate", "se"), c("est_prior_capacity", "se_prior_capacity"))
tmp[, variable := NULL]
auction_fe <- merge(auction_fe, tmp, by = "auction_id", all = TRUE)
setorder(auction_fe, auction_id)
first_positive_prior_time <- min(auction_data[
	auction_id %in% A_sb & prior_started_wind_capacity > 0, year_decimal
])
prior_capacity_by_auction <- auction_data[auction_id %in% A_sb
	, .(
		  auction_time = first(year_decimal)
		, positive_support = any(prior_started_wind_capacity > 0)
		, winning_bid_N = .N
	)
	, by = auction_id]
prior_capacity_by_auction[, `:=`(
	  estimable = auction_time >= first_positive_prior_time & positive_support
)]
write_snt(prior_capacity_by_auction[estimable == TRUE, .N]
	, "prior_capacity_estimation_auction_N.tex", 0)
write_snt(prior_capacity_by_auction[estimable == TRUE, sum(winning_bid_N)]
	, "prior_capacity_estimation_bidder_N.tex", 0)
# table
table <- data.table(
		  variable = variable_name
		, est_bidprice = c(
			  "-", "-"
			, est_format(auction_fe$est_bidprice, auction_fe$se_bidprice, d)
		, est_format(sqrt(params2_dist_result[variable == "sigma2_b"]$estimate)
			, deltamethod(~ sqrt(x1)
				, params2_dist_result[variable == "sigma2_b"]$estimate
				, params2_dist_result[variable == "sigma2_b"]$se^2)
			, d)
		)
		, est_capacity = c(
			  est_format(10*params2_dist_result[variable == "beta_C_b"]$estimate
				, 10*params2_dist_result[variable == "beta_C_b"]$se
				, d)
			, "-"
			, est_format(auction_fe$est_capacity, auction_fe$se_capacity, d)
		, est_format(sqrt(params2_dist_result[variable == "sigma2_C"]$estimate)
			, deltamethod(~ sqrt(x1)
				, params2_dist_result[variable == "sigma2_C"]$estimate
				, params2_dist_result[variable == "sigma2_C"]$se^2)
				, d)
			)
		, est_prior_capacity = c(
			  est_format(10*params2_dist_result[variable == "beta_K_b"]$estimate
				, 10*params2_dist_result[variable == "beta_K_b"]$se
				, d)
			, est_format(params2_dist_result[variable == "beta_K_logC"]$estimate
				, params2_dist_result[variable == "beta_K_logC"]$se
				, d)
			, est_format(auction_fe$est_prior_capacity
				, auction_fe$se_prior_capacity
				, d)
			, est_format(sqrt(params2_dist_result[variable == "sigma2_K"]$estimate)
				, deltamethod(~ sqrt(x1)
					, params2_dist_result[variable == "sigma2_K"]$estimate
					, params2_dist_result[variable == "sigma2_K"]$se^2)
				, d)
			)
)

# create latex file
kbl(table
	, format = "latex"
	, col.names = c(
		  "Parameter"
		, "Bid Price (\\$/MWh)"
		, "Log Capacity (MW)"
		, "Log [1+$K$ (100 MW)]"
	)
	, align = "lccc"
	, escape = FALSE
	, booktabs = TRUE
	, position = "htb!"
	, linesep = ""
	) %>%
	pack_rows(
	  group_label = "Auction Fixed Effects"
	, start_row = min(grep("Auction", variable_name))
	, end_row = max(grep("Auction", variable_name))
	, bold = FALSE
	) %>%
	row_spec(nrow(table) - 1, extra_latex_after = "\\addlinespace[0.3em]") %>%
	save_kable(file.path(tdir, "winprobparam.tex"))

# capacity distribution
beta_C_b <- params2_dist_result[variable == "beta_C_b"]$estimate
val <- exp(min(auction_fe$est_capacity + beta_C_b*auction_fe$est_bidprice, na.rm = TRUE))
write_snt(val, "capacity_mean_min.tex", 1)
val <- exp(max(auction_fe$est_capacity + beta_C_b*auction_fe$est_bidprice, na.rm = TRUE))
write_snt(val, "capacity_mean_max.tex", 1)

# winning probability function figure
# superimpose predicted winning probability to data by linear interpolation
for (a in A_sb_nomissing) {
	auction_data[auction_id == a, winprob := approx(
		  winprob_grid[auction_id == a]$bid_price
		, winprob_grid[auction_id == a]$winprob
		, bid_price, rule = 2)$y]
}
winprob_grid[, bid_price_min := min(bid_price), by = .(auction_id)]
winprob_grid[, bid_price_max := max(bid_price), by = .(auction_id)]
winprob_grid[, winprob_min := min(winprob), by = .(auction_id)]
winprob_grid[, winprob_max := max(winprob), by = .(auction_id)]
winprob_grid_tmp <- CJ(auction_id = A_sb_nomissing, bid_price = as.character(seq(20.0,70.0,grid_interval)))
winprob_grid[, bid_price := as.character(bid_price)]
winprob_grid <- merge(winprob_grid_tmp, winprob_grid, by = c("auction_id","bid_price"), all.x = TRUE)
winprob_grid[, bid_price := as.numeric(bid_price)]
winprob_grid[, bid_price_min := min(bid_price_min, na.rm = TRUE), by = .(auction_id)]
winprob_grid[, bid_price_max := max(bid_price_max, na.rm = TRUE), by = .(auction_id)]
winprob_grid[, winprob_min := min(winprob, na.rm = TRUE), by = .(auction_id)]
winprob_grid[, winprob_max := max(winprob, na.rm = TRUE), by = .(auction_id)]
winprob_grid[is.na(winprob) & (bid_price <= bid_price_min), winprob := winprob_max]
winprob_grid[is.na(winprob) & (bid_price >= bid_price_max), winprob := winprob_min]
levels <- c(A_sb_nomissing)
labels <- paste0("Auction ", A_sb_nomissing)
winprob_grid[, auction_label := factor(auction_id, levels=levels, labels=labels)]
auction_data[, auction_label := factor(auction_id, levels=levels, labels=labels)]
fig <- ggplot() +
	geom_line(aes(x=bid_price, y=winprob), winprob_grid) +
	geom_point(aes(x=bid_price, y=winprob), auction_data[auction_id %in% A_sb_nomissing]) +
	xlab("Bid Price ($/MWh)") +
	ylab("Winning Probability") +
	xlim(c(20,NA)) +
	facet_wrap( ~ auction_label, nrow=3) +
	theme_bw() +
	theme(legend.title=element_blank())
ggsave(file.path(fdir, "winprob.png")
	, plot = fig, width = 15, height = 10, units = "cm")

################################################################################
# risk aversion and cost distribution
################################################################################

# gamma_c_data
# gamma_c_dist
# gamma_c_winner_result

# summary statistics of the within-auction joint distribution
# Preserve the paper's capacity-weighted average of within-auction moments.
# For nonlinear statistics, transform each bootstrap draw before constructing
# percentile intervals.
auction_capacity_weights <- unique(auction_data[
	descending_clock == 0 & !is.na(number_participant),
	.(auction_id, total_capacity)
])

weighted_mean_complete <- function(x, w) {
	keep <- is.finite(x) & is.finite(w)
	if (!any(keep)) return(NA_real_)
	weighted.mean(x[keep], w[keep])
}

within_auction_dist_summary <- function(dt) {
	dt <- merge(copy(dt), auction_capacity_weights, by = "auction_id")
	c(
		  gamma_mean = weighted_mean_complete(dt$gamma_mean, dt$total_capacity)
		, gamma_sd = sqrt(weighted_mean_complete(dt$gamma_var, dt$total_capacity))
		, cost_mean = weighted_mean_complete(dt$c_mean, dt$total_capacity)
		, cost_sd = sqrt(weighted_mean_complete(dt$c_var, dt$total_capacity))
		, gamma_c_cor = weighted_mean_complete(dt$gamma_c_cor, dt$total_capacity)
	)
}

dist_summary <- within_auction_dist_summary(gamma_c_dist)
dist_summary_bs <- t(vapply(seq_len(bs_num), function(index_bs) {
	within_auction_dist_summary(fread(file.path(
		pathtodata, "intermediate", "params2_bs", "gamma_c_dist",
		paste0("bs", index_bs, ".csv")
	)))
}, numeric(length(dist_summary))))
colnames(dist_summary_bs) <- names(dist_summary)
dist_summary_ci <- apply(
	dist_summary_bs, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE
)

# row names in table
variable_name <- c(
	  "Mean [95\\% CI]"
	, "Standard Deviation [95\\% CI]"
	, "Correlation [95\\% CI]"
	)
table <- data.table(
	  variable = variable_name
	, est_gamma = c(
		  est_ci_format(dist_summary["gamma_mean"]
			, dist_summary_ci[1, "gamma_mean"], dist_summary_ci[2, "gamma_mean"], 2)
		, est_ci_format(dist_summary["gamma_sd"]
			, dist_summary_ci[1, "gamma_sd"], dist_summary_ci[2, "gamma_sd"], 2)
		, paste0("\\multicolumn{2}{c}{"
			, est_ci_format(dist_summary["gamma_c_cor"]
				, dist_summary_ci[1, "gamma_c_cor"]
				, dist_summary_ci[2, "gamma_c_cor"], 3)
			, "}")
		)
	, est_cost = c(
		  est_ci_format(dist_summary["cost_mean"]
			, dist_summary_ci[1, "cost_mean"], dist_summary_ci[2, "cost_mean"], 2)
		, est_ci_format(dist_summary["cost_sd"]
			, dist_summary_ci[1, "cost_sd"], dist_summary_ci[2, "cost_sd"], 2)
		, ""
		)
)
# create latex file
latex_table <- kbl(table
	, format = "latex"
	, col.names = c("Statistic","Risk Aversion","Cost (\\$/MWh)")
	, align = "lc"
	, escape = FALSE
	, booktabs = TRUE
	, position = "htb!"
	, linesep = c("", "\\addlinespace[0.3em]")
	)
latex_table <- gsub("(\\\\multicolumn\\{2\\}\\{c\\}\\{.*?\\}) &", "\\1", latex_table)
save_kable(latex_table, file.path(tdir, "gamma_cost_dist.tex"))

# gamma distribution
gamma_mean <- dist_summary["gamma_mean"]
gamma_sd <- dist_summary["gamma_sd"]
write_snt(gamma_sd, "gamma_sd_participant.tex")
# above 1 SD/below 1 SD ratio
val <- (gamma_mean+gamma_sd)/(gamma_mean-gamma_sd)
write_snt(val, "gamma_1sd_ratio.tex", 1)

# cost distribution
# cost difference between the participants and winners ($/MWh)
cost_diff <- gamma_c_dist_result[variable == "c_mean"]$estimate - gamma_c_winner_result[variable == "cost"]$estimate
write_snt(cost_diff, "cost_diff_participant_winner.tex")
# capacity-weighted average of within-auction participant-winner cost gaps,
# standardized by each auction's participant cost SD
winner_cost_a <- gamma_c_data[!is.na(cost),
	.(winner_cost_mean = weighted.mean(cost, capacity)), by = auction_id]
cost_gap_a <- merge(
	gamma_c_dist[, .(auction_id, participant_cost_mean = c_mean, participant_cost_sd = sqrt(c_var))],
	winner_cost_a, by = "auction_id"
)
cost_gap_a <- merge(cost_gap_a, auction_capacity_weights, by = "auction_id")
cost_gap_a[, standardized_cost_gap :=
	(participant_cost_mean - winner_cost_mean) / participant_cost_sd]
val <- weighted.mean(cost_gap_a$standardized_cost_gap, cost_gap_a$total_capacity)
write_snt(val, "cost_diff_within_auction_sd_ratio.tex", 1)
# average cost from 2010 to 2015
total_capacity_a <- copy(distinct(
	auction_data[descending_clock == 0 & !is.na(number_participant)
		, .(auction_id, year, total_capacity)]))
data <- merge(gamma_c_dist, total_capacity_a, by = "auction_id")
val <- data[year >= 2010 & year <= 2015, sum(c_mean*total_capacity)/sum(total_capacity)]
write_snt(val, "cost_mean_participant_post2010.tex")

# cost linear projection
# row order, names, and digits in table
lp_table_spec <- data.table(
	  variable = c("intercept", "t", "t2", "t3", "l", "t*l", "t2*l", "t3*l")
	, variable_name = c(
		  "Intercept"
		, "Auction Date (year)"
		, "Auction Date Square"
		, "Auction Date Cubic"
		, "Lead Time (year)"
		, paste0("Lead Time $", "\\", "times$ Auction Date")
		, paste0("Lead Time $", "\\", "times$ Auction Date Square")
		, paste0("Lead Time $", "\\", "times$ Auction Date Cubic")
		)
	, digits = c(2, 2, 2, 3, 2, 2, 2, 3)
	)
lp_table_result <- c_dist_lp_result[match(lp_table_spec$variable, c_dist_lp_result$variable)]
table <- data.table(
	  variable = lp_table_spec$variable_name
	, est_main = mapply(est_format
		, lp_table_result$main_estimate
		, lp_table_result$main_se
		, lp_table_spec$digits
		)
	, est_undiscounted = mapply(est_format
		, lp_table_result$undiscounted_estimate
		, lp_table_result$undiscounted_se
		, lp_table_spec$digits
		)
)
# create latex file
kbl(table
	, format = "latex"
	, col.names = c("Variable","Main","Undiscounted")
	, align = "lcc"
	, escape = FALSE
	, booktabs = TRUE
	, position = "htb!"
	, linesep = ""
	) %>%
	save_kable(file.path(tdir, "cost_dist_lp.tex"))

# capacity-weighted average risk premium, cost, markup
# risk premium
val <- gamma_c_winner_result[variable == "rp"]$estimate
write_snt(val, "rp_mean_winner.tex")
val <- gamma_c_winner_result[variable == "rp"]$ci_lower
write_snt(val, "rp_mean_winner_ci_lower.tex")
val <- gamma_c_winner_result[variable == "rp"]$ci_upper
write_snt(val, "rp_mean_winner_ci_upper.tex")
# cost
val <- gamma_c_winner_result[variable == "cost"]$estimate
write_snt(val, "cost_mean_winner.tex")
val <- gamma_c_winner_result[variable == "cost"]$ci_lower
write_snt(val, "cost_mean_winner_ci_lower.tex")
val <- gamma_c_winner_result[variable == "cost"]$ci_upper
write_snt(val, "cost_mean_winner_ci_upper.tex")
# markup
val <- gamma_c_winner_result[variable == "markup"]$estimate
write_snt(val, "markup_mean_winner.tex")
val <- gamma_c_winner_result[variable == "markup"]$ci_lower
write_snt(val, "markup_mean_winner_ci_lower.tex")
val <- gamma_c_winner_result[variable == "markup"]$ci_upper
write_snt(val, "markup_mean_winner_ci_upper.tex")
# risk premium/cost ratio
val <- gamma_c_winner_result[variable == "rp"]$estimate/gamma_c_winner_result[variable == "cost"]$estimate
write_snt(val, "rp_cost_ratio_mean_winner.tex", 1, 100, TRUE)

# overall impact
# capacity (GW)
val <- sum(gamma_c_data[!is.na(cost)]$capacity)/1e+3
write_snt(val, "total_capacity_sb_nomissing.tex", 1)
# compare to the overall capacity in Brazil as of December 2020 (Tolmasquim et al., 2021)
write_snt(val/177.16, "total_capacity_pct_sb_nomissing.tex", 1, 100, TRUE)
# risk premium (billion $)
val <- gamma_c_data[, sum(rp*capacity*24*365.25*T, na.rm = TRUE)/1e+9]
write_snt(val, "total_rp_sb_nomissing.tex", 1)
# cost (billion $)
val <- gamma_c_data[, sum(cost*capacity*24*365.25*T, na.rm = TRUE)/1e+9]
write_snt(val, "total_cost_sb_nomissing.tex", 1)

# Wholesale revenue associated with the winners included in the cost
# estimation.  The mean and variance parameters are for the average discounted
# wholesale price, so multiplying them by 20-year production gives discounted
# aggregate revenue.  Construct the cross-auction covariance matrix by aligning
# the annual innovations in the prespecified wholesale-price AR(1) process.
wholesale_revenue_a <- gamma_c_data[!is.na(cost), .(
	capacity = sum(capacity)
), by = auction_id]
wholesale_revenue_a <- merge(
	wholesale_revenue_a,
	unique(auction_data[, .(auction_id, year, lead_time, mc_operation)]),
	by = "auction_id"
)
wholesale_revenue_a[, `:=`(
	expected_wholesale_price = alpha_r[1] + alpha_r[2]*mc_operation,
	production = capacity*24*365.25*T,
	delta_tilde = delta^lead_time*delta_tilde_0
)]
val <- wholesale_revenue_a[, sum(expected_wholesale_price*production)]/1e+9
write_snt(val, "total_wholesale_revenue_mean_sb_nomissing.tex", 1)

# Innovation loadings for an integer lead time.  These are the loadings in
# func_spot_variance() in step1_spot.jl, divided by T because the wholesale
# price object is the average discounted price over the contract horizon.
rho_spot <- spot_result[variable == "rho"]$est
sigma2_xi_spot <- spot_result[variable == "var_xi"]$est
spot_innovation_loading <- function(lead_time_integer) {
	lead_time_integer <- as.integer(lead_time_integer)
	period <- seq_len(lead_time_integer + T - 1L)
	loading <- ifelse(
		period <= lead_time_integer,
		delta^lead_time_integer * rho_spot^(lead_time_integer - period) *
			(1 - delta^T*rho_spot^T)/(1 - delta*rho_spot),
		delta^period *
			(1 - delta^(lead_time_integer + T - period) *
				rho_spot^(lead_time_integer + T - period))/
			(1 - delta*rho_spot)
	)
	loading/T
}

# The prespecified marginal variances linearly interpolate the variances for
# floor and ceiling integer lead times.  Preserve that interpolation on the
# covariance diagonal and use the correspondingly interpolated common-shock
# loadings off the diagonal.
year_min <- min(wholesale_revenue_a$year)
innovation_N <- wholesale_revenue_a[, max(
	year - year_min + ceiling(lead_time) + T - 1L
)]
innovation_loading <- matrix(
	0,
	nrow = nrow(wholesale_revenue_a),
	ncol = innovation_N
)
wholesale_variance <- numeric(nrow(wholesale_revenue_a))
for (row in seq_len(nrow(wholesale_revenue_a))) {
	lead_time <- wholesale_revenue_a$lead_time[row]
	lead_time_floor <- floor(lead_time)
	lead_time_ceiling <- ceiling(lead_time)
	weight_ceiling <- lead_time - lead_time_floor
	weight <- c(1 - weight_ceiling, weight_ceiling)
	loading <- list(
		spot_innovation_loading(lead_time_floor),
		spot_innovation_loading(lead_time_ceiling)
	)
	year_offset <- wholesale_revenue_a$year[row] - year_min
	for (index in seq_along(loading)) {
		columns <- year_offset + seq_along(loading[[index]])
		innovation_loading[row, columns] <-
			innovation_loading[row, columns] + weight[index]*loading[[index]]
	}
	wholesale_variance[row] <- sigma2_xi_spot*sum(
		weight*c(sum(loading[[1]]^2), sum(loading[[2]]^2))
	)
}
wholesale_covariance <- sigma2_xi_spot*tcrossprod(innovation_loading)
diag(wholesale_covariance) <- wholesale_variance
stopifnot(isTRUE(all.equal(
	diag(wholesale_covariance),
	sigma2_r[wholesale_revenue_a$auction_id]
)))
val <- sqrt(drop(
	t(wholesale_revenue_a$production) %*%
		wholesale_covariance %*%
		wholesale_revenue_a$production
))/1e+9
write_snt(val, "total_wholesale_revenue_sd_sb_nomissing.tex", 1)

# Monthly residential-bill equivalent of aggregate policymaker net
# expenditure. First convert each auction's discounted 20-year totals to an
# undiscounted annual average. For the full purchase agreement, use the
# AR(1)-consistent covariance matrix above to preserve correlated wholesale
# price shocks across auctions. A uniform per-MWh surcharge across all end
# users assigns residential customers the residential share of consumption.
bill_cost_risk_a <- gamma_c_data[!is.na(cost), .(
	discounted_cost = sum(cost*capacity*24*365.25*T),
	discounted_risk_premium = sum(rp*capacity*24*365.25*T)
), by = auction_id]
bill_a <- merge(
	wholesale_revenue_a,
	bill_cost_risk_a,
	by = "auction_id"
)
stopifnot(all(bill_a$auction_id == wholesale_revenue_a$auction_id))
bill_a[, discounted_revenue := expected_wholesale_price*production]

fullpa_annual_net_expenditure <- bill_a[, sum(
	(discounted_cost - discounted_revenue)/(delta_tilde*T)
)]
subsidy_annual_net_expenditure <- bill_a[, sum(
	(discounted_cost + discounted_risk_premium - discounted_revenue)/
		(delta_tilde*T)
)]
annual_revenue_weight <- bill_a$production/(bill_a$delta_tilde*T)
fullpa_annual_net_expenditure_sd <- sqrt(drop(
	t(annual_revenue_weight) %*%
		wholesale_covariance %*%
		annual_revenue_weight
))

# EPE Statistical Yearbook of Electric Energy 2021, Tables 2.15, 3.3, and
# 3.34: 2019 residential tariff, consumption, and customer count. Monetary
# values use the paper's fixed R$5/US$ conversion.
residential_rate_brl_mwh <- 551.16
residential_consumption_gwh <- 142781
total_consumption_gwh <- 482226
residential_customer_N <- 73380396
brl_per_usd <- 5
residential_monthly_bill <-
	(residential_consumption_gwh*1e3/residential_customer_N)*
	(residential_rate_brl_mwh/brl_per_usd)/12
monthly_bill_multiplier <-
	(residential_consumption_gwh/total_consumption_gwh)/
	(residential_customer_N*12)
subsidy_monthly_surcharge <-
	subsidy_annual_net_expenditure*monthly_bill_multiplier
fullpa_monthly_bill_change <-
	fullpa_annual_net_expenditure*monthly_bill_multiplier
fullpa_monthly_bill_sd <-
	fullpa_annual_net_expenditure_sd*monthly_bill_multiplier
stopifnot(
	subsidy_monthly_surcharge > 0,
	fullpa_monthly_bill_change < 0,
	fullpa_monthly_bill_sd > 0
)
write_snt(
	residential_monthly_bill,
	"residential_monthly_bill_2019.tex", 2
)
write_snt(
	subsidy_monthly_surcharge,
	"subsidy_monthly_surcharge_cent.tex", 0, 100
)
write_snt(
	-fullpa_monthly_bill_change,
	"fullpa_monthly_bill_reduction_cent.tex", 0, 100
)
write_snt(
	fullpa_monthly_bill_sd,
	"fullpa_monthly_bill_sd_cent.tex", 0, 100
)
write_snt(
	fullpa_monthly_bill_sd/residential_monthly_bill,
	"fullpa_monthly_bill_sd_pct.tex", 1, 100, TRUE
)

# Additional wind capacity whose 20-year project cost could be covered by the
# aggregate risk-premium saving.  This is the current counterpart of the
# previous counterfactual calculation using loser costs: within each auction,
# value additional projects at the estimated participant mean cost.
risk_premium_a <- gamma_c_data[, .(
	risk_premium_total = sum(rp*capacity*24*365.25*T, na.rm = TRUE)
), by = auction_id]
additional_capacity_a <- merge(
	risk_premium_a,
	gamma_c_dist[, .(auction_id, participant_cost_mean = c_mean)],
	by = "auction_id"
)
val <- additional_capacity_a[, sum(
	risk_premium_total / (participant_cost_mean*24*365.25*T)
)] / 1e+3
write_snt(val, "additional_capacity_rp_sb_nomissing.tex", 1)

# cost distribution linear projection

# figure setting
t_grid <- seq(0, 5, by = 0.01)
year_grid <- 2010 + t_grid
design_lp <- cbind(1, t_grid, t_grid^2, t_grid^3)
levels_lp <- c("main", "undiscounted")
labels_lp <- c("Main", "Undiscounted")
color_lp <- c("Main" = "#000000", "Undiscounted" = "#0072B2")
linewidth_lp <- c("Main" = 0.6, "Undiscounted" = 1.2)
c_dist_lp_bs <- lapply(seq_len(bs_num), function(index_bs) {
	fread(file.path(
		pathtodata, "intermediate", "params2_bs", "c_dist_lp",
		paste0("bs", index_bs, ".csv")
	))
})

lp_path <- function(dt, column, terms) {
	coefficients <- dt[match(terms, variable), get(column)]
	stopifnot(length(coefficients) == length(terms), all(is.finite(coefficients)))
	as.vector(design_lp %*% coefficients)
}

lp_bootstrap_ci <- function(column, terms) {
	paths <- vapply(c_dist_lp_bs, function(dt) {
		lp_path(dt, column, terms)
	}, numeric(nrow(design_lp)))
	t(apply(paths, 1, quantile, probs = c(0.025, 0.975), names = FALSE))
}

# cost trend by auction date, setting lead time = 0
# figure
trend_terms <- c("intercept", "t", "t2", "t3")
trend_ci <- lp_bootstrap_ci("main", trend_terms)
plot_dt <- data.table(
	  year = year_grid
	, value = lp_path(c_dist_lp_result, "main_estimate", trend_terms)
	, ci_lower = trend_ci[, 1]
	, ci_upper = trend_ci[, 2]
)
fig <- ggplot(plot_dt, aes(x = year, y = value)) +
	geom_line(aes(y = ci_lower), color = color_lp[1], linewidth = linewidth_lp[1]
		, linetype = "dashed") +
	geom_line(aes(y = ci_upper), color = color_lp[1], linewidth = linewidth_lp[1]
		, linetype = "dashed") +
	geom_line(color = color_lp[1], linewidth = linewidth_lp[1]) +
	scale_x_continuous(breaks = 2010:2015) +
	labs(x = NULL, y = "$/MWh") +
	theme_bw() +
	theme(
		  legend.title = element_blank()
	)
ggsave(file.path(fdir, "c_dist_lp_t_l0.png")
	, plot = fig, width = 7.5, height = 5, units = "cm")

# coefficient on lead time by auction date
lead_terms <- c("l", "t*l", "t2*l", "t3*l")
main_lead_ci <- lp_bootstrap_ci("main", lead_terms)
undiscounted_lead_ci <- lp_bootstrap_ci("undiscounted", lead_terms)
plot_dt <- rbind(
	data.table(
		  year = year_grid
		, specification = "main"
		, value = lp_path(c_dist_lp_result, "main_estimate", lead_terms)
		, ci_lower = main_lead_ci[, 1]
		, ci_upper = main_lead_ci[, 2]
	)
	, data.table(
		  year = year_grid
		, specification = "undiscounted"
		, value = lp_path(c_dist_lp_result, "undiscounted_estimate", lead_terms)
		, ci_lower = undiscounted_lead_ci[, 1]
		, ci_upper = undiscounted_lead_ci[, 2]
	)
)
plot_dt[, specification := factor(specification, levels = levels_lp, labels = labels_lp)]
fig <- ggplot(plot_dt, aes(x = year, color = specification, linewidth = specification)) +
	geom_line(aes(y = ci_lower), linetype = "dashed", show.legend = FALSE) +
	geom_line(aes(y = ci_upper), linetype = "dashed", show.legend = FALSE) +
	geom_line(aes(y = value)) +
	scale_color_manual(values = color_lp) +
	scale_linewidth_manual(values = linewidth_lp) +
	scale_x_continuous(breaks = 2010:2015) +
	labs(x = NULL, y = "$/MWh per year") +
	theme_bw() +
	theme(
		  legend.title = element_blank()
		, legend.position = c(0.5, 0.02)
		, legend.justification = c(0.5, 0)
		, legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
	)
ggsave(file.path(fdir, "c_dist_lp_l_by_t.png")
	, plot = fig, width = 7.5, height = 5, units = "cm")

################################################################################
# goodness-of-fit
################################################################################

# Add the auction-level covariates used in the first-stage model.
gof_auction_covariates <- unique(auction_data[
	, .(auction_id, descending_clock, mc_operation)
])
bid_share_predicted <- merge(
	bid_share_predicted, gof_auction_covariates, by = "auction_id", all.x = TRUE
)

# Ensure levels and labels are set
levels <- c("bid_share", "bid_share_predicted")
labels <- c("Observed", "Predicted")

# Define a consistent color palette and line width for the plots
color <- c("Observed" = "#000000", "Predicted" = "#0072B2")
linewidth <- c("Observed" = 0.6, "Predicted" = 1.2)

# cumulative distribution function
plot_dt <- melt(
	bid_share_predictive_cdf,
	id.vars = "bid_share",
	measure.vars = c("observed_cdf", "predicted_cdf"),
	variable.name = "variable", value.name = "cdf"
)
plot_dt[, variable := factor(
	variable,
	levels = c("observed_cdf", "predicted_cdf"),
	labels = labels
)]
fig <- ggplot(plot_dt, aes(
	x = bid_share, y = cdf, color = variable, linewidth = variable
)) +
	geom_step() +
	scale_color_manual(values = color) +
	scale_linewidth_manual(values = linewidth) +
	labs(x = "Bid Share", y = "CDF") +
	xlim(0,1.01) +
	theme_bw() +
	theme(
		  legend.title = element_blank()
		, legend.position = c(0.02, 0.98) 
		, legend.justification = c(0, 1)
		, legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
	)
ggsave(file.path(fdir, "bidshare_cdf_gof.png")
	, plot = fig, width = 7.5, height = 5, units = "cm")

# ------------------------------------------------------------------------------
# predicted vs observed bid shares by first-stage covariates
# ------------------------------------------------------------------------------
# quantile-bin assignment and summaries
quantile_bin_id <- function(x, bins) {
	breaks <- unique(quantile(x, probs = seq(0, 1, length.out = bins + 1), na.rm = TRUE))
	if (length(breaks) < 2) return(rep(1L, length(x)))
	as.integer(cut(x, breaks = breaks, include.lowest = TRUE))
}
plot_group_dt <- function(group, x, x_name) {
	x_dt <- bid_share_predicted[!is.na(get(group))
		, .(mean_x = mean(get(x), na.rm = TRUE)), by = group]
	long_dt <- melt(
		bid_share_predicted, id.vars = group, measure.vars = levels
	)
	long_dt[, variable := factor(variable, levels = levels, labels = labels)]
	share_dt <- long_dt[!is.na(get(group))
		, .(mean_share = mean(value, na.rm = TRUE)), by = c(group, "variable")]
	setnames(x_dt, "mean_x", x_name)
	merge(share_dt, x_dt, by = group)
}

bid_share_predicted[, price_bin := quantile_bin_id(bid_price, 10)]
mc_auction_bins <- unique(bid_share_predicted[, .(auction_id, mc_operation)])
setorder(mc_auction_bins, mc_operation, auction_id)
mc_auction_bins[, mc_bin := pmin(5L, ceiling(5 * seq_len(.N) / .N))]
bid_share_predicted <- merge(
	bid_share_predicted, mc_auction_bins[, .(auction_id, mc_bin)],
	by = "auction_id", all.x = TRUE
)

# Put the mass at zero in its own group and split positive prior capacity into
# quartiles. This avoids treating the zero mass as several nominal deciles.
positive_prior_breaks <- unique(quantile(
	bid_share_predicted[prior_started_wind_capacity > 0, prior_started_wind_capacity],
	probs = seq(0, 1, length.out = 5), na.rm = TRUE
))
bid_share_predicted[, prior_capacity_bin := 1L]
bid_share_predicted[prior_started_wind_capacity > 0,
	prior_capacity_bin := 1L + as.integer(cut(
		prior_started_wind_capacity,
		breaks = positive_prior_breaks,
		include.lowest = TRUE
	))]

# by prior started wind capacity: zero and positive-capacity quartiles
prior_capacity_dt <- plot_group_dt(
	"prior_capacity_bin", "prior_started_wind_capacity", "mean_prior_capacity")
fig <- ggplot(prior_capacity_dt, aes(
	x = mean_prior_capacity, y = mean_share, color = variable,
	shape = variable, group = variable
)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2.5) +
	scale_color_manual(values = color) +
	labs(x = "Prior Started Wind Capacity (MW)", y = "Bid Share") +
	ylim(0,1) +
	theme_bw() +
	theme(
		  legend.title = element_blank()
		, legend.position = c(0.98, 0.02)
		, legend.justification = c(1, 0)
		, legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
		)
ggsave(file.path(fdir, "bidshare_prior_capacity_gof.png")
	, plot = fig, width = 8, height = 5, units = "cm")

# by price deciles
price_dt <- plot_group_dt("price_bin", "bid_price", "mean_bid_price")
fig <- ggplot(price_dt, aes(x = mean_bid_price, y = mean_share, color = variable, shape = variable, group = variable)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2.5) +
	scale_color_manual(values = color) +
	labs(x = "Purchase Agreement Price ($/MWh)", y = "Bid Share") +
	ylim(0,1) +
	theme_bw() +
	theme(
		  legend.title = element_blank()
		, legend.position = c(0.98, 0.02)
		, legend.justification = c(1, 0)
		, legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
	)
ggsave(file.path(fdir, "bidshare_price_gof.png")
	, plot = fig, width = 8, height = 5, units = "cm")

# by forecasted marginal cost quintiles
mc_dt <- plot_group_dt("mc_bin", "mc_operation", "mean_mc_operation")
fig <- ggplot(mc_dt, aes(
	x = mean_mc_operation, y = mean_share, color = variable,
	shape = variable, group = variable
)) +
	geom_line(linewidth = 0.8) +
	geom_point(size = 2.5) +
	scale_color_manual(values = color) +
	labs(x = "Forecasted Marginal Cost ($/MWh)", y = "Bid Share") +
	ylim(0, 1) +
	theme_bw() +
	theme(
		  legend.title = element_blank()
		, legend.position = c(0.98, 0.02)
		, legend.justification = c(1, 0)
		, legend.background = element_rect(fill = alpha("white", 0.7), color = NA)
	)
ggsave(file.path(fdir, "bidshare_mc_gof.png")
	, plot = fig, width = 8, height = 5, units = "cm")

################################################################################
# sensitivity analysis
################################################################################

# digits to the right of the decimal point to be rounded
d <- 2
# sensitivity-analysis results
sa_result <- rbindlist(lapply(1:sa_num, function(index_sa) {
	params1_sa <- fread(file.path(pathtodata, "intermediate", "params1_sa"
		, paste0("sa", index_sa, ".csv")))
	gamma_c_sa <- fread(file.path(pathtodata, "intermediate", "params2_sa", "gamma_c_data"
		, paste0("sa", index_sa, ".csv")))
	sa_info <- params1_sa[1, .(sa_type, sa_value)]
	sa_outcome <- gamma_c_sa[, lapply(.SD, weighted.mean, w = capacity, na.rm = TRUE)
		, .SDcols = c("rp", "cost", "markup")]
	cbind(data.table(index_sa = index_sa), sa_info, sa_outcome)
}))

format_sa_value <- function(sa_type_arg, sa_value_arg, variable, d) {
	if (is.na(sa_value_arg)) {
		return("")
	}
	value <- sa_result[sa_type == sa_type_arg & abs(sa_value - sa_value_arg) < 1e-8, get(variable)]
	if (length(value) != 1) {
		stop(paste0("Expected one sensitivity result for ", sa_type_arg, " = ", sa_value_arg
			, " and variable ", variable, ". Found ", length(value), "."))
	}
	return(sprintf(paste0("$%.", d, "f$"), value))
}
format_sa_row <- function(sa_type, sa_values, variable, d) {
	sapply(sa_values, function(x) format_sa_value(sa_type, x, variable, d))
}
save_kable_with_full_midrule <- function(kable_input, file, ...) {
	save_kable(kable_input, file, ...)
	latex <- readLines(file, warn = FALSE)
	latex <- gsub(
		  paste0("\\", "cmidrule(l{3pt}r{3pt}){2-6}")
		, paste0("\\", "midrule")
		, latex
		, fixed = TRUE
	)
	writeLines(latex, file, useBytes = TRUE)
}

# column values
change_sd <- seq(0.5, 1.5, 0.25)
change_mean <- seq(-2, 2, 1)
change_delta <- c(NA, 0.90, 0.95, 0.99, NA)
change_sd_name <- c("$-50$\\%", "$-25$\\%", "Main", "$+25$\\%", "$+50$\\%")
change_mean_name <- c("$-2$SD", "$-1$SD", "Main", "$+1$SD", "$+2$SD")
change_delta_name <- c("", "$0.90$", "Main", "$0.99$", "")

# sensitivity-analysis snippets
sa_result[, rp_cost_ratio := rp/cost]
write_snt(min(change_sd), "sa_sd_multiplier_min_pct.tex", 0, 100, TRUE)
write_snt(max(change_sd), "sa_sd_multiplier_max_pct.tex", 0, 100, TRUE)
write_snt(sqrt(spot_result[variable == "var_xi"]$est), "sa_sigma_xi_main.tex", 2)
write_snt(delta, "sa_delta_main.tex", 2)

# overall cost range across Panels A--C
sa_types_panels_abc <- c("sigma2_xi", "mu_p", "sigma2_p")
write_snt(sa_result[sa_type %in% sa_types_panels_abc, min(rp)]
	, "sa_rp_panels_abc_min.tex", d)
write_snt(sa_result[sa_type %in% sa_types_panels_abc, max(rp)]
	, "sa_rp_panels_abc_max.tex", d)
write_snt(sa_result[sa_type %in% sa_types_panels_abc, min(cost)]
	, "sa_cost_panels_abc_min.tex", d)
write_snt(sa_result[sa_type %in% sa_types_panels_abc, max(cost)]
	, "sa_cost_panels_abc_max.tex", d)

write_snt(sa_result[sa_type == "delta", min(rp_cost_ratio)], "sa_rp_cost_ratio_delta_min.tex", 1, 100, TRUE)
write_snt(sa_result[sa_type == "delta", max(rp_cost_ratio)], "sa_rp_cost_ratio_delta_max.tex", 1, 100, TRUE)

# risk premium, cost, and markup table
# row names in table
variable_name <- c(
	  "A. Wholesale Price Innovation SD"
	, "\\hspace{1em}Risk Premium (\\$/MWh)"
	, "\\hspace{1em}Cost (\\$/MWh)"
	, "\\hspace{1em}Markup (\\$/MWh)"
	, "B. Clearing Price Mean, $\\mu_{p^\\ast} $"
	, "\\hspace{1em}Risk Premium (\\$/MWh)"
	, "\\hspace{1em}Cost (\\$/MWh)"
	, "\\hspace{1em}Markup (\\$/MWh)"
	, "C. Clearing Price SD, $\\sigma_{p^\\ast} $"
	, "\\hspace{1em}Risk Premium (\\$/MWh)"
	, "\\hspace{1em}Cost (\\$/MWh)"
	, "\\hspace{1em}Markup (\\$/MWh)"
	, "D. Discount Factor, $\\delta $"
	, "\\hspace{1em}Risk Premium (\\$/MWh)"
	, "\\hspace{1em}Cost (\\$/MWh)"
	, "\\hspace{1em}Markup (\\$/MWh)"
	)
cols <- list()
for (k in 1:length(change_sd)) {
	cols[[k]] <- c(
		  change_sd_name[k]
		, format_sa_row("sigma2_xi", change_sd, "rp", d)[k]
		, format_sa_row("sigma2_xi", change_sd, "cost", d)[k]
		, format_sa_row("sigma2_xi", change_sd, "markup", d)[k]
		, change_mean_name[k]
		, format_sa_row("mu_p", change_mean, "rp", d)[k]
		, format_sa_row("mu_p", change_mean, "cost", d)[k]
		, format_sa_row("mu_p", change_mean, "markup", d)[k]
		, change_sd_name[k]
		, format_sa_row("sigma2_p", change_sd, "rp", d)[k]
		, format_sa_row("sigma2_p", change_sd, "cost", d)[k]
		, format_sa_row("sigma2_p", change_sd, "markup", d)[k]
		, change_delta_name[k]
		, format_sa_row("delta", change_delta, "rp", d)[k]
		, format_sa_row("delta", change_delta, "cost", d)[k]
		, format_sa_row("delta", change_delta, "markup", d)[k]
		)
}
# table
table <- data.table(variable = variable_name)
table[, (paste0("col",1:length(cols))) := cols]
# create latex file
kbl(table
	, format = "latex"
	, col.names = NULL
	, align = "lccccc"
	, escape = FALSE
	, booktabs = TRUE
	, position = "htb!"
	, linesep = rep(c("\\addlinespace[0.1em]", "", "", "\\midrule"), 4)
	) %>%
	add_header_above(
		c("", "Change from the Main Specification"=5)
	) %>%
	save_kable_with_full_midrule(file.path(tdir, "sa_rp_cost_markup.tex"))
