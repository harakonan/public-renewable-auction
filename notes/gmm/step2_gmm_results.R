# Step 2 distribution summaries, trends, and discount-factor normalization
#
# @author: Konan Hara
# @date: September 13, 2026
# @software: R version 4.5

setwd("~/Workspace/project/renewable-auction/")
rm(list = ls())

library(data.table)
library(ggplot2)

pathtodata <- readLines("data.txt", warn = FALSE)[1]
rdir <- file.path("notes", "gmm")
bs_num <- 200

read_bootstrap <- function(directory) {
	lapply(seq_len(bs_num), function(b) {
		fread(file.path(
			pathtodata, "intermediate", "params2_bs", directory,
			paste0("bs", b, ".csv")
		))
	})
}

percentile_summary <- function(estimate, values) {
	data.table(
		  estimate = estimate
		, bootstrap_se = sd(values)
		, p025 = quantile(values, 0.025, names = FALSE)
		, p975 = quantile(values, 0.975, names = FALSE)
	)
}

# Means and standard deviations of gamma and cost under each weighting.
summary_main <- fread(file.path("output", "result", "gamma_c_summary.csv"))
summary_bs <- read_bootstrap("gamma_c_summary")
weightings <- c("equal_auction", "participant", "capacity")
statistics <- data.table(
	  statistic = c("gamma_mean", "gamma_sd", "cost_mean", "cost_sd")
	, source = c("gamma_mean", "gamma_var", "cost_mean", "cost_var")
	, transform = c("identity", "sqrt", "identity", "sqrt")
)

distribution_summary <- rbindlist(lapply(weightings, function(weighting) {
	rbindlist(lapply(seq_len(nrow(statistics)), function(j) {
		source <- statistics$source[j]
		transform <- statistics$transform[j]
		point <- summary_main[variable == source, get(weighting)]
		values <- vapply(summary_bs, function(draw) {
			draw[variable == source, get(weighting)]
		}, numeric(1))
		if (transform == "sqrt") {
			point <- sqrt(max(point, 0))
			values <- sqrt(pmax(values, 0))
		}
		cbind(
			data.table(weighting = weighting, statistic = statistics$statistic[j]),
			percentile_summary(point, values)
		)
	}))
}))
fwrite(distribution_summary, file.path(rdir, "step2_gmm_distribution_summary.csv"))

# Equal-auction projections match the weighting used for the paper's existing
# cost-trend figure. Bands are pointwise bootstrap percentile intervals.
lp_main <- fread(file.path("output", "result", "gamma_c_lp.csv"))
lp_bs <- read_bootstrap("gamma_c_lp")
t_grid <- seq(0, 5, by = 0.01)
year_grid <- 2010 + t_grid

projection_value <- function(dt, prefix, column, derivative = FALSE) {
	terms <- if (derivative) {
		c("l", "t_l", "t2_l", "t3_l")
	} else {
		c("intercept", "t", "t2", "t3")
	}
	beta <- vapply(terms, function(term) {
		dt[get("variable") == paste0(prefix, "_", term), get(column)]
	}, numeric(1))
	beta[1] + beta[2] * t_grid + beta[3] * t_grid^2 + beta[4] * t_grid^3
}

pointwise_projection <- function(prefix, column, derivative = FALSE,
		gamma_scale = FALSE) {
	point <- projection_value(lp_main, prefix, column, derivative)
	draws <- vapply(lp_bs, function(draw) {
		projection_value(draw, prefix, column, derivative)
	}, numeric(length(t_grid)))
	if (gamma_scale) {
		base_point <- projection_value(lp_main, "loggamma", column, FALSE)
		base_draws <- vapply(lp_bs, function(draw) {
			projection_value(draw, "loggamma", column, FALSE)
		}, numeric(length(t_grid)))
		if (derivative) {
			point <- exp(base_point) * point
			draws <- exp(base_draws) * draws
		} else {
			point <- exp(point)
			draws <- exp(draws)
		}
	}
	data.table(
		  year = year_grid
		, estimate = point
		, p025 = apply(draws, 1, quantile, probs = 0.025, names = FALSE)
		, p975 = apply(draws, 1, quantile, probs = 0.975, names = FALSE)
	)
}

gamma_trend <- pointwise_projection(
	"loggamma", "equal_auction", gamma_scale = TRUE
)[, panel := "Auction-date trend (lead time = 0)"]
gamma_lead <- pointwise_projection(
	"loggamma", "equal_auction", derivative = TRUE, gamma_scale = TRUE
)[, panel := "Lead-time effect by auction date"]
gamma_plot_data <- rbind(gamma_trend, gamma_lead)

cost_trend <- pointwise_projection(
	"cost", "equal_auction"
)[, `:=`(panel = "Auction-date trend (lead time = 0)", specification = "Main")]
cost_lead_main <- pointwise_projection(
	"cost", "equal_auction", derivative = TRUE
)[, `:=`(panel = "Lead-time effect by auction date", specification = "Main")]
cost_lead_undiscount <- pointwise_projection(
	"cost", "equal_auction_undiscount", derivative = TRUE
)[, `:=`(panel = "Lead-time effect by auction date", specification = "Undiscounted")]
cost_plot_data <- rbind(cost_trend, cost_lead_main, cost_lead_undiscount)

fwrite(gamma_plot_data, file.path(rdir, "step2_gmm_gamma_trends.csv"))
fwrite(cost_plot_data, file.path(rdir, "step2_gmm_cost_trends.csv"))

theme_trend <- theme_bw() + theme(
	  legend.title = element_blank()
	, legend.position = "top"
	, panel.grid.minor = element_blank()
)

save_single_trend <- function(dt, path, y_label, show_legend = FALSE) {
	fig <- ggplot(dt, aes(x = year, y = estimate)) +
		geom_ribbon(aes(ymin = p025, ymax = p975), fill = "#56B4E9", alpha = 0.25) +
		geom_line(color = "#0072B2", linewidth = 0.8) +
		scale_x_continuous(breaks = 2010:2015) +
		labs(x = NULL, y = y_label) + theme_trend
	if (!show_legend) fig <- fig + theme(legend.position = "none")
	ggsave(path, plot = fig, width = 8, height = 5.5, units = "cm", dpi = 300)
}

save_single_trend(
	gamma_trend, file.path(rdir, "step2_gmm_gamma_t_l0.png"),
	"Risk aversion"
)
save_single_trend(
	gamma_lead, file.path(rdir, "step2_gmm_gamma_l_by_t.png"),
	"Effect per lead-time year"
)
save_single_trend(
	cost_trend, file.path(rdir, "step2_gmm_cost_t_l0.png"),
	"$/MWh"
)

fig <- ggplot(
	cost_plot_data[panel == "Lead-time effect by auction date"],
	aes(x = year, y = estimate, color = specification, fill = specification)
) +
	geom_ribbon(aes(ymin = p025, ymax = p975), alpha = 0.16, color = NA) +
	geom_line(linewidth = 0.8) +
	scale_color_manual(values = c("Main" = "#000000", "Undiscounted" = "#0072B2")) +
	scale_fill_manual(values = c("Main" = "#999999", "Undiscounted" = "#56B4E9")) +
	scale_x_continuous(breaks = 2010:2015) +
	labs(x = NULL, y = "$/MWh per lead-time year") + theme_trend
ggsave(
	file.path(rdir, "step2_gmm_cost_l_by_t.png"), plot = fig,
	width = 8, height = 5.5, units = "cm", dpi = 300
)

# Discount-factor sensitivity in comparable annualized price units.
auction_data <- fread(file.path(
	pathtodata, "processed", "generation_auction_wind20.csv"
))
auction_lead <- unique(auction_data[, .(auction_id, lead_time)])
delta_grid <- c(0.90, 0.95, 0.99)
contract_length <- 20
delta_summary <- rbindlist(lapply(seq_along(delta_grid), function(j) {
	delta <- delta_grid[j]
	dt <- fread(file.path(
		pathtodata, "intermediate", "params2_sa", "gamma_c_data",
		paste0("sa", 15 + j, ".csv")
	))
	dt <- merge(dt, auction_lead, by = "auction_id")
	dt[, delta_tilde := vapply(lead_time, function(lead) {
		mean(delta^seq(lead, lead + contract_length - 1, by = 1))
	}, numeric(1))]
	data.table(
		  delta = delta
		, delta_tilde = weighted.mean(dt$delta_tilde, dt$capacity)
		, risk_premium = weighted.mean(dt$rp, dt$capacity)
		, risk_premium_annualized = weighted.mean(dt$rp/dt$delta_tilde, dt$capacity)
		, cost = weighted.mean(dt$cost, dt$capacity)
		, risk_premium_cost_ratio = weighted.mean(dt$rp, dt$capacity) /
			weighted.mean(dt$cost, dt$capacity)
	)
}))
fwrite(delta_summary, file.path(rdir, "step2_gmm_delta_scaling.csv"))

# The corresponding Step 1 statistic averages over all retained winners.
step1_data <- auction_data[partial_bid_quantity == 0]
step1_data[, K_log := log1p(prior_started_wind_capacity/100)]
step1_delta_summary <- rbindlist(lapply(seq_along(delta_grid), function(j) {
	delta <- delta_grid[j]
	params <- fread(file.path(
		pathtodata, "intermediate", "params1_sa", paste0("sa", 15 + j, ".csv")
	))
	alpha_gamma0 <- params[variable == "alpha_gamma0", hetero]
	alpha_gammaK <- params[variable == "alpha_gammaK", hetero]
	spot_variance <- params[grepl("^var_r_a", variable), .(
		  auction_id = as.integer(sub("var_r_a", "", variable))
		, variance = hetero
	)]
	dt <- merge(step1_data, spot_variance, by = "auction_id")
	dt[, `:=`(
		  gamma = exp(alpha_gamma0 + alpha_gammaK * K_log)
		, delta_tilde = vapply(lead_time, function(lead) {
			mean(delta^seq(lead, lead + contract_length - 1, by = 1))
		}, numeric(1))
	)]
	dt[, risk_premium := gamma * variance/2]
	data.table(
		  delta = delta
		, delta_tilde = mean(dt$delta_tilde)
		, risk_premium = mean(dt$risk_premium)
		, risk_premium_annualized = mean(dt$risk_premium/dt$delta_tilde)
	)
}))
fwrite(step1_delta_summary, file.path(rdir, "step1_gmm_delta_scaling.csv"))
