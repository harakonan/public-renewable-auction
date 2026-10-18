# Unweighted first-stage GMM goodness-of-fit figures
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
fdir <- rdir
dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
figure_path <- function(stem) {
	file.path(fdir, paste0(stem, ".png"))
}

diagnostics <- fread(file.path(
	pathtodata, "intermediate", "bid_share_gmm_diagnostics.csv"
))
predictive_cdf <- fread(file.path(
	rdir, "step1_gmm_predictive_cdf.csv"
))
predictive_cdf[, `:=`(
	  observed_cdf = pmin(pmax(observed_cdf, 0), 1)
	, predicted_cdf = pmin(pmax(predicted_cdf, 0), 1)
)]

model_levels <- "hetero"
model_labels <- "Heterogeneous"
format_levels <- c("all", "sealed_bid", "descending_clock")
format_labels <- c("Combined", "Sealed bid", "Descending clock")
series_colors <- c("Observed" = "#000000", "Predicted" = "#0072B2")

theme_gof <- theme_bw() +
	theme(
	  legend.title = element_blank()
	, legend.position = "top"
	, panel.grid.minor = element_blank()
	, strip.text.y = element_text(size = 8)
	)

fit_summary <- function(dt, group_variables) {
	dt[, .(
		  observations = .N
		, observed_mean = mean(bid_share)
		, predicted_mean = mean(predicted_mean)
		, mean_error = mean(bid_share - predicted_mean)
		, rmse = sqrt(mean((bid_share - predicted_mean)^2))
		, observed_lower = mean(observed_lower)
		, predicted_lower = mean(predicted_lower)
		, observed_interior = mean(observed_interior)
		, predicted_interior = mean(predicted_interior)
		, observed_upper = mean(observed_upper)
		, predicted_upper = mean(predicted_upper)
		, category_l1_error =
			abs(mean(observed_lower) - mean(predicted_lower)) +
			abs(mean(observed_interior) - mean(predicted_interior)) +
			abs(mean(observed_upper) - mean(predicted_upper))
	), by = group_variables]
}

fit_by_format <- rbind(
	  fit_summary(diagnostics, "model")[, auction_format := "all"]
	, fit_summary(diagnostics, c("model", "auction_format"))
	, fill = TRUE
)
setcolorder(fit_by_format, c("model", "auction_format"))
fwrite(fit_by_format, file.path(
	rdir, "step1_gmm_fit.csv"
))

fit_by_auction <- fit_summary(
	  diagnostics
	, c("model", "auction_id", "auction_type", "auction_format", "contract_type")
)
fwrite(fit_by_auction, file.path(
	rdir, "step1_gmm_fit_by_auction.csv"
))

# Exact predictive histograms, including censoring masses and eta dispersion.
bin_upper <- c(0.25, seq(0.35, 0.95, by = 0.10), 1.00)
bin_labels <- c(
	  "[0, .25]"
	, paste0(
		"(", sprintf("%.2f", bin_upper[-length(bin_upper)]), ", ",
		sprintf("%.2f", bin_upper[-1]), "]"
	  )
)
cdf_format <- predictive_cdf[
	group_type == "all" | group_type == "auction_format"
]
cdf_format[, auction_format := fifelse(
	group_type == "all", "all", group_value
)]
hist_dt <- cdf_format[round(bid_share, 2) %in% round(bin_upper, 2)]
hist_dt[, bid_share := round(bid_share, 2)]
setorder(hist_dt, model, auction_format, bid_share)
hist_dt[, `:=`(
	  observed_probability = observed_cdf - shift(observed_cdf, fill = 0)
	, predicted_probability = predicted_cdf - shift(predicted_cdf, fill = 0)
), by = c("model", "auction_format")]
hist_dt[, bin := factor(
	bin_labels[match(round(bid_share, 2), round(bin_upper, 2))]
	, levels = bin_labels
)]
hist_plot <- melt(
	  hist_dt
	, id.vars = c("model", "auction_format", "bin")
	, measure.vars = c("observed_probability", "predicted_probability")
	, variable.name = "series"
	, value.name = "probability"
)
hist_plot[, `:=`(
	  model = factor(model, levels = model_levels, labels = model_labels)
	, auction_format = factor(
		auction_format, levels = format_levels, labels = format_labels
	  )
	, series = factor(
		series
		, levels = c("observed_probability", "predicted_probability")
		, labels = c("Observed", "Predicted")
	  )
)]

fig <- ggplot(hist_plot, aes(x = bin, y = probability, fill = series)) +
	geom_col(position = position_dodge(width = 0.8), width = 0.75) +
	facet_grid(model ~ auction_format) +
	scale_fill_manual(values = series_colors) +
	scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
	labs(x = "Bid share bin", y = "Winner share") +
	theme_gof +
	theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(
	  figure_path("bidshare_gmm_hist_gof")
	, plot = fig, width = 18, height = 14, units = "cm", dpi = 300
)

# Auction-level mean-share fit.
mean_auction <- melt(
	  fit_by_auction
	, id.vars = c("model", "auction_id", "auction_type", "auction_format")
	, measure.vars = c("observed_mean", "predicted_mean")
	, variable.name = "series"
	, value.name = "bid_share"
)
mean_auction[, `:=`(
	  model = factor(model, levels = model_levels, labels = model_labels)
	, series = factor(
		series
		, levels = c("observed_mean", "predicted_mean")
		, labels = c("Observed", "Predicted")
	  )
)]

fig <- ggplot(mean_auction, aes(
	x = auction_id, y = bid_share, color = series, shape = series
)) +
	geom_line(aes(group = series), linewidth = 0.6) +
	geom_point(size = 2) +
	geom_vline(xintercept = 16.5, color = "grey60", linetype = "dashed") +
	facet_wrap(~model, ncol = 1) +
	scale_color_manual(values = series_colors) +
	scale_x_continuous(breaks = sort(unique(mean_auction$auction_id))) +
	scale_y_continuous(limits = c(0, 1)) +
	labs(x = "Auction ID", y = "Mean bid share") +
	theme_gof
ggsave(
	  figure_path("bidshare_gmm_auction_mean_gof")
	, plot = fig, width = 17, height = 13, units = "cm", dpi = 300
)

# Descending-clock auction-level censoring fit.
dc_auction <- fit_by_auction[auction_format == "descending_clock"]
observed_dc <- melt(
	  dc_auction
	, id.vars = c("model", "auction_id")
	, measure.vars = c("observed_lower", "observed_upper")
	, variable.name = "category"
	, value.name = "probability"
)
observed_dc[, series := "Observed"]
predicted_dc <- melt(
	  dc_auction
	, id.vars = c("model", "auction_id")
	, measure.vars = c("predicted_lower", "predicted_upper")
	, variable.name = "category"
	, value.name = "probability"
)
predicted_dc[, series := "Predicted"]
dc_plot <- rbind(observed_dc, predicted_dc)
dc_plot[, `:=`(
	  model = factor(model, levels = model_levels, labels = model_labels)
	, category = factor(
		sub("^(observed|predicted)_", "", category)
		, levels = c("lower", "upper")
		, labels = c("Lower bound", "Upper bound")
	  )
	, series = factor(series, levels = c("Observed", "Predicted"))
)]

fig <- ggplot(dc_plot, aes(
	x = auction_id, y = probability, color = series, shape = series
)) +
	geom_line(aes(group = series), linewidth = 0.6) +
	geom_point(size = 2) +
	facet_grid(model ~ category) +
	scale_color_manual(values = series_colors) +
	scale_x_continuous(breaks = sort(unique(dc_plot$auction_id))) +
	scale_y_continuous(
		limits = c(0, 1), labels = scales::percent_format(accuracy = 1)
	) +
	labs(x = "Descending-clock auction ID", y = "Winner share") +
	theme_gof
ggsave(
	  figure_path("bidshare_gmm_dc_auction_bounds_gof")
	, plot = fig, width = 16, height = 13, units = "cm", dpi = 300
)
