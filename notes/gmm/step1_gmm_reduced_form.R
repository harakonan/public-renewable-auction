# Reduced-form evidence for prior started wind capacity
#
# @author: Konan Hara
# @date: September 13, 2026
# @software: R version 4.5

setwd("~/Workspace/project/renewable-auction/")
rm(list = ls())

library(data.table)
library(fixest)

pathtodata <- readLines("data.txt", warn = FALSE)[1]
rdir <- file.path("notes", "gmm")
auction_data <- fread(file.path(
	pathtodata, "processed", "generation_auction_wind20.csv"
))[partial_bid_quantity == 0]
auction_data[, prior_started_wind_capacity_100 :=
	prior_started_wind_capacity/100]

tidy_model <- function(model, specification) {
	result <- as.data.table(coeftable(model), keep.rownames = "term")
	setnames(
		result
		, c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
		, c("estimate", "std_error", "statistic", "p_value")
	)
	result[, `:=`(
		  specification = specification
		, observations = nobs(model)
		, auctions = uniqueN(auction_data$auction_id)
		, parents = uniqueN(auction_data$parent)
		, r_squared = unname(r2(model, "r2"))
	)]
	setcolorder(result, c(
		  "specification", "term", "estimate", "std_error", "statistic"
		, "p_value", "observations", "auctions", "parents", "r_squared"
	))
	result
}

models <- list(
	  baseline = feols(
		bid_share ~ bid_price + mc_operation
		, data = auction_data, cluster = ~auction_id + parent
	)
	, prior_started_wind_capacity = feols(
		bid_share ~ bid_price + mc_operation + prior_started_wind_capacity_100
		, data = auction_data, cluster = ~auction_id + parent
	)
)

fwrite(
	rbindlist(Map(tidy_model, models, names(models)))
	, file.path(rdir, "step1_gmm_reduced_form.csv")
)
