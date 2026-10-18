# Brazil new energy auction counterfactual results

# @author: Konan Hara

# @date: September 24, 2026

# @software: R version 4.5.3
# clear environment
rm(list = ls())

# package loading
library(here)
library(data.table)
library(dplyr)
library(kableExtra)
library(ggplot2)

# output paths
rdir <- here("output", "result")
fdir <- here("output", "figure")
tdir <- here("output", "table")
sdir <- here("output", "snt")

# create target output directories if they do not exist
dir.create(fdir, showWarnings = FALSE, recursive = TRUE)
dir.create(tdir, showWarnings = FALSE, recursive = TRUE)
dir.create(sdir, showWarnings = FALSE, recursive = TRUE)

# load compact paper-facing counterfactual outputs
primary <- fread(file.path(rdir, "counterfactual_primary.csv"))
matched_comparison <- fread(file.path(
	rdir, "counterfactual_matched_comparison.csv"
))
type_residuals <- fread(file.path(
	rdir, "counterfactual_type_residuals.csv"
))
metadata <- fread(file.path(rdir, "counterfactual_metadata.csv"))
gamma_cost_distribution <- fread(file.path(
	rdir, "gamma_c_dist_result.csv"
))

stopifnot(
	all(c(
		  "benchmark"
		, "cost_shift"
		, "format"
		, "lambda"
		, "gamma_pm"
		, "entry_rate"
		, "entrant_mean_gamma"
		, "winner_mean_gamma"
		, "winner_low_gamma_share"
		, "winner_mean_q"
		, "entrant_mean_markup"
		, "winner_information_rent"
		, "awarded_capacity"
		, "mean_capacity_overshoot"
		, "outside_built_capacity"
		, "outside_built_capacity_share"
	) %in% names(primary)),
	all(c(
		  "loggamma_sd"
		, "cost_shift"
		, "rho_target"
		, "format"
		, "lambda"
		, "marginal_value_residual_low"
		, "marginal_value_residual_high"
	) %in% names(type_residuals)),
	all(c(
		  "rho_target"
		, "loggamma_sd"
		, "kappa"
		, "expected_objective_capacity"
		, "expected_participant_aggregate_capacity"
	) %in% names(metadata)),
	nrow(metadata) == 1
)

# formatting helpers
write_snt <- function(x, filename, d = 2, s = 1, pct = FALSE) {
	value <- sprintf(paste0("%.", d, "f"), x * s)
	value <- if (pct) paste0(value, "\\%%") else paste0(value, "%")
	cat(value, file = file.path(sdir, filename))
}

fmt_num <- function(x, d = 2) {
	sprintf(paste0("$%.", d, "f$"), x)
}

fmt_pct <- function(x, d = 1) {
	sprintf(paste0("$%.", d, "f$"), 100 * x)
}

# The compact primary file contains the factual-correlation, point-dispersion
# policy grid. Outcomes other than policymaker welfare are duplicated over
# gamma_pm; retain one copy using the risk-neutral policymaker rows.
policy <- copy(primary[benchmark == "standard" & gamma_pm == 0])
stopifnot(policy[, uniqueN(cost_shift)] == 11)

# Reference primitive specification.
rho_reference <- metadata$rho_target[[1]]
loggamma_sd_reference <- metadata$loggamma_sd[[1]]
loggamma_mean_reference <- gamma_cost_distribution[
	variable == "loggamma_mean", estimate
]
loggamma_var_reference <- gamma_cost_distribution[
	variable == "loggamma_var", estimate
]
stopifnot(
	length(loggamma_mean_reference) == 1,
	length(loggamma_var_reference) == 1
)
gamma_mean_reference <- exp(
	loggamma_mean_reference + 0.5 * loggamma_var_reference
)

matched_reference <- copy(matched_comparison[
	abs(rho_target - rho_reference) < 1e-10 &
	abs(loggamma_sd - loggamma_sd_reference) < 1e-10
])

gamma_pm_reference <- primary[
	benchmark == "standard" & gamma_pm > 0, unique(gamma_pm)]
stopifnot(length(gamma_pm_reference) == 1)

# Policymaker CE is computed per MWh of objective generation.  Full-project
# awards create a small capacity overshoot, so objective capacity is awarded
# capacity net of that overshoot.  Multiplying by its 20-year generation gives
# aggregate CE without attributing the overshoot to the policymaker's target.
contract_years <- 20
hours_per_year <- 8760
objective_capacity_mw <- policy[
	, median(awarded_capacity - mean_capacity_overshoot, na.rm = TRUE)
]
aggregate_100m <- function(per_mwh) {
	per_mwh * objective_capacity_mw * hours_per_year * contract_years / 1e8
}

################################################################################
# sentence files
################################################################################

choice <- copy(policy[format == "share_choice"])

get_choice <- function(shift, variable) {
	choice[cost_shift == shift, get(variable)]
}

write_snt(rho_reference, "cf_rho_reference.tex", 3)
write_snt(gamma_mean_reference, "cf_gamma_mean_reference.tex", 3)
write_snt(loggamma_sd_reference, "cf_loggamma_sd_reference.tex", 3)

write_snt_integer <- function(x, filename) {
	value <- format(round(x), big.mark = ",", scientific = FALSE, trim = TRUE)
	cat(paste0(value, "%"), file = file.path(sdir, filename))
}

expected_objective_capacity <- metadata$expected_objective_capacity[[1]]
expected_participant_capacity <-
	metadata$expected_participant_aggregate_capacity[[1]]
outside_capacity_m8 <- get_choice(-8, "outside_built_capacity")
outside_capacity_m10 <- get_choice(-10, "outside_built_capacity")

write_snt_integer(
	expected_objective_capacity, "cf_expected_objective_capacity.tex")
write_snt_integer(
	expected_participant_capacity, "cf_expected_participant_capacity.tex")
write_snt_integer(
	outside_capacity_m8, "cf_choice_outside_capacity_m8.tex")
write_snt_integer(
	outside_capacity_m10, "cf_choice_outside_capacity_m10.tex")
write_snt(
	outside_capacity_m8 / expected_participant_capacity,
	"cf_choice_outside_capacity_participant_pct_m8.tex", 0, 100, TRUE)
write_snt(
	outside_capacity_m10 / expected_participant_capacity,
	"cf_choice_outside_capacity_participant_pct_m10.tex", 0, 100, TRUE)

primary_kappa <- metadata$kappa[[1]]
stopifnot(is.finite(primary_kappa), primary_kappa > 0)
kappa_exponent <- floor(log10(primary_kappa))
kappa_mantissa <- primary_kappa / 10^kappa_exponent
cat(
	sprintf("%.2f\\times10^{%d}", kappa_mantissa, kappa_exponent),
	file = file.path(sdir, "cf_kappa_primary.tex")
)

write_snt(policy[, min(entry_rate)], "cf_entry_rate_min.tex", 1, 100, TRUE)
write_snt(policy[, max(entry_rate)], "cf_entry_rate_max.tex", 1, 100, TRUE)

for (shift in c(0, -5, -10)) {
	suffix <- paste0("m", abs(shift))
	write_snt(
		get_choice(shift, "winner_mean_q"),
		paste0("cf_choice_winner_q_", suffix, ".tex"),
		1, 100, TRUE
	)
}

for (shift in c(0, -4, -10)) {
	suffix <- paste0("m", abs(shift))
	write_snt(
		get_choice(shift, "winner_information_rent"),
		paste0("cf_choice_winner_markup_", suffix, ".tex"),
		2
	)
}

################################################################################
# matched-share policymaker cost and construction table
################################################################################

policy_reductions <- seq(0, 10, 2)
policy_shifts <- -policy_reductions

matched_policy_rn <- matched_reference[
	gamma_pm == 0 & cost_shift %in% policy_shifts,
	.(
		  cost_shift
		, choice_q = share_choice_winner_mean_q
		, choice_ce_rn = share_choice_pm_ce
		, matched_ce_rn = matched_fixed_pm_ce
		, choice_capacity = share_choice_total_built_capacity
		, matched_capacity = matched_fixed_total_built_capacity
	)
]
matched_markup <- primary[
	benchmark == "standard" & gamma_pm == 0 &
		cost_shift %in% policy_shifts &
		format %in% c("share_choice", "matched_fixed_share"),
	.(
		markup = winner_information_rent
	), by = .(cost_shift, format)
]
matched_markup <- dcast(
	matched_markup, cost_shift ~ format, value.var = "markup"
)
setnames(
	matched_markup,
	c("share_choice", "matched_fixed_share"),
	c("choice_markup", "matched_markup")
)
matched_policy_ra <- matched_reference[
	abs(gamma_pm - gamma_pm_reference) < 1e-10 &
		cost_shift %in% policy_shifts,
	.(
		  cost_shift
		, choice_ce_ra = share_choice_pm_ce
		, matched_ce_ra = matched_fixed_pm_ce
	)
]

policy_table <- merge(matched_policy_rn, matched_policy_ra, by = "cost_shift")
policy_table <- merge(policy_table, matched_markup, by = "cost_shift")
setorder(policy_table, -cost_shift)
policy_table[, `:=`(
	choice_ce_rn = aggregate_100m(choice_ce_rn),
	matched_ce_rn = aggregate_100m(matched_ce_rn),
	choice_ce_ra = aggregate_100m(choice_ce_ra),
	matched_ce_ra = aggregate_100m(matched_ce_ra),
	choice_markup = aggregate_100m(choice_markup),
	matched_markup = aggregate_100m(matched_markup)
)]

# The matched-share certainty-equivalent gap peaks at the $6/MWh cost
# reduction for both policymaker risk attitudes. The two gaps differ by less
# than $1 million, so report their rounded common magnitude in the paper.
matched_ce_gap_m6_million <- policy_table[cost_shift == -6, c(
	choice_ce_rn - matched_ce_rn,
	choice_ce_ra - matched_ce_ra
)] * 100
stopifnot(
	length(matched_ce_gap_m6_million) == 2,
	all(matched_ce_gap_m6_million >= 0),
	diff(range(matched_ce_gap_m6_million)) < 1
)
write_snt(
	mean(matched_ce_gap_m6_million),
	"cf_matched_ce_gap_m6_million.tex", 0
)

matched_markup_gap_m6_million <- policy_table[
	cost_shift == -6, choice_markup - matched_markup] * 100
choice_winner_cost_m6_billion <- policy[
	cost_shift == -6 & format == "share_choice", winner_mean_cost
] * objective_capacity_mw * hours_per_year * contract_years / 1e9
stopifnot(
	length(matched_markup_gap_m6_million) == 1,
	matched_markup_gap_m6_million > 0,
	length(choice_winner_cost_m6_billion) == 1,
	choice_winner_cost_m6_billion > 0
)
write_snt(
	matched_markup_gap_m6_million,
	"cf_matched_markup_gap_m6_million.tex", 0
)
write_snt(
	choice_winner_cost_m6_billion,
	"cf_choice_winner_cost_m6_billion.tex", 1
)

policy_value <- function(variable, shift, digits) {
	value <- policy_table[cost_shift == -shift, get(variable)]
	stopifnot(length(value) == 1)
	fmt_num(value, digits)
}

policy_table_output <- data.table(
	Outcome = c(
		  "Winner mean bid share"
		, "Bidder share choice"
		, "Matched fixed share"
		, "Bidder share choice"
		, "Matched fixed share"
		, "Bidder share choice"
		, "Matched fixed share"
		, "Bidder share choice"
		, "Matched fixed share"
	)
)
for (reduction in policy_reductions) {
	policy_table_output[, (as.character(reduction)) := c(
		  policy_value("choice_q", reduction, 2)
		, policy_value("choice_capacity", reduction, 1)
		, policy_value("matched_capacity", reduction, 1)
		, policy_value("choice_markup", reduction, 2)
		, policy_value("matched_markup", reduction, 2)
		, policy_value("choice_ce_rn", reduction, 2)
		, policy_value("matched_ce_rn", reduction, 2)
		, policy_value("choice_ce_ra", reduction, 2)
		, policy_value("matched_ce_ra", reduction, 2)
	)]
}

kbl(
	policy_table_output,
	format = "latex",
	align = paste0("l", strrep("c", length(policy_reductions))),
	escape = FALSE,
	booktabs = TRUE,
	linesep = ""
) %>%
	add_header_above(c(
		  " " = 1
		, "Physical Cost Reduction ($/MWh)" = length(policy_reductions)
	), escape = TRUE) %>%
	pack_rows("Total capacity built (MW)", 2, 3,
		bold = FALSE, italic = TRUE, escape = FALSE) %>%
	pack_rows("Markup (\\$100 million)", 4, 5,
		bold = FALSE, italic = TRUE, escape = FALSE) %>%
	pack_rows("Risk-neutral policymaker CE (\\$100 million)", 6, 7,
		bold = FALSE, italic = TRUE, escape = FALSE) %>%
	pack_rows("Risk-averse policymaker CE (\\$100 million)", 8, 9,
		bold = FALSE, italic = TRUE, escape = FALSE) %>%
	save_kable(file.path(tdir, "counterfactual_policy.tex"))

################################################################################
# type-specific common-support approximation diagnostics
################################################################################

residual_primary <- copy(type_residuals[
	abs(rho_target - rho_reference) < 1e-10 &
		abs(loggamma_sd - loggamma_sd_reference) < 1e-10
])
low_type_value_above_kappa <- exp(
	max(residual_primary$marginal_value_residual_low)
) - 1
high_type_value_below_kappa <- 1 - exp(
	min(residual_primary$marginal_value_residual_high)
)
write_snt(
	low_type_value_above_kappa,
	"cf_common_support_low_value_above_max_pct.tex", 1, 100, TRUE
)
write_snt(
	high_type_value_below_kappa,
	"cf_common_support_high_value_below_max_pct.tex", 1, 100, TRUE
)
################################################################################
# plotting helpers
################################################################################

inside_legend_theme <- function(
	show_legend, position = c(0.02, 0.98), justification = c(0, 1)
) {
	if (!show_legend) {
		return(theme(legend.position = "none"))
	}
	theme(
		  legend.position = "inside"
		, legend.position.inside = position
		, legend.justification = justification
		, legend.background = element_rect(
			fill = scales::alpha("white", 0.72), color = NA
		)
		, legend.title = element_text(size = 8)
		, legend.text = element_text(size = 8)
		, legend.key.height = grid::unit(0.35, "cm")
		, legend.key.width = grid::unit(0.60, "cm")
	)
}

counterfactual_theme <- theme_bw() +
	theme(panel.grid.minor = element_blank())

# Compare aggregate winner markups under bidder share choice with a fixed
# share set equal to the corresponding mean winning share under bidder choice.
matched_markup_paths <- melt(
	policy_table[, .(
		  cost_reduction = -cost_shift
		, `Bidder Share Choice` = choice_markup
		, `Matched Fixed Share` = matched_markup
	)],
	id.vars = "cost_reduction",
	variable.name = "mechanism",
	value.name = "markup_100m"
)
matched_markup_levels <- c("Bidder Share Choice", "Matched Fixed Share")
matched_markup_paths[, mechanism := factor(
	mechanism, levels = matched_markup_levels
)]
matched_markup_color <- setNames(
	c("black", hcl.colors(5, palette = "Blues 3")[2]),
	matched_markup_levels
)
matched_markup_linewidth <- setNames(c(1.05, 0.70), matched_markup_levels)
matched_markup_linetype <- setNames(c("solid", "dashed"), matched_markup_levels)

fig_matched_markup <- ggplot(
	matched_markup_paths,
	aes(
		  x = cost_reduction
		, y = markup_100m
		, color = mechanism
		, linewidth = mechanism
		, linetype = mechanism
	)
) +
	geom_line() +
	geom_point(size = 0.9, show.legend = FALSE) +
	scale_x_continuous(breaks = seq(0, 10, 2)) +
	scale_y_continuous(
		limits = c(0, 1.2), breaks = seq(0, 1.2, 0.2),
		labels = sprintf("%.1f", seq(0, 1.2, 0.2)),
		expand = expansion(mult = c(0.015, 0.025))
	) +
	scale_color_manual(values = matched_markup_color) +
	scale_linewidth_manual(values = matched_markup_linewidth) +
	scale_linetype_manual(values = matched_markup_linetype) +
	xlab("Physical Cost Reduction ($/MWh)") +
	ylab("Winner Markup ($100 million)") +
	labs(color = NULL, linewidth = NULL, linetype = NULL) +
	counterfactual_theme +
	inside_legend_theme(TRUE, c(0.02, 0.02), c(0, 0))

# Main mechanism comparison: bidder share choice is conceptually distinct and
# is drawn first as a thick black reference path.  Selected fixed shares use
# ordered blue dashed lines and remain visible when paths overlap.
mechanism_paths <- copy(policy[
	format == "share_choice" |
	format == "subsidy" |
	(format == "fixed_share" & lambda %in% c(0.5, 1))
])
mechanism_paths[, `:=`(
	  cost_reduction = -cost_shift
	, entrant_high_gamma_share = 1 - entrant_low_gamma_share
	, winner_high_gamma_share = 1 - winner_low_gamma_share
	, outside_builder_high_gamma_share = 1 - outside_builder_low_gamma_share
	, mechanism = fcase(
		format == "share_choice", "Bidder Share Choice",
		format == "subsidy", "λ = 0",
		format == "fixed_share" & abs(lambda - 0.5) < 1e-10, "λ = 0.5",
		format == "fixed_share" & abs(lambda - 1) < 1e-10, "λ = 1"
	)
)]
mechanism_levels <- c(
	"Bidder Share Choice", "λ = 0", "λ = 0.5", "λ = 1"
)
mechanism_paths[, mechanism := factor(mechanism, levels = mechanism_levels)]
mechanism_color <- setNames(
	c("black", hcl.colors(5, palette = "Blues 3")[3:1]),
	mechanism_levels
)
mechanism_linewidth <- setNames(c(1.05, 0.45, 0.60, 0.75), mechanism_levels)
mechanism_linetype <- setNames(c("solid", "dashed", "dashed", "dashed"),
	mechanism_levels)

mechanism_path_panel <- function(
	data, variable, y_label, y_limits = NULL, y_breaks = waiver(),
	legend_position = c(0.98, 0.98), legend_justification = c(1, 1),
	y_labels = waiver()
) {
	data <- data[is.finite(get(variable))]
	ggplot(
		data,
		aes(
			  x = cost_reduction
			, y = .data[[variable]]
			, color = mechanism
			, linewidth = mechanism
			, linetype = mechanism
		)
	) +
		geom_line() +
		geom_point(size = 0.9, show.legend = FALSE) +
		scale_x_continuous(breaks = seq(0, 10, 2)) +
		scale_y_continuous(
			limits = y_limits, breaks = y_breaks, labels = y_labels,
			expand = expansion(mult = c(0.015, 0.025))
		) +
		scale_color_manual(values = mechanism_color) +
		scale_linewidth_manual(values = mechanism_linewidth) +
		scale_linetype_manual(values = mechanism_linetype) +
		xlab("Physical Cost Reduction ($/MWh)") +
		ylab(y_label) +
		labs(color = NULL, linewidth = NULL, linetype = NULL) +
		counterfactual_theme +
		inside_legend_theme(TRUE, legend_position, legend_justification)
}

choice_path_panel <- function(
	variable, y_label, y_limits = NULL, y_breaks = waiver(),
	legend_position = c(0.98, 0.98), legend_justification = c(1, 1)
) {
	mechanism_path_panel(
		mechanism_paths[mechanism == "Bidder Share Choice"],
		variable, y_label, y_limits, y_breaks,
		legend_position, legend_justification
	)
}

main_panels <- list(
	participant_choice_primary = choice_path_panel(
		"entrant_high_gamma_share",
		"Pr(High Risk Aversion)",
		c(0.45, 0.55), c(0.45, 0.50, 0.55), c(0.98, 0.98)
	) + theme(legend.position = "none"),
	markup_mechanisms = mechanism_path_panel(
		mechanism_paths, "winner_information_rent",
		"Winner Markup ($/MWh)", c(0, 1.5), seq(0, 1.5, 0.5),
		c(0.02, 0.02), c(0, 0), sprintf("%.1f", seq(0, 1.5, 0.5))
	),
	choice_q_primary = choice_path_panel(
		"winner_mean_q", "Winner Mean Share",
		c(0, 1), seq(0, 1, 0.2), c(0.98, 0.98)
	),
	matched_markup = fig_matched_markup,
	winner_high_gamma_mechanisms = mechanism_path_panel(
		mechanism_paths, "winner_high_gamma_share",
		"Pr(High Risk Aversion)",
		c(0, 1), seq(0, 1, 0.2), c(0.02, 0.98), c(0, 1)
	),
	outside_capacity_mechanisms = choice_path_panel(
		"outside_built_capacity",
		"Capacity Without Support (MW)", NULL, waiver(),
		c(0.02, 0.98), c(0, 1)
	)
)

for (panel_name in names(main_panels)) {
	ggsave(
		file.path(fdir, paste0("counterfactual_", panel_name, ".png")),
		plot = main_panels[[panel_name]],
		width = 8.8, height = 6.4, units = "cm", dpi = 300
	)
}

################################################################################
# high-risk-aversion type shares along the fixed-share policy grid
################################################################################

fixed_selection <- copy(policy[
	cost_shift %in% c(0, -5, -10) &
	(format == "fixed_share" | format == "subsidy")
])
fixed_selection[, `:=`(
	  fixed_share = fifelse(format == "subsidy", 0, lambda)
	, entrant_high_gamma_share = 1 - entrant_low_gamma_share
	, winner_high_gamma_share = 1 - winner_low_gamma_share
)]
fixed_selection[, cost_reduction := factor(
	-cost_shift,
	levels = c(0, 5, 10),
	labels = c("0", "5", "10")
)]
fixed_color <- setNames(
	hcl.colors(5, palette = "Blues 3")[3:1],
	levels(fixed_selection$cost_reduction)
)
fixed_linewidth <- setNames(
	c(0.40, 0.60, 0.80),
	levels(fixed_selection$cost_reduction)
)

fixed_type_panel <- function(
	data, variable, y_label,
	y_limits = c(0.45, 0.65), y_breaks = seq(0.45, 0.65, 0.05)
) {
	ggplot(
		data,
		aes(
			  x = fixed_share
			, y = .data[[variable]]
			, color = cost_reduction
			, linewidth = cost_reduction
		)
	) +
		geom_line() +
		geom_point(size = 0.9, show.legend = FALSE) +
		scale_x_continuous(breaks = seq(0, 1, 0.25)) +
		scale_y_continuous(
			limits = y_limits, breaks = y_breaks,
			expand = expansion(mult = c(0.015, 0.025))
		) +
		scale_color_manual(values = fixed_color) +
		scale_linewidth_manual(values = fixed_linewidth) +
		xlab(expression(paste("Fixed Share (", lambda, ")"))) +
		ylab(y_label) +
		labs(
			color = "Cost Reduction ($/MWh)",
			linewidth = "Cost Reduction ($/MWh)"
		) +
		counterfactual_theme +
		guides(
			color = guide_legend(
				title.position = "top", title.hjust = 0.5, ncol = 1
			)
		) +
		inside_legend_theme(TRUE, c(0.98, 0.02), c(1, 0))
}

fig_participant_fixed_primary <- fixed_type_panel(
	fixed_selection, "entrant_high_gamma_share",
	"Pr(High Risk Aversion)",
	c(0.45, 0.55), c(0.45, 0.50, 0.55)
)
ggsave(
	file.path(fdir, "counterfactual_participant_fixed_primary.png"),
	plot = fig_participant_fixed_primary,
	width = 8.8, height = 6.4, units = "cm", dpi = 300
)

################################################################################
# policymaker outcomes across share choice and selected fixed shares
################################################################################

fixed_welfare <- copy(primary[
	benchmark == "standard" &
	(format == "share_choice" |
	 format == "subsidy" |
	 (format == "fixed_share" & lambda %in% c(0.5, 1)))
])
fixed_welfare[, `:=`(
	  cost_reduction = -cost_shift
	, pm_ce_100m = aggregate_100m(pm_ce)
	, mechanism = fcase(
		format == "share_choice", "Bidder Share Choice",
		format == "subsidy", "λ = 0",
		format == "fixed_share" & abs(lambda - 0.5) < 1e-10, "λ = 0.5",
		format == "fixed_share" & abs(lambda - 1) < 1e-10, "λ = 1"
	)
)]
fixed_welfare[, mechanism := factor(mechanism, levels = mechanism_levels)]

pm_ce_limits <- c(-8, 8)
pm_ce_breaks <- seq(-8, 8, 2)

fig_pm_ce_rn <- mechanism_path_panel(
	fixed_welfare[gamma_pm == 0], "pm_ce_100m",
	"Policymaker CE ($100 million)",
	pm_ce_limits, pm_ce_breaks, c(0.98, 0.98)
)
fig_pm_ce_ra <- mechanism_path_panel(
	fixed_welfare[abs(gamma_pm - gamma_pm_reference) < 1e-10],
	"pm_ce_100m", "Policymaker CE ($100 million)",
	pm_ce_limits, pm_ce_breaks, c(0.98, 0.98)
)

ggsave(
	file.path(fdir, "counterfactual_pm_ce_rn_mechanisms.png"),
	plot = fig_pm_ce_rn, width = 8.8, height = 6.4, units = "cm", dpi = 300
)
ggsave(
	file.path(fdir, "counterfactual_pm_ce_ra_mechanisms.png"),
	plot = fig_pm_ce_ra, width = 8.8, height = 6.4, units = "cm", dpi = 300
)
