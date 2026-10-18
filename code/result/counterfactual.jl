#=
Create compact, paper-facing counterfactual outputs from the complete
production results stored in Dropbox.
=#

using CSV
using DataFrames

cd(dirname(dirname(@__DIR__)))
data_directory = strip(read("data.txt", String))
source_directory = joinpath(
	data_directory, "intermediate", "counterfactual_output")
result_directory = joinpath(pwd(), "output", "result")
mkpath(result_directory)

read_output(name) = CSV.read(joinpath(
	source_directory, "counterfactual_$(name).csv"), DataFrame)
select_existing(table, columns) = select(
	table, [column for column in columns if hasproperty(table, column)])

results = read_output("results")
matched_comparison = read_output("matched_comparison")
metadata = read_output("metadata")

primary_rho = 0.13284679368741764
primary_loggamma_sd = 0.0383370135
primary = results[
	isapprox.(results.rho_target, primary_rho; atol=1e-12) .&
	isapprox.(results.loggamma_sd, primary_loggamma_sd; atol=1e-12), :]
primary_metadata = metadata[
	isapprox.(metadata.rho_target, primary_rho; atol=1e-12) .&
	isapprox.(metadata.loggamma_sd, primary_loggamma_sd; atol=1e-12) .&
	metadata.primary_calibration, :]
@assert nrow(primary_metadata) == 1
CSV.write(joinpath(result_directory, "counterfactual_metadata.csv"),
	select_existing(primary_metadata, [
		:rho_target, :loggamma_sd, :kappa, :K_anchor,
		:expected_objective_capacity,
		:expected_participant_aggregate_capacity,
		:registered_entry_rate, :registered_moment_rmse]))
primary_columns = [
	:benchmark, :solve_pass, :cost_shift, :format, :lambda, :gamma_pm,
	:pm_ce, :pm_expected_net_cost, :pm_variance,
	:entry_rate, :n_participants, :entrant_mean_bid, :entrant_mean_markup,
	:potential_mean_gamma, :entrant_mean_gamma, :winner_mean_gamma,
	:entrant_mean_cost, :winner_mean_cost,
	:potential_low_gamma_share, :entrant_low_gamma_share,
	:winner_low_gamma_share, :winner_mean_q,
	:awarded_capacity, :mean_capacity_overshoot,
	:potential_project_mean_q, :auction_supported_capacity,
	:outside_built_capacity, :total_built_capacity,
	:outside_built_capacity_share, :winner_information_rent, :cost_loss,
	:nonentrant_project_share, :auction_loser_project_share,
	:auction_winner_project_share, :nonentrant_outside_build_share,
	:loser_outside_build_share, :outside_builder_nonentrant_share,
	:outside_capacity_nonentrant_share,
	:nonentrant_mean_gamma, :nonentrant_low_gamma_share,
	:nonentrant_mean_cost, :loser_mean_gamma, :loser_low_gamma_share,
	:loser_mean_cost, :outside_builder_mean_gamma,
	:outside_builder_low_gamma_share, :outside_builder_mean_cost,
	:marginal_winning_probability,
	:marginal_conditional_surplus_low, :marginal_conditional_surplus_high,
	:marginal_auction_value_low, :marginal_auction_value_high,
	:marginal_value_gap_low, :marginal_value_gap_high,
	:marginal_value_to_surplus_low, :marginal_value_to_surplus_high,
	:criterion, :joint_converged, :entry_residual,
	:entry_rmse_residual, :entry_max_residual]
CSV.write(joinpath(result_directory, "counterfactual_primary.csv"),
	select_existing(primary, primary_columns))

# Type-specific common-support residuals for the mechanisms reported in the
# appendix diagnostic table. Retain every cost shift and rho value so each
# table cell can report the full range over the production cost path.
diagnostic_cost_shifts = collect(0.0:-1.0:-10.0)
diagnostic_mechanism =
	(results.format .== "share_choice") .|
	(results.format .== "subsidy") .|
	((results.format .== "fixed_share") .&
		(in.(results.lambda, Ref([0.5, 1.0]))))
type_residuals = results[
	(results.benchmark .== "standard") .&
	(results.gamma_pm .== 0.0) .&
	in.(results.cost_shift, Ref(diagnostic_cost_shifts)) .&
	diagnostic_mechanism, :]
type_residual_columns = [
	:loggamma_sd, :cost_shift, :rho_target, :format, :lambda,
	:marginal_value_residual_low, :marginal_value_residual_high,
	:entry_rmse_residual, :entry_max_residual]
CSV.write(joinpath(
	result_directory, "counterfactual_type_residuals.csv"),
	select_existing(type_residuals, type_residual_columns))

matched_columns = [
	:loggamma_sd, :cost_shift, :rho_target, :gamma_pm,
	:share_choice_winner_mean_q, :matched_fixed_lambda,
	:matched_share_gap, :share_choice_pm_ce, :matched_fixed_pm_ce,
	:share_choice_minus_matched_fixed_pm_ce,
	:share_choice_entry_rate, :matched_fixed_entry_rate,
	:share_choice_minus_matched_fixed_entry_rate,
	:share_choice_entrant_mean_gamma, :matched_fixed_entrant_mean_gamma,
	:share_choice_winner_mean_gamma, :matched_fixed_winner_mean_gamma,
	:share_choice_entrant_mean_cost, :matched_fixed_entrant_mean_cost,
	:share_choice_winner_mean_cost, :matched_fixed_winner_mean_cost,
	:share_choice_total_built_capacity,
	:matched_fixed_total_built_capacity,
	:share_choice_minus_matched_fixed_total_built_capacity]
CSV.write(joinpath(
	result_directory, "counterfactual_matched_comparison.csv"),
	select_existing(matched_comparison, matched_columns))

println("Compact counterfactual outputs written to $result_directory")
