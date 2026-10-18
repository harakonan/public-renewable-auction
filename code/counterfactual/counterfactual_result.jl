# Aggregate the common-support primitive scenarios and create policy
# comparisons used by code/result/counterfactual.jl.

include("counterfactual_prep.jl")

function prefixed_row(choice, alternative, alternative_name)
	base = (
		loggamma_sd=choice.loggamma_sd,
		cost_shift=choice.cost_shift,
		rho_target=choice.rho_target,
		gamma_pm=choice.gamma_pm,
		share_choice_winner_mean_q=choice.winner_mean_q)
	variables = [
		:pm_ce, :entry_rate, :entrant_mean_gamma, :winner_mean_gamma,
		:entrant_mean_cost, :winner_mean_cost,
		:entrant_low_gamma_share, :winner_low_gamma_share,
		:total_built_capacity]
	result = Dict{Symbol, Any}(pairs(base))
	if alternative_name == "best_fixed"
		result[:best_fixed_lambda] = alternative.lambda
	else
		result[:matched_fixed_lambda] = alternative.lambda
		result[:matched_share_gap] = alternative.lambda - choice.winner_mean_q
	end
	for variable in variables
		choice_name = Symbol("share_choice_", variable)
		alternative_variable = Symbol(alternative_name, "_", variable)
		difference_name = Symbol(
			"share_choice_minus_", alternative_name, "_", variable)
		result[choice_name] = getproperty(choice, variable)
		result[alternative_variable] = getproperty(alternative, variable)
		result[difference_name] = getproperty(choice, variable) -
			getproperty(alternative, variable)
	end
	return (; result...)
end

function main()
	task_root = get(ENV, "CF_TASK_OUTPUT_ROOT", intermediate_directory)
	grid_count = length(parse_grid("CF_RHO_GRID", join([
		rho_sensitivity_lower, 0.0, representative_primitive.factual_rho,
		rho_sensitivity_upper], ","))) * length(parse_grid(
		"CF_LOGGAMMA_SD_GRID", "0.0383370135,0.0884369676"))
	frames = DataFrame[]
	metadata_frames = DataFrame[]
	for task in 1:grid_count
		directory = joinpath(task_root,
			@sprintf("counterfactual_task_%03d", task))
		result_path = joinpath(directory, "results.csv")
		metadata_path = joinpath(directory, "metadata.csv")
		isfile(result_path) || error("Missing task result: $result_path")
		isfile(metadata_path) || error("Missing task metadata: $metadata_path")
		push!(frames, CSV.read(result_path, DataFrame))
		push!(metadata_frames, CSV.read(metadata_path, DataFrame))
	end
	results = vcat(frames...; cols=:union)
	metadata = vcat(metadata_frames...; cols=:union)
	sort!(results, [:loggamma_sd, :rho_target, :cost_shift,
		:format, :lambda, :gamma_pm])

	comparison = NamedTuple[]
	matched_comparison = NamedTuple[]
	keys = unique(select(results,
		[:loggamma_sd, :rho_target, :cost_shift, :gamma_pm]))
	for key in eachrow(keys)
		cell = results[
			isapprox.(results.loggamma_sd, key.loggamma_sd; atol=1e-12) .&
			isapprox.(results.rho_target, key.rho_target; atol=1e-12) .&
			isapprox.(results.cost_shift, key.cost_shift; atol=1e-12) .&
			isapprox.(results.gamma_pm, key.gamma_pm; atol=1e-12), :]
		choice = only(eachrow(cell[cell.format .== "share_choice", :]))
		fixed = cell[in.(cell.format, Ref(["fixed_share", "subsidy"])), :]
		best = fixed[argmin(fixed.pm_ce), :]
		matched = only(eachrow(cell[cell.format .== "matched_fixed_share", :]))
		push!(comparison, prefixed_row(choice, best, "best_fixed"))
		push!(matched_comparison,
			prefixed_row(choice, matched, "matched_fixed"))
	end
	comparison_table = DataFrame(comparison)
	matched_table = DataFrame(matched_comparison)

	equilibria = results[results.gamma_pm .== 0.0, :]
	diagnostics = combine(groupby(equilibria,
		[:loggamma_sd, :format, :lambda]),
		nrow => :n_scenarios,
		:criterion => median => :median_criterion,
		:criterion => maximum => :max_criterion,
		:entry_rmse_residual => maximum => :max_entry_residual,
		:entry_max_residual => maximum => :max_entry_distribution_residual,
		:ode_probability_residual => maximum => :max_ode_probability_residual,
		:joint_converged => sum => :n_joint_converged,
		:ode_probability_monotone => all => :all_ode_probability_monotone,
		:ode_bid_monotone => all => :all_ode_bid_monotone,
		:ode_K_positive => all => :all_ode_K_positive)
	insertcols!(diagnostics, 2, :benchmark => fill("standard", nrow(diagnostics)))

	fixed = equilibria[in.(equilibria.format,
		Ref(["fixed_share", "subsidy"])), :]
	lambda_diagnostic_rows = NamedTuple[]
	for cell in groupby(fixed, [:loggamma_sd, :rho_target, :cost_shift])
		ordered = sort(cell, :lambda)
		winner_change = diff(ordered.winner_low_gamma_share)
		entry_change = diff(ordered.entry_rate)
		log_K_change = diff(log.(ordered.upper_K))
		push!(lambda_diagnostic_rows, (
			loggamma_sd=ordered.loggamma_sd[1],
			rho_target=ordered.rho_target[1],
			cost_shift=ordered.cost_shift[1],
			n_lambda=nrow(ordered),
			max_adjacent_entry_change=maximum(abs.(entry_change)),
			max_adjacent_winner_low_gamma_change=
				maximum(abs.(winner_change)),
			max_winner_low_gamma_reversal=
				max(0.0, maximum(winner_change)),
			max_adjacent_log_K_change=maximum(abs.(log_K_change)),
			max_entry_distribution_residual=
				maximum(ordered.entry_max_residual)))
	end
	lambda_diagnostics = DataFrame(lambda_diagnostic_rows)

	output_directory = get(ENV, "CF_AGGREGATE_OUTPUT_ROOT", joinpath(
		dirname(intermediate_directory), "counterfactual_output"))
	mkpath(output_directory)
	CSV.write(joinpath(output_directory, "counterfactual_results.csv"), results)
	CSV.write(joinpath(output_directory, "counterfactual_entry.csv"), select(
		results, [:loggamma_sd, :rho_target, :cost_shift, :format, :lambda,
			:gamma_pm, :entry_rate, :entrant_mean_cost,
			:entrant_mean_gamma, :entrant_low_gamma_share]))
	CSV.write(joinpath(output_directory, "counterfactual_comparison.csv"),
		comparison_table)
	CSV.write(joinpath(output_directory,
		"counterfactual_matched_comparison.csv"), matched_table)
	CSV.write(joinpath(output_directory, "counterfactual_matched_fixed.csv"),
		results[results.format .== "matched_fixed_share", :])
	CSV.write(joinpath(output_directory, "counterfactual_diagnostics.csv"),
		diagnostics)
	CSV.write(joinpath(output_directory,
		"counterfactual_lambda_diagnostics.csv"), lambda_diagnostics)
	CSV.write(joinpath(output_directory, "counterfactual_metadata.csv"), metadata)
	println("Aggregated $grid_count primitive scenarios into $output_directory")
	println("Common-support pooled marginal-value solves converged: ",
		sum(equilibria.joint_converged), " / ", nrow(equilibria))
end

main()
