# Solve one primitive scenario in the two-risk-type common-support model.

include("counterfactual_model.jl")

using Serialization

function scenario_grid()
	rhos = parse_grid("CF_RHO_GRID", join([
		rho_sensitivity_lower, 0.0, representative_primitive.factual_rho,
		rho_sensitivity_upper], ","))
	sds = parse_grid("CF_LOGGAMMA_SD_GRID",
		"0.0383370135,0.0884369676")
	return [(rho=rho, loggamma_sd=sd) for sd in sds for rho in rhos]
end

function solve_path(format, lambda_path, shifts, weights, cache, map;
		initial_prediction=nothing)
	rows = NamedTuple[]
	boundaries = Any[]
	prediction = initial_prediction
	for (index, shift) in enumerate(shifts)
		draws, achieved_rho, outside = two_type_draws(
			weights.loggamma_sd, weights.rho_target, shift)
		primitives = primitive_distribution(draws, weights.potential_weight)
		lambda = lambda_path[index]
		if isnothing(prediction)
			prediction = [cache.upper_bid, log(cache.K_anchor)]
		end
		boundary = regularized_boundary(primitives, format, lambda, map,
			cache.kappa, prediction)
		path = integrate_strategy(boundary, primitives, format, lambda, map;
			maximum_bid_step=parse(Float64,
				get(ENV, "CF_ODE_MAXIMUM_BID_STEP", "0.10")),
			maximum_probability_step=parse(Float64,
				get(ENV, "CF_ODE_MAXIMUM_PROBABILITY_STEP", "0.002")))
		selection = selection_probability(draws, primitives, boundary)
		entry_weight = weights.potential_weight .* selection
		inversion = strategy_inversion(
			draws, primitives, boundary, path, format, lambda)
		entry = participation_statistics(inversion, draws, environment;
			entry_weight=entry_weight,
			potential_weight=weights.potential_weight,
			low_gamma_cutoff=mean(draws.gamma))
		allocation = allocation_statistics(
			(details=(implied=inversion,),), format, lambda, draws,
			environment, controls;
			entry_weight=entry_weight,
			potential_weight=weights.potential_weight,
			outside_active=outside.active,
			inversion=inversion,
			n_participants=environment.n_registered * sum(entry_weight),
			low_gamma_cutoff=mean(draws.gamma))
		entry_rate = sum(entry_weight)
		entrant_mean_cost = sum(entry_weight .* draws.cost) / entry_rate
		potential_mean_cost = sum(
			weights.potential_weight .* draws.cost) / sum(weights.potential_weight)
		type_residual = log.(boundary.value ./ cache.kappa)
		value_to_surplus = boundary.value ./ boundary.surplus
		for gamma_pm in (0.0, mean(draws.gamma))
			push!(rows, merge((
				benchmark="standard", solve_pass="production",
				loggamma_sd=weights.loggamma_sd,
				cost_shift=shift, rho_target=weights.rho_target,
				rho_achieved=achieved_rho, format=String(format),
				lambda=format == :share_choice ? NaN : lambda,
				kappa=cache.kappa,
				gamma_pm=gamma_pm,
				pm_ce=allocation.pm_expected_net_cost +
					0.5 * gamma_pm * allocation.pm_variance,
				entry_rate=entry_rate,
				n_participants=environment.n_registered * entry_rate,
				potential_mean_cost=potential_mean_cost,
				entrant_mean_cost=entrant_mean_cost,
				upper_bid=boundary.bid, upper_K=boundary.K,
				cutoff_low=boundary.x[1], cutoff_high=boundary.x[2],
				marginal_winning_probability=boundary.winning_probability,
				marginal_conditional_surplus_low=boundary.surplus[1],
				marginal_conditional_surplus_high=boundary.surplus[2],
				marginal_auction_value_low=boundary.value[1],
				marginal_auction_value_high=boundary.value[2],
				marginal_value_gap_low=boundary.value[1] - cache.kappa,
				marginal_value_gap_high=boundary.value[2] - cache.kappa,
				marginal_value_to_surplus_low=value_to_surplus[1],
				marginal_value_to_surplus_high=value_to_surplus[2],
				marginal_value_residual_low=type_residual[1],
				marginal_value_residual_high=type_residual[2],
				pooled_marginal_value_residual=boundary.value_residual,
				continuation_distance=boundary.continuation_distance,
				criterion=boundary.objective,
				joint_converged=abs(boundary.value_residual) <= 5e-4,
				entry_residual=boundary.value_residual,
				entry_rmse_residual=sqrt(mean(type_residual.^2)),
				entry_max_residual=maximum(abs.(type_residual)),
				ode_probability_monotone=all(diff([
					state.probability for state in path.state]) .>= -1e-8),
				ode_probability_residual=maximum(abs.([
					state.probability for state in path.state] .-
					path.probability)),
				ode_bid_monotone=all(diff(path.bid) .>= -1e-10),
				ode_K_positive=minimum(path.K) > 0), entry, allocation))
		end
		push!(boundaries, boundary)
		prediction = boundary.parameter
	end
	return rows, boundaries
end

function main()
	cache_path = get(ENV, "CF_BASELINE_FILE", joinpath(
		intermediate_directory, "counterfactual_baseline.jls"))
	isfile(cache_path) || error("Baseline cache not found: $cache_path")
	cache = deserialize(cache_path)
	cache.cache_kind ==
		"two_type_common_support_marginal_participant_bid_anchor" || error(
		"The baseline cache does not use the marginal-participant anchor. " *
		"Run the production pipeline with CF_FRESH=1.")
	grid = scenario_grid()
	task = parse(Int, get(ENV, "CF_TASK_ID", "1"))
	1 <= task <= length(grid) || error("CF_TASK_ID is outside the scenario grid.")
	scenario = grid[task]
	map = allocation_map()
	draws0, _, _ = two_type_draws(
		scenario.loggamma_sd, scenario.rho, 0.0)
	primary = isapprox(scenario.rho, cache.factual_rho; atol=1e-12) &&
		isapprox(scenario.loggamma_sd, cache.primary_loggamma_sd; atol=1e-12)
	calibration = calibrate_registered_pool(
		draws0, map, cache.upper_bid_anchor;
		alpha_start=cache.registered_tilt_alpha)
	# Treat each sensitivity specification as the factual primitive model when
	# calibrating entry.  Its implied kappa is then held fixed over every policy
	# format and cost shift within that specification.
	scenario_cache = merge(cache, (
		kappa=calibration.kappa,
		upper_bid=calibration.boundary.bid,
		K_anchor=calibration.boundary.K))
	potential_weight = calibration.potential_weight
	weights = (rho_target=scenario.rho,
		loggamma_sd=scenario.loggamma_sd,
		potential_weight=potential_weight)
	shifts = sort(parse_grid("CF_COST_SHIFT_GRID",
		"0,-1,-2,-3,-4,-5,-6,-7,-8,-9,-10");
		rev=true)
	rows = NamedTuple[]
	share_rows, share_boundaries = solve_path(
		:share_choice, fill(NaN, length(shifts)), shifts, weights,
		scenario_cache, map;
		initial_prediction=calibration.boundary.parameter)
	append!(rows, share_rows)
	share_q = [only([row.winner_mean_q for row in share_rows
		if row.cost_shift == shift && row.gamma_pm == 0.0]) for shift in shifts]

	lambdas = sort(parse_grid("CF_LAMBDA_GRID",
		join(collect(0.0:0.05:1.0), ",")); rev=true)
	for lambda in lambdas
		format = lambda == 0.0 ? :subsidy : :fixed_share
		fixed_rows, _ = solve_path(format, fill(lambda, length(shifts)),
			shifts, weights, scenario_cache, map;
			initial_prediction=share_boundaries[1].parameter)
		append!(rows, fixed_rows)
	end
	matched_rows, _ = solve_path(:fixed_share, share_q, shifts,
		weights, scenario_cache, map;
		initial_prediction=share_boundaries[1].parameter)
	matched_rows = [merge(row, (format="matched_fixed_share",))
		for row in matched_rows]
	append!(rows, matched_rows)

	output_root = get(ENV, "CF_TASK_OUTPUT_ROOT", intermediate_directory)
	output_directory = joinpath(output_root,
		@sprintf("counterfactual_task_%03d", task))
	mkpath(output_directory)
	CSV.write(joinpath(output_directory, "results.csv"), DataFrame(rows))
	metadata = DataFrame(
		task=[task], rho_target=[scenario.rho],
		loggamma_sd=[scenario.loggamma_sd],
		expected_objective_capacity=[environment.objective_capacity],
		expected_participant_aggregate_capacity=[
			environment.n_participants * mean_capacity],
		registered_entry_rate=[calibration.entry_rate],
		registered_entry_rate_residual=[calibration.entry_residual],
		registered_moment_rmse=[calibration.moment_rmse],
		registered_moment_maximum=[calibration.moment_maximum],
		kappa=[calibration.kappa], K_anchor=[calibration.K_anchor],
		primary_calibration=[primary])
	CSV.write(joinpath(output_directory, "metadata.csv"), metadata)
	println("Completed scenario $task / $(length(grid)): ",
		"rho=$(scenario.rho), loggamma_sd=$(scenario.loggamma_sd), ",
		"rows=$(length(rows))")
end

main()
