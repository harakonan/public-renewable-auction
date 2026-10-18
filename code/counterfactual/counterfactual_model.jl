# Two-risk-type common-upper-support production model.

include("counterfactual_prep.jl")

using Serialization

struct CFKernelDistribution
	x::Vector{Float64}
	weight::Vector{Float64}
	bandwidth::Float64
end

function two_type_draws(loggamma_sd, rho, cost_shift)
	n = controls.n_aux
	iseven(n) || error("CF_N_AUX must be even for the two-type model.")
	class_score = repeat([-1.0, 1.0], n ÷ 2)
	shuffle!(MersenneTwister(controls.seed), class_score)
	distance = loggamma_sd * sqrt((n - 1) / n)
	log_center = log(counterfactual_gamma_mean) - log(cosh(distance))
	loggamma = log_center .+ distance .* class_score
	gamma = exp.(loggamma)
	z_loggamma = (loggamma .- mean(loggamma)) ./ std(loggamma)
	z_gamma = (gamma .- mean(gamma)) ./ std(gamma)
	noise = standardized_draws(MersenneTwister(controls.seed + 1), n)
	noise .-= dot(noise, z_gamma) / dot(z_gamma, z_gamma) .* z_gamma
	noise = (noise .- mean(noise)) ./ std(noise)
	cost_score = rho .* z_gamma .+ sqrt(max(1 - rho^2, 0.0)) .* noise
	cost_score = (cost_score .- mean(cost_score)) ./ std(cost_score)
	cost = representative_primitive.cost_mean .+
		sqrt(representative_primitive.cost_var) .* cost_score
	capacity_index = systematic_weighted_indices(
		MersenneTwister(controls.seed + 3), primitive_weight, n)
	capacity = Float64.(primitive_table.capacity[capacity_index])
	draws = CFDraws(gamma, cost, capacity, z_loggamma)
	outside = apply_wholesale_outside_option(
		draws.cost, draws.gamma, cost_shift, environment)
	updated = CFDraws(
		draws.gamma, outside.cost, draws.capacity, draws.z_gamma)
	return updated, cor(updated.gamma, updated.cost), outside
end

two_type_class(draws::CFDraws) = ifelse.(draws.z_gamma .< 0, 1, 2)

function interpolate_flat(x, y, value)
	length(x) == length(y) || error("Interpolation inputs differ in length.")
	order = sortperm(x)
	x_sorted = Float64.(x[order])
	y_sorted = Float64.(y[order])
	unique_x = Float64[]
	unique_y = Float64[]
	for j in eachindex(x_sorted)
		if isempty(unique_x) || x_sorted[j] > unique_x[end] + 1e-10
			push!(unique_x, x_sorted[j])
			push!(unique_y, y_sorted[j])
		else
			unique_y[end] = y_sorted[j]
		end
	end
	length(unique_x) == 1 && return unique_y[1]
	value <= unique_x[1] && return unique_y[1]
	value >= unique_x[end] && return unique_y[end]
	j = searchsortedlast(unique_x, value)
	weight = (value - unique_x[j]) / (unique_x[j + 1] - unique_x[j])
	return (1 - weight) * unique_y[j] + weight * unique_y[j + 1]
end

function kernel_distribution(x, weight)
	weight = Float64.(weight)
	weight ./= sum(weight)
	mean_x = sum(weight .* x)
	sd_x = sqrt(sum(weight .* (x .- mean_x).^2))
	effective_n = inv(sum(weight.^2))
	bandwidth = max(1.06 * sd_x * effective_n^(-1 / 5), 0.03)
	return CFKernelDistribution(Float64.(x), weight, bandwidth)
end

function kernel_cdf(distribution::CFKernelDistribution, value)
	total = 0.0
	for j in eachindex(distribution.x)
		total += distribution.weight[j] * cdf(Normal(),
			(value - distribution.x[j]) / distribution.bandwidth)
	end
	return clamp(total, 1e-10, 1 - 1e-10)
end

function kernel_pdf(distribution::CFKernelDistribution, value)
	total = 0.0
	for j in eachindex(distribution.x)
		z = (value - distribution.x[j]) / distribution.bandwidth
		total += distribution.weight[j] * pdf(Normal(), z) /
			distribution.bandwidth
	end
	return max(total, 1e-12)
end

function primitive_distribution(draws::CFDraws, potential_weight)
	weight = potential_weight ./ sum(potential_weight)
	class = two_type_class(draws)
	floor = wholesale_outside_floor(
		draws.gamma, environment.mu_r, environment.sigma2_r)
	x = max.(draws.cost .- floor, 0.0)
	distribution = [kernel_distribution(
		x[class .== g], weight[class .== g]) for g in 1:2]
	class_probability = [sum(weight[class .== g]) for g in 1:2]
	gamma = [mean(draws.gamma[class .== g]) for g in 1:2]
	floor_by_class = wholesale_outside_floor(
		gamma, environment.mu_r, environment.sigma2_r)
	return (distribution=distribution,
		class_probability=class_probability, gamma=gamma,
		floor=floor_by_class, x=x, class=class)
end

function allocation_map()
	n_simulation = controls.n_winprob
	n_competitor = max(environment.n_registered - 1, 1)
	index = systematic_weighted_indices(
		MersenneTwister(controls.seed + 7101), primitive_weight,
		n_simulation * n_competitor)
	capacity = reshape(Float64.(primitive_table.capacity[index]),
		n_simulation, n_competitor)
	objective = objective_capacity_draws(
		environment, n_simulation, MersenneTwister(controls.seed + 7102))
	conditional = zeros(n_competitor + 1)
	for s in 1:n_simulation
		cumulative = 0.0
		conditional[1] += cdf(Normal(), objective[s] / controls.h_capacity)
		for k in 1:n_competitor
			cumulative += capacity[s, k]
			conditional[k + 1] += cdf(Normal(),
				(objective[s] - cumulative) / controls.h_capacity)
		end
	end
	conditional ./= n_simulation
	conditional = decreasing_pava(conditional, ones(length(conditional)))
	function value_derivative(probability)
		p = clamp(probability, 1e-10, 1 - 1e-10)
		binomial = Binomial(n_competitor, p)
		value = 0.0
		derivative = 0.0
		for k in 0:n_competitor
			mass = pdf(binomial, k)
			score = k / p - (n_competitor - k) / (1 - p)
			value += mass * conditional[k + 1]
			derivative += mass * score * conditional[k + 1]
		end
		return clamp(value, 1e-300, 1 - 1e-12),
			min(derivative, -1e-300)
	end
	return (conditional_probability=conditional,
		value_derivative=value_derivative)
end

function format_share_derivative(bid, gamma, format, lambda)
	if format == :share_choice
		raw = 1 - (environment.mu_r_share - environment.delta_tilde * bid) /
			(gamma * environment.sigma2_r)
		if raw <= environment.q_lower
			return environment.q_lower, 0.0
		elseif raw >= 1
			return 1.0, 0.0
		end
		return raw, environment.delta_tilde /
			(gamma * environment.sigma2_r)
	elseif format == :subsidy
		return 0.0, 0.0
	elseif format == :fixed_share
		return lambda, 0.0
	end
	error("Unknown auction format: $format")
end

function inverse_state(bid, K, primitives, format, lambda)
	K > 0 || error("K must be positive.")
	x = zeros(2)
	surplus = zeros(2)
	q = zeros(2)
	x_bid = zeros(2)
	x_K = zeros(2)
	for g in 1:2
		gamma = primitives.gamma[g]
		q[g], q_derivative = format_share_derivative(
			bid, gamma, format, lambda)
		marginal_payment = format == :share_choice ?
			q[g] * environment.delta_tilde : environment.delta_tilde
		gross = if format == :share_choice
			q[g] * environment.delta_tilde * bid +
				(1 - q[g]) * environment.mu_r -
				0.5 * gamma * (1 - q[g])^2 * environment.sigma2_r
		else
			environment.delta_tilde * bid +
				(1 - q[g]) * environment.mu_r -
				0.5 * gamma * (1 - q[g])^2 * environment.sigma2_r
		end
		increment = gamma * K * marginal_payment
		surplus[g] = log1p(increment) / gamma
		x[g] = gross - surplus[g] - primitives.floor[g]
		gross_bid = marginal_payment
		surplus_bid = format == :share_choice ?
			K * environment.delta_tilde * q_derivative / (1 + increment) : 0.0
		x_bid[g] = gross_bid - surplus_bid
		x_K[g] = -marginal_payment / (1 + increment)
	end
	probability = 0.0
	probability_bid = 0.0
	probability_K = 0.0
	for g in 1:2
		F = kernel_cdf(primitives.distribution[g], x[g])
		f = kernel_pdf(primitives.distribution[g], x[g])
		mass = primitives.class_probability[g]
		probability += mass * F
		probability_bid += mass * f * x_bid[g]
		probability_K += mass * f * x_K[g]
	end
	return (x=x, surplus=surplus, q=q,
		probability=clamp(probability, 1e-8, 1 - 1e-8),
		probability_bid=probability_bid,
		probability_K=probability_K)
end

function boundary_values(bid, K, primitives, format, lambda, map)
	state = inverse_state(bid, K, primitives, format, lambda)
	winning_probability, _ = map.value_derivative(state.probability)
	value = [lottery_ce(winning_probability, state.surplus[g],
		primitives.gamma[g]) for g in 1:2]
	return merge(state, (bid=bid, K=K,
		winning_probability=winning_probability, value=value))
end

function pooled_log_value(boundary, primitives)
	weight = primitives.class_probability ./ sum(primitives.class_probability)
	return dot(weight, log.(boundary.value))
end

function regularized_boundary(primitives, format, lambda, map, kappa,
		prediction; maxeval=300)
	value_tolerance = parse(Float64,
		get(ENV, "CF_COMMON_SUPPORT_VALUE_TOLERANCE", "1e-4"))
	regularization = parse(Float64,
		get(ENV, "CF_COMMON_SUPPORT_REGULARIZATION", "1e-4"))
	bid_scale = parse(Float64,
		get(ENV, "CF_COMMON_SUPPORT_BID_SCALE", "2.0"))
	log_K_scale = parse(Float64,
		get(ENV, "CF_COMMON_SUPPORT_LOGK_SCALE", "0.5"))
	best = Ref{Any}(nothing)
	function evaluate(parameter)
		boundary = boundary_values(parameter[1], exp(parameter[2]),
			primitives, format, lambda, map)
		value_residual = pooled_log_value(boundary, primitives) - log(kappa)
		distance = ((parameter[1] - prediction[1]) / bid_scale)^2 +
			((parameter[2] - prediction[2]) / log_K_scale)^2
		objective = (value_residual / value_tolerance)^2 +
			regularization * distance
		candidate = merge(boundary, (parameter=copy(parameter),
			value_residual=value_residual,
			continuation_distance=distance, objective=objective))
		if isnothing(best[]) || objective < best[].objective
			best[] = candidate
		end
		return objective
	end
	for bid_offset in (-2.0, 0.0, 2.0), log_K_offset in (-1.0, 0.0, 1.0)
		opt = Opt(:LN_BOBYQA, 2)
		lower_bounds!(opt, [-40.0, log(1e-10)])
		upper_bounds!(opt, [200.0, log(1e12)])
		initial_step!(opt, [0.5, 0.15])
		maxeval!(opt, maxeval)
		xtol_rel!(opt, 1e-8)
		min_objective!(opt, (parameter, _) -> evaluate(parameter))
		try
			optimize(opt, prediction .+ [bid_offset, log_K_offset])
		catch
		end
	end
	isnothing(best[]) && error("Common-support boundary solve failed.")
	return best[]
end

function project_log_K(probability, bid, log_K, primitives, format, lambda;
		tolerance=1e-10, maximum_iterations=20)
	value = clamp(log_K, log(1e-12), log(1e14))
	for _ in 1:maximum_iterations
		K = exp(value)
		state = inverse_state(bid, K, primitives, format, lambda)
		residual = state.probability - probability
		abs(residual) <= tolerance && return value, state, residual
		derivative = state.probability_K * K
		abs(derivative) > 1e-14 || break
		value -= clamp(residual / derivative, -1.0, 1.0)
		value = clamp(value, log(1e-12), log(1e14))
	end
	state = inverse_state(bid, exp(value), primitives, format, lambda)
	return value, state, state.probability - probability
end

function probability_parameter_derivative(probability, bid, log_K,
		primitives, format, lambda, map)
	projected_log_K, state, residual = project_log_K(
		probability, bid, log_K, primitives, format, lambda)
	abs(residual) <= 1e-7 || return nothing
	K = exp(projected_log_K)
	W, dW = map.value_derivative(probability)
	bid_derivative = -dW * K / W
	abs(state.probability_K) > 1e-14 || return nothing
	K_derivative = (1 - state.probability_bid * bid_derivative) /
		state.probability_K
	log_K_derivative = K_derivative / K
	all(isfinite, (bid_derivative, log_K_derivative)) || return nothing
	bid_derivative > 0 || return nothing
	return (derivative=(bid_derivative, log_K_derivative),
		log_K=projected_log_K, state=state)
end

function integrate_strategy(boundary, primitives, format, lambda, map;
		maximum_bid_step=0.10, maximum_probability_step=0.002,
		minimum_probability=0.0005)
	probability = Float64[boundary.probability]
	bid = Float64[boundary.bid]
	log_K = Float64[log(boundary.K)]
	state = Any[boundary]
	step_count = 0
	while probability[end] > minimum_probability + 1e-12
		step_count += 1
		step_count <= 100_000 || error(
			"Probability-parameterized ODE exceeded its step limit.")
		current = probability_parameter_derivative(
			probability[end], bid[end], log_K[end], primitives,
			format, lambda, map)
		isnothing(current) && error(
			"Probability-parameterized ODE failed at an accepted node.")
		log_K[end] = current.log_K
		maximum_probability_change = min(maximum_probability_step,
			maximum_bid_step / max(current.derivative[1], 1e-12),
			probability[end] - minimum_probability)
		step = -maximum_probability_change
		accepted = false
		for _ in 1:20
			P0, b0, L0 = probability[end], bid[end], log_K[end]
			k1 = current.derivative
			s2 = probability_parameter_derivative(
				P0 + step / 2, b0 + step * k1[1] / 2,
				L0 + step * k1[2] / 2, primitives, format, lambda, map)
			if isnothing(s2)
				step /= 2
				continue
			end
			k2 = s2.derivative
			s3 = probability_parameter_derivative(
				P0 + step / 2, b0 + step * k2[1] / 2,
				L0 + step * k2[2] / 2, primitives, format, lambda, map)
			if isnothing(s3)
				step /= 2
				continue
			end
			k3 = s3.derivative
			s4 = probability_parameter_derivative(
				P0 + step, b0 + step * k3[1],
				L0 + step * k3[2], primitives, format, lambda, map)
			if isnothing(s4)
				step /= 2
				continue
			end
			k4 = s4.derivative
			new_probability = P0 + step
			new_bid = b0 + step *
				(k1[1] + 2k2[1] + 2k3[1] + k4[1]) / 6
			new_log_K_guess = L0 + step *
				(k1[2] + 2k2[2] + 2k3[2] + k4[2]) / 6
			new_log_K, new_inverse_state, residual = project_log_K(
				new_probability, new_bid, new_log_K_guess,
				primitives, format, lambda)
			valid = isfinite(new_bid) && isfinite(new_log_K) &&
				new_bid < b0 && b0 - new_bid <= 1.5maximum_bid_step &&
				abs(residual) <= 1e-7
			if !valid
				step /= 2
				continue
			end
			new_state = boundary_values(new_bid, exp(new_log_K),
				primitives, format, lambda, map)
			push!(probability, new_probability)
			push!(bid, new_bid)
			push!(log_K, new_log_K)
			push!(state, new_state)
			accepted = true
			break
		end
		accepted || error(
			"Probability-parameterized ODE could not take a valid step.")
	end
	order = reverse(eachindex(bid))
	return (probability=probability[order], bid=bid[order],
		K=exp.(log_K[order]), state=state[order])
end

function selection_probability(draws, primitives, boundary)
	class = primitives.class
	return [cdf(Normal(), (boundary.x[class[j]] - primitives.x[j]) /
		primitives.distribution[class[j]].bandwidth)
		for j in eachindex(class)]
end

function strategy_inversion(draws, primitives, boundary, path, format, lambda)
	bid = zeros(length(draws.gamma))
	for g in 1:2
		nodes = [state.x[g] for state in path.state]
		order = sortperm(nodes)
		nodes = -decreasing_pava(-nodes[order], ones(length(nodes)))
		bid_nodes = path.bid[order]
		for j in eachindex(bid)
			primitives.class[j] == g || continue
			bid[j] = interpolate_flat(nodes, bid_nodes, primitives.x[j])
		end
	end
	q = [format_share_derivative(
		bid[j], draws.gamma[j], format, lambda)[1] for j in eachindex(bid)]
	risk = 0.5 .* draws.gamma .* (1 .- q).^2 .* environment.sigma2_r
	revenue = format == :share_choice ?
		q .* environment.delta_tilde .* bid .+ (1 .- q) .* environment.mu_r :
		environment.delta_tilde .* bid .+ (1 .- q) .* environment.mu_r
	markup = revenue .- risk .- draws.cost
	return (bid=bid, q=q, cost=draws.cost, markup=markup,
		win_probability=fill(NaN, length(bid)),
		valid=trues(length(bid)))
end

"""
Hold the representative marginal-participant bid fixed and solve its common
K value so the implied participant share equals the observed baseline share.
This avoids extrapolating the empirical winning-probability slope to an
unobserved bid.
"""
function entry_boundary_at_bid(bid, primitives, map;
		target_entry=environment.entry_rate)
	function evaluate(log_K)
		boundary = boundary_values(
			bid, exp(log_K), primitives, :share_choice, NaN, map)
		return boundary.probability - target_entry, boundary
	end
	lower = log(1e-12)
	upper = log(1e12)
	fl, _ = evaluate(lower)
	fu, _ = evaluate(upper)
	(isfinite(fl) && fl >= 0 && isfinite(fu) && fu <= 0) || error(
		"The marginal-participant bid does not bracket the entry-rate target.")
	for _ in 1:100
		middle = (lower + upper) / 2
		fm, _ = evaluate(middle)
		abs(fm) <= 1e-12 && begin
			_, boundary = evaluate(middle)
			return merge(boundary, (
				parameter=[bid, middle], value_residual=0.0,
				continuation_distance=0.0, objective=0.0))
		end
		if fm >= 0
			lower, fl = middle, fm
		else
			upper, fu = middle, fm
		end
	end
	log_K = (lower + upper) / 2
	residual, boundary = evaluate(log_K)
	return merge(boundary, (
		parameter=[bid, log_K], value_residual=0.0,
		continuation_distance=0.0, objective=residual^2))
end

function calibrate_registered_pool(draws, map, bid_anchor;
		alpha_start=nothing)
	components = registered_tilt_components(draws)
	dimension = 5
	best = Ref{Any}(nothing)
	function evaluate(alpha)
		weight = exponential_tilt_weight(components.basis, alpha)
		primitives = primitive_distribution(draws, weight)
		boundary = try
			entry_boundary_at_bid(bid_anchor, primitives, map)
		catch
			return 1e10 + sum(abs2, alpha)
		end
		selection = selection_probability(draws, primitives, boundary)
		entry_weight = weight .* selection
		entry_rate = sum(entry_weight)
		moment_residual = selected_moments(components, entry_weight) .-
			components.target
		entry_residual = entry_rate - environment.entry_rate
		objective = sum(moment_residual.^2) + 1e-8 * sum(abs2, alpha)
		candidate = (alpha=copy(alpha), potential_weight=weight,
			primitives=primitives, boundary=boundary,
			selection=selection, entry_weight=entry_weight,
			entry_rate=entry_rate, moment_residual=moment_residual,
			entry_residual=entry_residual, objective=objective)
		if isnothing(best[]) || objective < best[].objective
			best[] = candidate
		end
		return objective
	end
	opt = Opt(:LN_NELDERMEAD, dimension)
	lower_bounds!(opt, fill(-8.0, dimension))
	upper_bounds!(opt, fill(8.0, dimension))
	maxeval!(opt, controls.registered_maxeval)
	xtol_rel!(opt, 1e-7)
	min_objective!(opt, (alpha, _) -> evaluate(alpha))
	start = zeros(dimension)
	if !isnothing(alpha_start)
		copy_count = min(length(alpha_start), dimension)
		start[1:copy_count] .= Float64.(alpha_start[1:copy_count])
	end
	try
		optimize(opt, start)
	catch error_message
		isnothing(best[]) && rethrow(error_message)
	end
	result = best[]
	isnothing(result) && error("Registered-pool calibration failed.")
	calibrated_kappa = exp(pooled_log_value(
		result.boundary, result.primitives))
	return merge(result, (kappa=calibrated_kappa,
		K_anchor=result.boundary.K, bid_anchor=bid_anchor,
		moment_rmse=sqrt(mean(result.moment_residual.^2)),
		moment_maximum=maximum(abs.(result.moment_residual))))
end
