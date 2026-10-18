# Shared data, numerical, entry, and allocation utilities for the production
# two-risk-type common-support counterfactual.

using CSV
using DataFrames
using Distributions
using LinearAlgebra
using NLopt
using Random
using Statistics

Base.@kwdef struct CFEnvironment
	n_registered::Int
	n_participants::Int
	entry_rate::Float64
	objective_capacity::Float64
	objective_capacity_draws::Vector{Float64} = Float64[]
	objective_log_mean::Float64 = NaN
	objective_log_sd::Float64 = 0.0
	q_lower::Float64
	quantity_minimum_unit::Float64 = 0.0
	delta_tilde::Float64
	mu_r::Float64
	mu_r_share::Float64
	sigma2_r::Float64
end

Base.@kwdef struct CFControls
	n_aux::Int
	n_winprob::Int
	n_alloc::Int
	h_capacity::Float64
	registered_maxeval::Int = 2000
	seed::Int = 20260916
end

struct CFDraws
	gamma::Vector{Float64}
	cost::Vector{Float64}
	capacity::Vector{Float64}
	z_gamma::Vector{Float64}
end

function objective_capacity_draws(environment::CFEnvironment, n::Int,
		rng::AbstractRNG)
	if !isempty(environment.objective_capacity_draws)
		return environment.objective_capacity_draws[
			rand(rng, 1:length(environment.objective_capacity_draws), n)]
	elseif isfinite(environment.objective_log_mean)
		return exp.(environment.objective_log_mean .+
			environment.objective_log_sd .* randn(rng, n))
	end
	return fill(environment.objective_capacity, n)
end

function weighted_indices(rng::AbstractRNG, weights::AbstractVector, n::Int)
	probability = cumsum(weights ./ sum(weights))
	return searchsortedfirst.(Ref(probability), rand(rng, n))
end

function systematic_weighted_indices(rng::AbstractRNG,
		weights::AbstractVector, n::Int)
	probability = cumsum(weights ./ sum(weights))
	u = ((0:(n - 1)) .+ rand(rng)) ./ n
	index = searchsortedfirst.(Ref(probability), u)
	shuffle!(rng, index)
	return index
end

wholesale_outside_floor(gamma, mu_r, sigma2_r) =
	mu_r .- 0.5 .* gamma .* sigma2_r

function apply_wholesale_outside_option(cost, gamma, cost_shift,
		environment::CFEnvironment)
	cost_shift <= 0 || error(
		"Positive shifts require an effective-cost decomposition.")
	floor = wholesale_outside_floor(
		gamma, environment.mu_r, environment.sigma2_r)
	unfloored = cost .+ cost_shift
	return (cost=max.(unfloored, floor), floor=floor,
		active=unfloored .<= floor .+ 1e-10, unfloored=unfloored)
end

function standardized_draws(rng::AbstractRNG, n::Int)
	draw = randn(rng, n)
	return (draw .- mean(draw)) ./ std(draw)
end

function weighted_standard_deviation(x, weight)
	normalized_weight = weight ./ sum(weight)
	center = sum(normalized_weight .* x)
	return sqrt(sum(normalized_weight .* (x .- center).^2))
end

function participation_statistics(inversion, draws::CFDraws,
		environment::CFEnvironment; entry_weight, potential_weight,
		low_gamma_cutoff)
	entry_rate = sum(entry_weight)
	potential_mass = sum(potential_weight)
	q = inversion.q
	entrant_mean_gamma = sum(entry_weight .* draws.gamma) / entry_rate
	potential_mean_gamma = sum(potential_weight .* draws.gamma) / potential_mass
	low_gamma = Float64.(draws.gamma .<= low_gamma_cutoff)
	entrant_low_gamma_share = sum(entry_weight .* low_gamma) / entry_rate
	potential_low_gamma_share = sum(potential_weight .* low_gamma) / potential_mass
	risk_premium = 0.5 .* draws.gamma .* (1 .- q).^2 .* environment.sigma2_r
	valid_weight = entry_weight .* Float64.(inversion.valid .&
		isfinite.(inversion.markup))
	valid_mass = sum(valid_weight)
	return (
		entrant_mean_q=sum(entry_weight .* q) / entry_rate,
		potential_project_mean_q=sum(entry_weight .* q) / potential_mass,
		entrant_share_q_zero=sum(entry_weight .* Float64.(q .<= 1e-10)) /
			entry_rate,
		entrant_share_q_one=sum(entry_weight .* Float64.(q .>= 1 - 1e-10)) /
			entry_rate,
		entrant_mean_gamma=entrant_mean_gamma,
		potential_mean_gamma=potential_mean_gamma,
		entrant_gamma_selection_ratio=entrant_mean_gamma / potential_mean_gamma,
		entrant_low_gamma_share=entrant_low_gamma_share,
		potential_low_gamma_share=potential_low_gamma_share,
		entrant_low_gamma_selection_ratio=entrant_low_gamma_share /
			potential_low_gamma_share,
		entrant_mean_risk_premium=sum(entry_weight .* risk_premium) / entry_rate,
		entrant_mean_bid=sum(entry_weight .* inversion.bid) / entry_rate,
		entrant_bid_sd=weighted_standard_deviation(inversion.bid, entry_weight),
		entrant_mean_markup=sum(valid_weight .* ifelse.(inversion.valid,
			inversion.markup, 0.0)) / valid_mass,
		valid_entry_share=valid_mass / entry_rate)
end

function lottery_ce(probability, surplus, gamma)
	probability = clamp(probability, 0.0, 1.0)
	return -log1p(probability * expm1(-gamma * surplus)) / gamma
end

function decreasing_pava(value, weight)
	block_value = Float64[]
	block_weight = Float64[]
	block_start = Int[]
	block_end = Int[]
	for j in eachindex(value)
		push!(block_value, value[j]); push!(block_weight, weight[j])
		push!(block_start, j); push!(block_end, j)
		while length(block_value) >= 2 && block_value[end - 1] < block_value[end]
			combined_weight = block_weight[end - 1] + block_weight[end]
			block_value[end - 1] = (
				block_weight[end - 1] * block_value[end - 1] +
				block_weight[end] * block_value[end]) / combined_weight
			block_weight[end - 1] = combined_weight
			block_end[end - 1] = block_end[end]
			pop!(block_value); pop!(block_weight)
			pop!(block_start); pop!(block_end)
		end
	end
	fitted = similar(Float64.(value))
	for k in eachindex(block_value)
		fitted[block_start[k]:block_end[k]] .= block_value[k]
	end
	return fitted
end

function registered_tilt_components(draws::CFDraws)
	z_loggamma = (log.(draws.gamma) .- mean(log.(draws.gamma))) ./
		std(log.(draws.gamma))
	z_cost = (draws.cost .- mean(draws.cost)) ./ std(draws.cost)
	z_gamma_level = (draws.gamma .- mean(draws.gamma)) ./ std(draws.gamma)
	tail = max.(z_cost, 0.0).^3
	tail = (tail .- mean(tail)) ./ std(tail)
	basis = hcat(z_loggamma, z_cost, z_loggamma.^2 .- 1,
		z_cost.^2 .- 1, z_loggamma .* z_cost, tail)
	target = [0.0, 1.0, 0.0, 1.0, mean(z_gamma_level .* z_cost)]
	return (basis=basis, z_loggamma=z_loggamma, z_cost=z_cost,
		z_gamma_level=z_gamma_level, target=target)
end

function exponential_tilt_weight(basis, alpha)
	log_weight = basis[:, eachindex(alpha)] * alpha
	log_weight .-= maximum(log_weight)
	weight = exp.(log_weight)
	return weight ./ sum(weight)
end

function selected_moments(components, entrant_weight)
	weight = entrant_weight ./ sum(entrant_weight)
	return [sum(weight .* components.z_loggamma),
		sum(weight .* components.z_loggamma.^2),
		sum(weight .* components.z_cost),
		sum(weight .* components.z_cost.^2),
		sum(weight .* components.z_gamma_level .* components.z_cost)]
end

function full_award_winners(order, capacity, objective_capacity)
	cumulative = cumsum(capacity[order])
	winner_count = min(
		searchsortedfirst(cumulative, objective_capacity), length(order))
	return order[1:winner_count]
end

function allocation_statistics(solution, format, lambda, draws::CFDraws,
		environment::CFEnvironment, controls::CFControls;
		entry_weight=ones(length(draws.gamma)),
		potential_weight=fill(inv(length(draws.gamma)), length(draws.gamma)),
		outside_active=falses(length(draws.gamma)),
		inversion=solution.details.implied,
		n_participants=environment.n_participants,
		low_gamma_cutoff=median(draws.gamma))
	rng = MersenneTwister(controls.seed + 701)
	valid_weight = entry_weight .* Float64.(inversion.valid)
	totals = Dict(name => 0.0 for name in (
		:net_cost, :variance, :rent, :cost_loss, :winner_count,
		:awarded_capacity, :q, :winner_gamma, :winner_cost, :winner_bid,
		:winner_risk_premium, :winner_low_gamma, :winner_q_zero,
		:winner_q_one, :outside_capacity, :loser_count, :loser_gamma,
		:loser_cost, :loser_low_gamma, :loser_outside_count,
		:loser_outside_gamma, :loser_outside_cost,
		:loser_outside_low_gamma, :loser_outside_capacity,
		:overshoot, :relative_overshoot))
	objective_met_count = 0
	insufficient_capacity_count = 0
	objective_draws = objective_capacity_draws(environment, controls.n_alloc, rng)
	potential_probability = potential_weight ./ sum(potential_weight)
	entry_probability = entry_weight ./ sum(potential_weight)
	nonentrant_probability = potential_probability .- entry_probability
	minimum(nonentrant_probability) >= -1e-12 || error(
		"Entry weights cannot exceed potential weights.")
	nonentrant_probability = max.(nonentrant_probability, 0.0)
	low_gamma = Float64.(draws.gamma .<= low_gamma_cutoff)
	nonentrant_count = environment.n_registered * sum(nonentrant_probability)
	nonentrant_gamma_sum = environment.n_registered * sum(
		nonentrant_probability .* draws.gamma)
	nonentrant_cost_sum = environment.n_registered * sum(
		nonentrant_probability .* draws.cost)
	nonentrant_low_gamma_count = environment.n_registered * sum(
		nonentrant_probability .* low_gamma)
	nonentrant_outside_probability = nonentrant_probability .*
		Float64.(outside_active)
	nonentrant_outside_count = environment.n_registered * sum(
		nonentrant_outside_probability)
	nonentrant_outside_gamma_sum = environment.n_registered * sum(
		nonentrant_outside_probability .* draws.gamma)
	nonentrant_outside_cost_sum = environment.n_registered * sum(
		nonentrant_outside_probability .* draws.cost)
	nonentrant_outside_low_gamma_count = environment.n_registered * sum(
		nonentrant_outside_probability .* low_gamma)
	nonentrant_outside_capacity = environment.n_registered * sum(
		nonentrant_outside_probability .* draws.capacity)
	sampled_participant_count = round(Int, n_participants)
	participant_count_scale = sampled_participant_count == 0 ? 0.0 :
		n_participants / sampled_participant_count
	for allocation_index in 1:controls.n_alloc
		objective_capacity = objective_draws[allocation_index]
		index = weighted_indices(rng, valid_weight, sampled_participant_count)
		bid = inversion.bid[index]
		capacity = draws.capacity[index]
		order = sortperm(bid)
		winner = full_award_winners(order, capacity, objective_capacity)
		q = inversion.q[index]
		cost = inversion.cost[index]
		gamma = draws.gamma[index]
		markup = inversion.markup[index]
		policy_cost = if format == :share_choice
			q .* (environment.delta_tilde .* bid .- environment.mu_r)
		elseif format == :subsidy
			environment.delta_tilde .* bid
		else
			environment.delta_tilde .* bid .- q .* environment.mu_r
		end
		totals[:net_cost] += sum(capacity[winner] .* policy_cost[winner]) /
			objective_capacity
		exposure = sum(capacity[winner] .* q[winner]) / objective_capacity
		totals[:variance] += exposure^2 * environment.sigma2_r
		totals[:rent] += sum(capacity[winner] .* markup[winner]) /
			objective_capacity
		first_best = full_award_winners(sortperm(cost), capacity,
			objective_capacity)
		totals[:cost_loss] += (sum(capacity[winner] .* cost[winner]) -
			sum(capacity[first_best] .* cost[first_best])) / objective_capacity
		totals[:winner_count] += length(winner)
		awarded_capacity = sum(capacity[winner])
		totals[:awarded_capacity] += awarded_capacity
		winner_indicator = falses(length(index)); winner_indicator[winner] .= true
		loser = .!winner_indicator
		loser_outside = loser .& outside_active[index]
		totals[:loser_count] += sum(loser)
		totals[:loser_gamma] += sum(gamma[loser])
		totals[:loser_cost] += sum(cost[loser])
		totals[:loser_low_gamma] += sum(Float64.(
			gamma[loser] .<= low_gamma_cutoff))
		totals[:loser_outside_count] += sum(loser_outside)
		totals[:loser_outside_gamma] += sum(gamma[loser_outside])
		totals[:loser_outside_cost] += sum(cost[loser_outside])
		totals[:loser_outside_low_gamma] += sum(Float64.(
			gamma[loser_outside] .<= low_gamma_cutoff))
		loser_outside_capacity = sum(capacity[loser_outside])
		totals[:loser_outside_capacity] += loser_outside_capacity
		totals[:outside_capacity] +=
			nonentrant_outside_capacity + loser_outside_capacity
		if awarded_capacity >= objective_capacity
			overshoot = awarded_capacity - objective_capacity
			totals[:overshoot] += overshoot
			totals[:relative_overshoot] += overshoot / objective_capacity
			objective_met_count += 1
		else
			insufficient_capacity_count += 1
		end
		winner_weight = capacity[winner]
		winner_mass = sum(winner_weight)
		totals[:q] += sum(winner_weight .* q[winner]) / winner_mass
		totals[:winner_gamma] += sum(winner_weight .* gamma[winner]) / winner_mass
		totals[:winner_cost] += sum(winner_weight .* cost[winner]) / winner_mass
		totals[:winner_bid] += sum(winner_weight .* bid[winner]) / winner_mass
		risk_premium = 0.5 .* gamma .* (1 .- q).^2 .* environment.sigma2_r
		totals[:winner_risk_premium] += sum(
			winner_weight .* risk_premium[winner]) / winner_mass
		totals[:winner_low_gamma] += sum(winner_weight .* Float64.(
			gamma[winner] .<= low_gamma_cutoff)) / winner_mass
		totals[:winner_q_zero] += sum(winner_weight .* Float64.(
			q[winner] .<= 1e-10)) / winner_mass
		totals[:winner_q_one] += sum(winner_weight .* Float64.(
			q[winner] .>= 1 - 1e-10)) / winner_mass
	end
	n = controls.n_alloc
	auction_capacity = totals[:awarded_capacity] / n
	outside_capacity = totals[:outside_capacity] / n
	loser_count = totals[:loser_count] / n * participant_count_scale
	loser_outside_count = totals[:loser_outside_count] / n *
		participant_count_scale
	loser_outside_capacity = totals[:loser_outside_capacity] / n
	outside_builder_count = nonentrant_outside_count + loser_outside_count
	outside_builder_gamma_sum = nonentrant_outside_gamma_sum +
		totals[:loser_outside_gamma] / n * participant_count_scale
	outside_builder_cost_sum = nonentrant_outside_cost_sum +
		totals[:loser_outside_cost] / n * participant_count_scale
	outside_builder_low_gamma_count = nonentrant_outside_low_gamma_count +
		totals[:loser_outside_low_gamma] / n * participant_count_scale
	total_capacity = auction_capacity + outside_capacity
	safe_ratio(numerator, denominator) = denominator <= 0 ? NaN :
		numerator / denominator
	return (
		pm_expected_net_cost=totals[:net_cost] / n,
		pm_variance=totals[:variance] / n,
		winner_information_rent=totals[:rent] / n,
		cost_loss=totals[:cost_loss] / n,
		winner_count=totals[:winner_count] / n,
		awarded_capacity=auction_capacity,
		auction_supported_capacity=auction_capacity,
		outside_built_capacity=outside_capacity,
		nonentrant_outside_built_capacity=nonentrant_outside_capacity,
		loser_outside_built_capacity=loser_outside_capacity,
		total_built_capacity=total_capacity,
		outside_built_capacity_share=total_capacity <= 0 ? NaN :
			outside_capacity / total_capacity,
		nonentrant_project_share=nonentrant_count / environment.n_registered,
		auction_loser_project_share=loser_count / environment.n_registered,
		auction_winner_project_share=(totals[:winner_count] / n *
			participant_count_scale) /
			environment.n_registered,
		nonentrant_outside_build_share=safe_ratio(
			nonentrant_outside_count, nonentrant_count),
		loser_outside_build_share=safe_ratio(
			loser_outside_count, loser_count),
		outside_builder_nonentrant_share=safe_ratio(
			nonentrant_outside_count, outside_builder_count),
		outside_capacity_nonentrant_share=safe_ratio(
			nonentrant_outside_capacity, outside_capacity),
		nonentrant_mean_gamma=safe_ratio(
			nonentrant_gamma_sum, nonentrant_count),
		nonentrant_low_gamma_share=safe_ratio(
			nonentrant_low_gamma_count, nonentrant_count),
		nonentrant_mean_cost=safe_ratio(
			nonentrant_cost_sum, nonentrant_count),
		loser_mean_gamma=safe_ratio(totals[:loser_gamma], totals[:loser_count]),
		loser_low_gamma_share=safe_ratio(
			totals[:loser_low_gamma], totals[:loser_count]),
		loser_mean_cost=safe_ratio(totals[:loser_cost], totals[:loser_count]),
		nonentrant_outside_mean_gamma=safe_ratio(
			nonentrant_outside_gamma_sum, nonentrant_outside_count),
		nonentrant_outside_low_gamma_share=safe_ratio(
			nonentrant_outside_low_gamma_count, nonentrant_outside_count),
		nonentrant_outside_mean_cost=safe_ratio(
			nonentrant_outside_cost_sum, nonentrant_outside_count),
		loser_outside_mean_gamma=safe_ratio(
			totals[:loser_outside_gamma], totals[:loser_outside_count]),
		loser_outside_low_gamma_share=safe_ratio(
			totals[:loser_outside_low_gamma], totals[:loser_outside_count]),
		loser_outside_mean_cost=safe_ratio(
			totals[:loser_outside_cost], totals[:loser_outside_count]),
		outside_builder_mean_gamma=safe_ratio(
			outside_builder_gamma_sum, outside_builder_count),
		outside_builder_low_gamma_share=safe_ratio(
			outside_builder_low_gamma_count, outside_builder_count),
		outside_builder_mean_cost=safe_ratio(
			outside_builder_cost_sum, outside_builder_count),
		winner_mean_capacity=totals[:awarded_capacity] /
			totals[:winner_count],
		winner_mean_q=totals[:q] / n,
		winner_mean_gamma=totals[:winner_gamma] / n,
		winner_mean_cost=totals[:winner_cost] / n,
		winner_mean_bid=totals[:winner_bid] / n,
		winner_mean_risk_premium=totals[:winner_risk_premium] / n,
		winner_low_gamma_share=totals[:winner_low_gamma] / n,
		winner_share_q_zero=totals[:winner_q_zero] / n,
		winner_share_q_one=totals[:winner_q_one] / n,
		mean_capacity_overshoot=objective_met_count == 0 ? NaN :
			totals[:overshoot] / objective_met_count,
		mean_relative_capacity_overshoot=objective_met_count == 0 ? NaN :
			totals[:relative_overshoot] / objective_met_count,
		insufficient_capacity_probability=insufficient_capacity_count / n)
end
