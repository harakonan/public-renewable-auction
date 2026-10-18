#=
Preparation for synthetic counterfactual auctions.

Large simulated objects are written to the timestamp-free
Dropbox/renewable-auction/intermediate/counterfactual_task directory. Full
aggregated outputs are written separately to counterfactual_output.
=#

include("counterfactual_functions.jl")

using Printf

cd(dirname(dirname(@__DIR__)))

path_to_data = strip(read("data.txt", String))
intermediate_directory = joinpath(
	path_to_data, "intermediate", "counterfactual_task")
mkpath(intermediate_directory)

primitive_table = CSV.read(joinpath(
	intermediate_directory, "participant_primitives.csv"), DataFrame)
gamma_cost_distribution = CSV.read(joinpath(
	"output", "result", "gamma_c_dist_result.csv"), DataFrame)

function distribution_estimate(variable)
	value = gamma_cost_distribution[
		gamma_cost_distribution.variable .== variable, :estimate]
	length(value) == 1 || error("Distribution statistic $variable is not unique.")
	return Float64(only(value))
end

# Table 5 reports total-capacity-weighted averages of within-auction moments.
# These, rather than moments of the raw mixture across auctions, define the
# representative primitive distribution used in the counterfactual.

representative_primitive = (
	loggamma_mean=distribution_estimate("loggamma_mean"),
	loggamma_var=distribution_estimate("loggamma_var"),
	gamma_mean=distribution_estimate("gamma_mean"),
	gamma_var=distribution_estimate("gamma_var"),
	cost_mean=distribution_estimate("c_mean"),
	cost_var=distribution_estimate("c_var"),
	factual_rho=distribution_estimate("gamma_c_cor"))
rho_sensitivity_lower = parse(Float64,
	get(ENV, "CF_RHO_LOWER", "-0.658"))
rho_sensitivity_upper = parse(Float64,
	get(ENV, "CF_RHO_UPPER", "0.454"))
auction_data = CSV.read(joinpath(
	path_to_data, "processed", "generation_auction_wind20.csv"), DataFrame)
auction_data = auction_data[(auction_data.partial_bid_quantity .== 0) .&
	(auction_data.descending_clock .== 0), :]

active_auction_ids = sort(unique(primitive_table.auction_id))
auction_data = auction_data[in.(auction_data.auction_id,
	Ref(Set(active_auction_ids))), :]
auction_rows = unique(auction_data[:, [
	:auction_id, :year_decimal, :lead_time, :number_registered,
	:number_participant, :total_capacity_registered, :total_capacity,
	:procurement_quantity,
	:bid_share_lb_auction, :quantity_minimum_unit, :mc_operation]])

participant_count = Dict(
	row.auction_id => Float64(row.number_participant)
	for row in eachrow(auction_rows))
primitive_weight = [participant_count[id] for id in primitive_table.auction_id]

params1 = CSV.read(joinpath("output", "result", "params1_hat.csv"), DataFrame)
params2 = CSV.read(joinpath(
	"output", "result", "params2_dist_result.csv"), DataFrame)
params_specified = CSV.read(joinpath(
	"output", "result", "params_specified.csv"), DataFrame)
alpha_r0 = only(params1[params1.variable .== "alpha_r0", :hetero])
alpha_r1 = only(params1[params1.variable .== "alpha_r1", :hetero])
spot_result = CSV.read(joinpath("output", "result", "spot_result.csv"), DataFrame)
spot_variance = spot_result[occursin.(r"var_r_a", spot_result.variable), :est]

row_weight = Float64.(auction_rows.number_participant)
weighted_mean(x) = sum(row_weight .* Float64.(x)) / sum(row_weight)

delta = 0.95
contract_length = 20
delta_tilde = [sum(delta^t for t in row.lead_time:1.0:
	(row.lead_time + contract_length - 1)) / contract_length
	for row in eachrow(auction_rows)]
mu_r = alpha_r0 + alpha_r1 * weighted_mean(auction_rows.mc_operation)
mean_lead_discount = weighted_mean(delta .^ Float64.(auction_rows.lead_time))
sigma2_r = sum(row_weight .* [spot_variance[row.auction_id]
	for row in eachrow(auction_rows)]) / sum(row_weight)

entry_rate = sum(Float64.(auction_rows.number_participant)) /
	sum(Float64.(auction_rows.number_registered))
registrant_mean_nameplate_capacity =
	sum(Float64.(auction_rows.total_capacity_registered)) /
	sum(Float64.(auction_rows.number_registered))
participant_mean_nameplate_capacity =
	sum(Float64.(auction_rows.total_capacity)) /
	sum(Float64.(auction_rows.number_participant))
observed_capacity_selection_ratio = participant_mean_nameplate_capacity /
	registrant_mean_nameplate_capacity

mode = Symbol(get(ENV, "CF_MODE", "small"))
n_registered_full = round(Int, median(Float64.(auction_rows.number_registered)))
n_registered = mode == :small ? min(n_registered_full,
	parse(Int, get(ENV, "CF_N_REGISTERED", "120"))) : n_registered_full
controls = CFControls(
	n_aux=parse(Int, get(ENV, "CF_N_AUX", mode == :small ? "800" : "5000")),
	n_winprob=parse(Int, get(
		ENV, "CF_N_WINPROB", mode == :small ? "80" : "1000")),
	n_alloc=parse(Int, get(
		ENV, "CF_N_ALLOC", mode == :small ? "200" : "2000")),
	h_capacity=parse(Float64, get(ENV, "CF_H_CAPACITY", "15.0")),
	registered_maxeval=parse(Int, get(ENV, "CF_REGISTERED_MAXEVAL",
		mode == :small ? "800" : "2000")))
controls.n_aux >= 100 || error("CF_N_AUX must be at least 100.")
controls.n_winprob >= 50 || error("CF_N_WINPROB must be at least 50.")
controls.n_alloc >= 1 || error("CF_N_ALLOC must be positive.")
controls.registered_maxeval >= 100 ||
	error("CF_REGISTERED_MAXEVAL must be at least 100.")

n_participants = max(2, round(Int, n_registered * entry_rate))
mean_capacity = sum(primitive_weight .* primitive_table.capacity) /
	sum(primitive_weight)

function parameter_value(table, variable, column)
	value = table[table.variable .== variable, column]
	length(value) == 1 || error("Parameter $variable is not unique.")
	return Float64(only(value))
end

"""
Construct the representative baseline marginal-participant bid anchor.

For each baseline sealed-bid auction, approximate the marginal participant's
bid by the simulated expectation of the highest of N draws from the fitted
participant bid distribution. The representative anchor is the
participant-count-weighted mean of the auction-specific marginal-participant
bids.
"""
function marginal_participant_bid_anchor()
	n_simulation = parse(Int, get(
		ENV, "CF_MARGINAL_BID_DRAWS", "5000"))
	n_simulation >= 500 || error(
		"CF_MARGINAL_BID_DRAWS must be at least 500.")
	sigma_b = sqrt(parameter_value(params2, "sigma2_b", :estimate))

	diagnostic = DataFrame(
		auction_id=Int[], number_participant=Int[],
		maximum_winning_bid=Float64[], marginal_participant_bid=Float64[])
	for row in eachrow(auction_rows)
		auction_id = Int(row.auction_id)
		n_participant = Int(row.number_participant)
		n_participant >= 2 || continue
		rng = MersenneTwister(controls.seed + 8_000 + auction_id)

		mu_b = parameter_value(
			params2, "beta_b_a$(auction_id)", :estimate)
		bid_draw = mu_b .+ sigma_b .* randn(
			rng, n_simulation, n_participant)
		marginal_bid = mean(vec(maximum(bid_draw; dims=2)))
		maximum_winning_bid = maximum(Float64.(
			auction_data[auction_data.auction_id .== auction_id, :bid_price]))
		push!(diagnostic, (
			auction_id=auction_id,
			number_participant=n_participant,
			maximum_winning_bid=maximum_winning_bid,
			marginal_participant_bid=marginal_bid))
	end

	nrow(diagnostic) > 0 || error(
		"No valid auction-specific marginal-participant bid anchors.")
	weight = Float64.(diagnostic.number_participant)
	anchor = sum(weight .* diagnostic.marginal_participant_bid) / sum(weight)
	return (anchor=anchor, diagnostic=diagnostic)
end

"""
Simulate the objective-capacity cutoff used by the counterfactual auctions.

For a factual auction draw, procurement capacity is filled in units of
q_i * Capacity_i. The marginal project is prorated, which maps that factual
procurement draw into Dtilde measured in full installed capacity. Draws are
then rescaled by expected total installed capacity so the pooled factual
cutoff distribution applies to the representative synthetic auction.
"""
function simulate_objective_capacity_draws(n_draws, representative_n,
		representative_mean_capacity, seed)
	rng = MersenneTwister(seed)
	row_probability = Float64.(auction_rows.number_participant)
	row_probability ./= sum(row_probability)
	row_index = rand(rng, Categorical(row_probability), n_draws)

	beta_C_b = parameter_value(params2, "beta_C_b", :estimate)
	sigma2_C = parameter_value(params2, "sigma2_C", :estimate)
	beta_K_b = parameter_value(params2, "beta_K_b", :estimate)
	beta_K_logC = parameter_value(params2, "beta_K_logC", :estimate)
	sigma2_K = parameter_value(params2, "sigma2_K", :estimate)
	sigma2_b = parameter_value(params2, "sigma2_b", :estimate)
	alpha_gamma0 = parameter_value(params1, "alpha_gamma0", :hetero)
	alpha_gammaK = parameter_value(params1, "alpha_gammaK", :hetero)
	beta_D = [parameter_value(
		params_specified, "beta_D$(j)_sb", :value) for j in 0:3]
	sigma2_D = parameter_value(params_specified, "sigma2_D_sb", :value)

	draws = Vector{Float64}(undef, n_draws)
	for s in eachindex(draws)
		row = auction_rows[row_index[s], :]
		auction_id = Int(row.auction_id)
		n = Int(row.number_participant)
		mu_b = parameter_value(
			params2, "beta_b_a$(auction_id)", :estimate)
		beta_C_a = parameter_value(
			params2, "beta_C_a$(auction_id)", :estimate)
		beta_K_a = parameter_value(
			params2, "beta_K_a$(auction_id)", :estimate)
		bid = mu_b .+ sqrt(sigma2_b) .* randn(rng, n)
		capacity = exp.(beta_C_b .* bid .+ beta_C_a .+
			sqrt(sigma2_C) .* randn(rng, n))
		K_log = isfinite(beta_K_a) ? max.(0.0,
			beta_K_a .+ beta_K_b .* bid .+
			beta_K_logC .* log.(capacity) .+
			sqrt(sigma2_K) .* randn(rng, n)) : zeros(n)
		gamma = exp.(alpha_gamma0 .+ alpha_gammaK .* K_log)
		delta_tilde_a = sum(delta^t for t in row.lead_time:1.0:
			(row.lead_time + contract_length - 1)) / contract_length
		mu_r_a = alpha_r0 + alpha_r1 * row.mc_operation
		q_unconstrained = 1 .- (delta^row.lead_time * mu_r_a .-
			delta_tilde_a .* bid) ./ (gamma .* spot_variance[auction_id])
		unit = Float64(row.quantity_minimum_unit)
		q_lower = ceil.(capacity .* Float64(row.bid_share_lb_auction) ./
			unit .- eps(1.0)) .* unit ./ capacity
		q_upper = floor.(capacity ./ unit .+ eps(1.0)) .* unit ./ capacity
		q = clamp.(q_unconstrained, q_lower, q_upper)

		mu_D = beta_D[1] + beta_D[2] *
			(row.year_decimal + row.lead_time) +
			beta_D[3] * row.mc_operation + beta_D[4] * n
		procurement_capacity = exp(mu_D + sqrt(sigma2_D) * randn(rng))
		order = sortperm(bid)
		contracted = q[order] .* capacity[order]
		cumulative = cumsum(contracted)
		marginal = searchsortedfirst(cumulative, procurement_capacity)
		Dtilde = if marginal > n
			sum(capacity)
		else
			previous_contracted = marginal == 1 ? 0.0 : cumulative[marginal - 1]
			previous_capacity = marginal == 1 ? 0.0 :
				sum(@view capacity[order[1:(marginal - 1)]])
			fraction = (procurement_capacity - previous_contracted) /
				contracted[marginal]
			previous_capacity + clamp(fraction, 0.0, 1.0) *
				capacity[order[marginal]]
		end

		expected_capacity = exp(beta_C_a + beta_C_b * mu_b +
			0.5 * (sigma2_C + beta_C_b^2 * sigma2_b))
		draws[s] = Dtilde / (n * expected_capacity) *
			representative_n * representative_mean_capacity
	end
	return draws
end

n_objective_draws = parse(Int, get(ENV, "CF_N_OBJECTIVE_DRAWS",
	mode == :small ? "500" : "2000"))
objective_draws = simulate_objective_capacity_draws(
	n_objective_draws, n_participants, mean_capacity, controls.seed + 31)
objective_capacity = mean(objective_draws)
objective_share = objective_capacity / (n_participants * mean_capacity)
counterfactual_q_lower = parse(Float64, get(ENV, "CF_Q_LOWER", "0.0"))
0.0 <= counterfactual_q_lower < 1.0 ||
	error("CF_Q_LOWER must lie in [0,1).")

environment = CFEnvironment(
	n_registered=n_registered,
	n_participants=n_participants,
	entry_rate=entry_rate,
	objective_capacity=objective_capacity,
	objective_capacity_draws=objective_draws,
	q_lower=counterfactual_q_lower,
	delta_tilde=sum(row_weight .* delta_tilde) / sum(row_weight),
	mu_r=mu_r,
	mu_r_share=mean_lead_discount * mu_r,
	sigma2_r=sigma2_r)

function parse_grid(name, default)
	return parse.(Float64, split(get(ENV, name, default), ","))
end

# Hold the level mean of risk aversion fixed as dispersion changes.  Defining
# the target from the primary lognormal moments preserves the primary
# distribution exactly even though separately averaged level and log moments
# need not satisfy the lognormal identity after aggregation across auctions.
counterfactual_gamma_mean = exp(
	representative_primitive.loggamma_mean +
	0.5 * representative_primitive.loggamma_var)
