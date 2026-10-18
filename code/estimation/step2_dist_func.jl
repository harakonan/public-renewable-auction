#=
Brazil new energy auction structural model estimation
sealed-bid auction concomitant and bid-price distributions
functions

@author: Konan Hara

@date: September 13, 2026

@software: Julia version 1.12

=#

# Weighted Gaussian model for log C conditional on bid and auction. Common
# slopes are identified from within-auction variation; unsupported auction
# effects remain NaN rather than reverting to original-sample estimates.
function weighted_capacity_distribution(sb_mat, weights)
	b = Float64.(sb_mat[:, 2])
	logC = log.(Float64.(sb_mat[:, 3]))
	active = [sum(weights[a_sb .== a]) > 0 for a in a_sb_list]

	b_mean = fill(NaN, A_sb)
	logC_mean = fill(NaN, A_sb)
	for a in a_sb_list[active]
		mask = a_sb .== a
		weight_a = sum(weights[mask])
		b_mean[a] = sum(weights[mask] .* b[mask]) / weight_a
		logC_mean[a] = sum(weights[mask] .* logC[mask]) / weight_a
	end

	numerator = 0.0
	denominator = 0.0
	for a in a_sb_list[active]
		mask = a_sb .== a
		b_dev = b[mask] .- b_mean[a]
		logC_dev = logC[mask] .- logC_mean[a]
		numerator += sum(weights[mask] .* b_dev .* logC_dev)
		denominator += sum(weights[mask] .* b_dev .^ 2)
	end
	denominator > 0 || error("No within-auction bid variation remains for f(C | b*,a)")
	βC_b = numerator / denominator

	βC_a = fill(NaN, A_sb)
	for a in a_sb_list[active]
		mask = a_sb .== a
		βC_a[a] = sum(weights[mask] .* (logC[mask] .- βC_b .* b[mask])) /
			sum(weights[mask])
	end
	positive_weight = weights .> 0
	residual = logC[positive_weight] .- βC_b .* b[positive_weight] .-
		βC_a[a_sb[positive_weight]]
	σ2C = sum(weights[positive_weight] .* residual .^ 2) /
		sum(weights[positive_weight])
	return [βC_b; βC_a], σ2C, active
end

# Censored-normal model for K_log=log(1+K/100), where K is prior-started wind
# capacity in MW. The observed transformed capacity is max(0,K_log*). Auctions
# before any previously awarded wind contract could have started have K_log=0
# structurally.
function weighted_prior_distribution(sb_mat, weights;
	βK_init=missing, σ2K_start=σ2K_init)
	b = Float64.(sb_mat[:, 2])
	logC = log.(Float64.(sb_mat[:, 3]))
	K = prior_assurance_100_sb
	y = log1p.(K)
	active = [sum(weights[a_sb .== a]) > 0 for a in a_sb_list]
	positive_support = [sum(weights[(a_sb .== a) .& (K .> 0)]) > 0 for a in a_sb_list]
	estimable = active .& .!structural_zero_prior_sb_a .& positive_support
	estimable_auctions = a_sb_list[estimable]
	isempty(estimable_auctions) && error("No positive prior-capacity support remains")

	keep = (weights .> 0) .& estimable[a_sb]
	b_kept = b[keep]
	logC_kept = logC[keep]
	y_kept = y[keep]
	weights_kept = Float64.(weights[keep])
	weight_total = sum(weights_kept)
	position = Dict(a => k for (k, a) in enumerate(estimable_auctions))
	a_index = [position[a] for a in a_sb[keep]]
	zero_index = findall(iszero, y_kept)
	positive_index = findall(>(0), y_kept)

	model = Model(Ipopt.Optimizer)
	set_silent(model)
	set_optimizer_attribute(model, "max_iter", max_iter)
	set_optimizer_attribute(model, "tol", tol)
	set_optimizer_attribute(model, "hessian_approximation", "limited-memory")
	set_optimizer_attribute(model, "acceptable_tol", 1e-5)
	set_optimizer_attribute(model, "acceptable_iter", 5)
	func_logPhi_zero(μ, σ2) = logcdf(Normal(), -μ / sqrt(σ2))
	JuMP.register(model, :func_logPhi_zero, 2, func_logPhi_zero, autodiff=true)

	if ismissing(βK_init)
		@variable(model, βK_b)
		@variable(model, βK_logC)
		@variable(model, βK_a[1:length(estimable_auctions)])
	else
		@variable(model, βK_b, start=βK_init[1])
		@variable(model, βK_logC, start=βK_init[2])
		@variable(model, βK_a[k=1:length(estimable_auctions)],
			start=βK_init[2 + estimable_auctions[k]])
	end
	@variable(model, σ2K_lb <= σ2K <= σ2K_ub, start=σ2K_start)
	@NLexpression(model, μK[i=eachindex(y_kept)],
		βK_a[a_index[i]] + βK_b * b_kept[i] + βK_logC * logC_kept[i])
	@NLobjective(model, Min,
		(sum(weights_kept[i] * (-func_logPhi_zero(μK[i], σ2K)) for i in zero_index) +
		 sum(weights_kept[i] * ((1 / 2) * log(σ2K) +
			(1 / 2) * (y_kept[i] - μK[i])^2 / σ2K) for i in positive_index)) /
		weight_total)

	JuMP.optimize!(model)
	termination = JuMP.termination_status(model)
	retried = false
	if !(termination in [JuMP.MOI.LOCALLY_SOLVED, JuMP.MOI.ALMOST_LOCALLY_SOLVED]) &&
		JuMP.has_values(model)
		retried = true
		retry_βK_b = value(βK_b)
		retry_βK_logC = value(βK_logC)
		retry_βK_a = value.(βK_a)
		retry_σ2K = value(σ2K)
		set_start_value(βK_b, retry_βK_b)
		set_start_value(βK_logC, retry_βK_logC)
		set_start_value.(βK_a, retry_βK_a)
		set_start_value(σ2K, retry_σ2K)
		set_optimizer_attribute(model, "acceptable_tol", 1e-4)
		set_optimizer_attribute(model, "acceptable_iter", 3)
		JuMP.optimize!(model)
		termination = JuMP.termination_status(model)
	end
	JuMP.has_values(model) || error("Prior-capacity Tobit failed: $termination")
	βK_a_est = fill(NaN, A_sb)
	βK_a_est[estimable_auctions] = value.(βK_a)
	return (
		  βK_est=vcat(value(βK_b), value(βK_logC), βK_a_est)
		, σ2K_est=value(σ2K)
		, active=active
		, estimable=estimable
		, positive_support=positive_support
		, raw_status=JuMP.raw_status(model)
		, termination_status=string(termination)
		, primal_status=string(JuMP.primal_status(model))
		, objective_value=JuMP.objective_value(model)
		, retried=retried
	)
end

# Weighted order-statistic likelihood for participant bids. Auction effects
# with no bootstrap support are returned as NaN and are not imputed.
function weighted_bid_distribution(sb_mat, N_sb, weights, βb_init, σ2b_init)
	b_all = Float64.(sb_mat[:, 2])
	br_all = Array{Int64}(undef, length(b_all))
	for a in a_sb_list
		br_all[a_sb .== a] = invperm(sortperm(b_all[a_sb .== a]))
	end

	with_N = .!ismissing.(N_sb)
	eligible_auctions = a_sb_list[.!ismissing.(N_sb_a)]
	active = [sum(weights[(a_sb .== a) .& with_N]) > 0 for a in eligible_auctions]
	active_auctions = eligible_auctions[active]
	keep = with_N .& (weights .> 0) .& in.(a_sb, Ref(Set(active_auctions)))
	any(keep) || error("No weighted winner observations remain for f(b*)")

	b = b_all[keep]
	br = br_all[keep]
	N = Int.(N_sb[keep])
	w = Float64.(weights[keep])
	weight_total = sum(w)
	position = Dict(a => k for (k, a) in enumerate(active_auctions))
	a_index = [position[a] for a in a_sb[keep]]

	model = Model(Ipopt.Optimizer)
	set_silent(model)
	set_optimizer_attribute(model, "max_iter", max_iter)
	set_optimizer_attribute(model, "tol", tol)
	func_Fb(b, μ, σ2) = min(max(cdf(Normal(μ, sqrt(σ2)), b), small_num), 1 - small_num)
	JuMP.register(model, :func_Fb, 3, func_Fb, autodiff=true)

	if ismissing(βb_init)
		@variable(model, βb[1:length(active_auctions)])
	else
		active_positions = indexin(active_auctions, eligible_auctions)
		@variable(model, βb[k=1:length(active_auctions)], start=βb_init[active_positions[k]])
	end
	@variable(model, σ2b_lb <= σ2b <= σ2b_ub, start=σ2b_init)
	@NLexpression(model, μb[i=eachindex(b)], βb[a_index[i]])
	@NLexpression(model, Fb[i=eachindex(b)], func_Fb(b[i], μb[i], σ2b))
	@NLexpression(model, nll[i=eachindex(b)],
		(1 / 2) * log(σ2b) + (1 / 2) * (b[i] - μb[i])^2 / σ2b -
		(br[i] - 1) * log(Fb[i]) - (N[i] - br[i]) * log(1 - Fb[i]))
	@NLobjective(model, Min, sum(w[i] * nll[i] for i in eachindex(b)) / weight_total)

	JuMP.optimize!(model)
	termination = JuMP.termination_status(model)
	JuMP.has_values(model) || error("Bid-price distribution optimization failed: $termination")
	βb_est = fill(NaN, length(eligible_auctions))
	βb_est[indexin(active_auctions, eligible_auctions)] = value.(βb)
	return (
		  βb_est=βb_est
		, σ2b_est=value(σ2b)
		, active=active
		, raw_status=JuMP.raw_status(model)
		, termination_status=string(termination)
		, primal_status=string(JuMP.primal_status(model))
		, objective_value=JuMP.objective_value(model)
	)
end

function func_dist(sb_mat, a_sb_list, N_sb, obs_sb, βb_init, σ2b_init;
	weights=ones(obs_sb), βK_init=missing, σ2K_start=σ2K_init)
	length(weights) == obs_sb || error("Distribution weights must have length $obs_sb")
	all(weights .>= 0) || error("Distribution weights must be nonnegative")
	sum(weights) > 0 || error("Distribution weights must have positive total weight")

	βC_est, σ2C_est, active_capacity =
		weighted_capacity_distribution(sb_mat, Float64.(weights))
	prior_result = weighted_prior_distribution(
		sb_mat, Float64.(weights); βK_init=βK_init, σ2K_start=σ2K_start)
	bid_result = weighted_bid_distribution(
		sb_mat, N_sb, Float64.(weights), βb_init, σ2b_init)

	return (
		  βC_est=βC_est
		, σ2C_est=σ2C_est
		, βK_est=prior_result.βK_est
		, σ2K_est=prior_result.σ2K_est
		, βb_est=bid_result.βb_est
		, σ2b_est=bid_result.σ2b_est
		, raw_status_K=prior_result.raw_status
		, termination_status_K=prior_result.termination_status
		, primal_status_K=prior_result.primal_status
		, objective_value_K=prior_result.objective_value
		, retried_K=prior_result.retried
		, raw_status_b=bid_result.raw_status
		, termination_status_b=bid_result.termination_status
		, primal_status_b=bid_result.primal_status
		, objective_value_b=bid_result.objective_value
		, capacity_active_auctions=count(active_capacity)
		, prior_active_auctions=count(prior_result.active)
		, prior_estimable_auctions=count(prior_result.estimable)
		, bid_active_auctions=count(bid_result.active)
		, weight_sum=sum(weights)
	)
end
