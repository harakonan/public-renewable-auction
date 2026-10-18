#=
Brazil new energy auction structural model estimation
first stage

@author: Konan Hara

@date: September 13, 2026

@software: Julia version 1.12
=#

# preparation code
include("prep.jl")
using NLopt

# Auction capacity defines the portfolio share and legal bounds. The transformed
# firm history K_log=log(1+K/100), with K measured in MW, is stored in column 3.
function build_step1_bidder_matrix(auction_data_local)
	return hcat(
		  Float64.(auction_data_local.bid_share)
		, Float64.(auction_data_local.bid_price)
		, log1p.(Float64.(auction_data_local.prior_started_wind_capacity) ./ 100.0)
		, Float64.(auction_data_local.bid_share_lb)
		, Float64.(auction_data_local.bid_share_ub)
	)
end

sb_mat = build_step1_bidder_matrix(auction_data_sb)
dc_mat = build_step1_bidder_matrix(auction_data_dc)

# Simple bootstrap control
# SLURM_ARRAY_TASK_ID maps directly to 'b'. If empty/missing, b = 0 (main estimation).
b = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "0"))
definitions_only = get(ENV, "STEP1_DEFINITIONS_ONLY", "false") == "true"

result_path(stem) = joinpath("output", "result", string(stem, ".csv"))
intermediate_path(stem) = joinpath(pathtodata, "intermediate", string(stem, ".csv"))

# Remove this draw's prior outputs before loading or calculating anything. A
# failed rerun must not leave stale files that appear to be current results.
if b > 0 && !definitions_only
	for output_dir in ("params1_bs", "params1_bs_status")
		rm(joinpath(
			pathtodata, "intermediate", output_dir, string("bs", b, ".csv")),
			force=true)
	end
end

# inverse cdf of descending-clock auction bid price truncated from above at purchase agreement price
inv_cdf_b_dc = deserialize(joinpath(pathtodata, "intermediate", "inv_cdf_b_dc.dat"))

# load specified parameters
params_specified = CSV.read(joinpath("output", "result", "params_specified.csv"), DataFrame)

# unpack parameters
μphat = filter(row -> occursin(r"mu_p_", row.variable), params_specified)[:,:value]
σ2phat = filter(row -> row.variable == "sigma2_p", params_specified)[:,:value][1]

# sealed-bid auction FOC (returns bid share)
# b: bid price
# (γ,μr): structural parameter
foc_sb = Array{Function}(undef, A_sb)
for a=1:A_sb
	function foc_sb_a(b, γ, μr)
		(1 - ((δ^l_sb_a[a]) * μr - δtilde_sb_a[a] * b) / (γ * σ2r_sb_a[a]))
	end
	foc_sb[a] = foc_sb_a
end

# descending-clock auction FOC (returns equality)
# q: bid share
# b: bid price
# (γ,μr): structural parameter
foc_dc = Array{Function}(undef, A_dc)
for a = 1:A_dc
	function foc_dc_a(q, b, γ, μr)
		z = (μphat[a] - q * γ * δtilde_dc_a[a] * σ2phat - b) / sqrt(σ2phat)
        μ_tilde = μphat[a] + sqrt(σ2phat) * (pdf(Normal(), z) / cdf(Normal(), z))
		(1 / (1 + (δtilde_dc_a[a]^2) * σ2phat / σ2r_dc_a[a])) *
		(1 - ((δ^l_dc_a[a]) * μr - δtilde_dc_a[a] * μ_tilde) / (γ * σ2r_dc_a[a])) - q
	end
	foc_dc[a] = foc_dc_a
end

# Gauss–Legendre nodes/weights on (0,1)
# Gauss–Legendre nodes/weights mapped from [-1,1] to (0,1). sum(w) ≈ 1.
# S: number of nodes/weights
function uniform_quadrature_legendre(S::Int)
	x, w = gausslegendre(S)         # nodes, weights on [-1,1]
	u = @. 0.5 * (x + 1.0)          # map to (0,1)
	w = @. 0.5 * w                  # rescale for (0,1) integral
	return u, w
end

# Compute lower-spike, upper-spike, and interior-mean moments.
# sb_mat, dc_mat: bidder variables
# b_dc_s: simulated descending clock bids
# IV: instrumental variables
# αγ, αr, αη: structural parameters
# w_b_dc_s: weights over draws
function func_g_obs(
	  sb_mat
	, dc_mat
	, b_dc_s
	, IV
	, αγ
	, αr
	, αη;
	  w_b_dc_s::AbstractVector
	, obs_weights=nothing
)
	# Bidder variables
	q_sb = sb_mat[:, 1]
	q_dc = dc_mat[:, 1]
	b_sb = sb_mat[:, 2]
	K_log_sb = sb_mat[:, 3]
	K_log_dc = dc_mat[:, 3]
	qlb_sb = sb_mat[:, 4]
	qlb_dc = dc_mat[:, 4]
	qub_sb = sb_mat[:, 5]
	qub_dc = dc_mat[:, 5]
	weights = isnothing(obs_weights) ? nothing : Float64.(obs_weights)
	isnothing(weights) || length(weights) == obs || throw(DimensionMismatch("obs_weights must match the number of observations"))

	# γ per observation
	if length(αγ) == 1
		# Homogeneous risk aversion
		γ_sb = fill(exp(αγ), obs_sb)
		γ_dc = fill(exp(αγ), obs_dc)
	else
		# Heterogeneous risk aversion
		γ_sb = [exp(αγ[1] + αγ[2] * K_log_sb[i]) for i = 1:obs_sb]
		γ_dc = [exp(αγ[1] + αγ[2] * K_log_dc[i]) for i = 1:obs_dc]
	end

	# μ_r per observation
	if length(αr) == 1
		μr_sb = fill(αr[1], obs_sb)
		μr_dc = fill(αr[1], obs_dc)
	elseif length(αr) == 2
		μr_sb = [αr[1] + αr[2] * m_sb[i] for i = 1:obs_sb]
		μr_dc = [αr[1] + αr[2] * m_dc[i] for i = 1:obs_dc]
	else
		throw(ArgumentError("αr must contain one or two parameters"))
	end

	# ση
	ση = exp(αη)
	stdN = Normal()

	# SB: lower/upper probabilities and conditional interior mean
	pL_sb = zeros(Float64, obs_sb)
	pU_sb = zeros(Float64, obs_sb)
	Eint_sb = zeros(Float64, obs_sb)

	for i in 1:obs_sb
		!isnothing(weights) && weights[i] == 0 && continue
		q_tilde = foc_sb[a_sb[i]](b_sb[i], γ_sb[i], μr_sb[i])
		α = (qlb_sb[i] - q_tilde) / ση
		β = (qub_sb[i] - q_tilde) / ση
		Φα = cdf(stdN, α)
		Φβ = cdf(stdN, β)
		ϕα = pdf(stdN, α)
		ϕβ = pdf(stdN, β)
		pL_sb[i] = Φα
		pU_sb[i] = 1.0 - Φβ
		pI = max(Φβ - Φα, small_num)
		μI = ση * (ϕα - ϕβ) / pI
		Eint_sb[i] = q_tilde + μI
	end

	# DC: lower/upper probabilities and conditional interior mean
	pL_dc = zeros(Float64, obs_dc)
	pU_dc = zeros(Float64, obs_dc)
	Eint_dc = zeros(Float64, obs_dc)

	# QUADRATURE mode (b_dc_s: obs_dc × S)
	S_local = size(b_dc_s, 2)
	w = w_b_dc_s

	for i = 1:obs_dc
		!isnothing(weights) && weights[obs_sb + i] == 0 && continue
		sum_pL = 0.0
		sum_pU = 0.0

		# Mixture over b for the interior:
		sum_pI = 0.0                 # Σ w_s pI_{i,s}
		sum_pI_mu = 0.0              # Σ w_s pI_{i,s} μ_{i,s}

		@inbounds for s = 1:S_local
			q_tilde = solve(
				  ZeroProblem(x -> foc_dc[a_dc[i]-A_sb](x, b_dc_s[i, s], γ_dc[i], μr_dc[i]), x0)
				, Order0()
			)

			α = (qlb_dc[i] - q_tilde) / ση
			β = (qub_dc[i] - q_tilde) / ση

			Φα = cdf(stdN, α)
			Φβ = cdf(stdN, β)
			ϕα = pdf(stdN, α)
			ϕβ = pdf(stdN, β)

			pI = max(Φβ - Φα, small_num)
			μI = ση * (ϕα - ϕβ) / pI
			μ_s = q_tilde + μI

			# Spikes: unconditional over b
			sum_pL += w[s] * Φα
			sum_pU += w[s] * (1.0 - Φβ)

			# Interior mixture: condition on Interior via pI-weights
			sum_pI += w[s] * pI
			sum_pI_mu += w[s] * pI * μ_s
		end

		pL_dc[i] = sum_pL
		pU_dc[i] = sum_pU

		Eint_dc[i] = sum_pI_mu / sum_pI
	end

	# Moments: lower bound, upper bound, and conditional interior mean
	K = size(IV, 2)
	gL = Array{Float64}(undef, obs, K)
	gU = Array{Float64}(undef, obs, K)
	gI = Array{Float64}(undef, obs, K)

	# SB indicators (tolerant comparisons)
	isL_sb = (q_sb .< (qlb_sb .+ small_num))
	isU_sb = (q_sb .> (qub_sb .- small_num))
	isI_sb = .!(isL_sb .| isU_sb)

	# SB blocks
	for i in 1:obs_sb, k in 1:K
		gL[i, k] = IV[i, k] * ((isL_sb[i] ? 1.0 : 0.0) - pL_sb[i])
		gU[i, k] = IV[i, k] * ((isU_sb[i] ? 1.0 : 0.0) - pU_sb[i])
		gI[i, k] = isI_sb[i] ? IV[i, k] * (q_sb[i] - Eint_sb[i]) : 0.0
	end

	# DC indicators
	isL_dc = (q_dc .< (qlb_dc .+ small_num))
	isU_dc = (q_dc .> (qub_dc .- small_num))
	isI_dc = .!(isL_dc .| isU_dc)

	# DC blocks
	for i in 1:obs_dc, k in 1:K
		row = obs_sb + i
		gL[row, k] = IV[row, k] * ((isL_dc[i] ? 1.0 : 0.0) - pL_dc[i])
		gU[row, k] = IV[row, k] * ((isU_dc[i] ? 1.0 : 0.0) - pU_dc[i])
		gI[row, k] = isI_dc[i] ? IV[row, k] * (q_dc[i] - Eint_dc[i]) : 0.0
	end

	return hcat(gL, gU, gI)
end

# Compute the sample mean of the observation-level moments.
function mean_g_obs(
	  IV
	, αγ
	, αr::AbstractVector{<:Real}
	, αη::Real;
	  b_dc_s
	, w_b_dc_s::AbstractVector
	, obs_weights=nothing
	, sb_mat_local=sb_mat
	, dc_mat_local=dc_mat
)
	g = func_g_obs(
		  sb_mat_local
		, dc_mat_local
		, b_dc_s
		, IV
		, αγ
		, αr
		, αη;
		  w_b_dc_s=w_b_dc_s
		, obs_weights=obs_weights
	)
	weights = isnothing(obs_weights) ? ones(Float64, size(g, 1)) : Float64.(obs_weights)
	length(weights) == size(g, 1) || throw(DimensionMismatch("obs_weights must match the number of observations"))
	weight_sum = sum(weights)
	weight_sum > 0 || throw(ArgumentError("obs_weights must have positive total weight"))
	return vec(sum(g .* reshape(weights, :, 1), dims=1) ./ weight_sum)
end

# Build GMM weight matrix at given parameter values
# IV: instrumental variables
# αγ, αr, αη: structural parameters used to evaluate moments
# b_dc_s: simulated descending clock bids
# w_b_dc_s: weights over draws
function build_W(
	  IV
	, αγ
	, αr::AbstractVector{<:Real}
	, αη::Real;
	  b_dc_s
	, w_b_dc_s::AbstractVector
	, obs_weights=nothing
	, sb_mat_local=sb_mat
	, dc_mat_local=dc_mat
)
	g = func_g_obs(
		  sb_mat_local
		, dc_mat_local
		, b_dc_s
		, IV
		, αγ
		, αr
		, αη;
		  w_b_dc_s=w_b_dc_s
		, obs_weights=obs_weights
	)
	weights = isnothing(obs_weights) ? ones(Float64, size(g, 1)) : Float64.(obs_weights)
	length(weights) == size(g, 1) || throw(DimensionMismatch("obs_weights must match the number of observations"))
	weight_sum = sum(weights)
	weight_sum > 0 || throw(ArgumentError("obs_weights must have positive total weight"))
	weighted_g = g .* reshape(weights, :, 1)
	gbar = vec(sum(weighted_g, dims=1) ./ weight_sum)
	g_centered = g .- reshape(gbar, 1, :)
	Σ = (g_centered' * (g_centered .* reshape(weights, :, 1))) / weight_sum
	return pinv(Symmetric(Σ))
end

# Estimate parameters by GMM using bounded Nelder-Mead.
# IV: instrumental variables
# αγ_init, αr_init, αη_init: initial parameters
# b_dc_s: simulated descending clock bids
# w_b_dc_s: weights over draws
# lower, upper: lower and upper bounds for parameters
# W: optional pre-calculated weight matrix
# sb_mat_local, dc_mat_local: bidder matrices used in moments
function estimate_gmm!(
      IV
    , αγ_init
    , αr_init::AbstractVector{<:Real}
    , αη_init::Real;
      b_dc_s
    , w_b_dc_s::AbstractVector
    , lower
    , upper
	, W = nothing
	, recenter_gbar = nothing
	, obs_weights=nothing
	, sb_mat_local=sb_mat
	, dc_mat_local=dc_mat
)

    # dimensions and arity
    αγ_len = isa(αγ_init, Number) ? 1 : length(αγ_init)
    αr_len = length(αr_init)
    arity = αγ_len + αr_len + 1

    # Determine Weight Matrix
    W_gmm = if W === nothing
		build_W(
			  IV
			, αγ_init
			, αr_init
			, αη_init;
			  b_dc_s=b_dc_s
			, w_b_dc_s=w_b_dc_s
			, obs_weights=obs_weights
			, sb_mat_local=sb_mat_local
			, dc_mat_local=dc_mat_local
		)
    else
        W
    end
    recenter_use = isnothing(recenter_gbar) ? zeros(Float64, size(W_gmm, 1)) : recenter_gbar
	weights_use = isnothing(obs_weights) ? ones(Float64, size(IV, 1)) : Float64.(obs_weights)
	length(weights_use) == size(IV, 1) || throw(DimensionMismatch("obs_weights must match the number of observations"))
	weight_sum = sum(weights_use)
	weight_sum > 0 || throw(ArgumentError("obs_weights must have positive total weight"))

    # inline objective function
    function gmm_obj_func(θ::AbstractVector)
        # inline unpacking
        αγ = (αγ_len == 1) ? θ[1] : θ[1:αγ_len]
        αr = θ[(αγ_len+1):(αγ_len+αr_len)]
        αη = θ[end]

        g = func_g_obs(
			  sb_mat_local
			, dc_mat_local
			, b_dc_s
			, IV
			, αγ
			, αr
			, αη;
			  w_b_dc_s=w_b_dc_s
			, obs_weights=obs_weights
		)
        ḡ = vec(sum(g .* reshape(weights_use, :, 1), dims=1) ./ weight_sum)
        ḡ_recentered = ḡ - recenter_use
        return dot(ḡ_recentered, W_gmm * ḡ_recentered)
    end

	# starting point θ0
	θ0 = vcat((αγ_len == 1 ? [αγ_init] : αγ_init), αr_init, αη_init)
	nlopt = NLopt.Opt(:LN_NELDERMEAD, arity)
	NLopt.lower_bounds!(nlopt, Float64.(lower))
	NLopt.upper_bounds!(nlopt, Float64.(upper))
	NLopt.min_objective!(nlopt, (θ, gradient) -> gmm_obj_func(θ))
	NLopt.initial_step!(nlopt, 0.01 .* (Float64.(upper) .- Float64.(lower)))
	NLopt.xtol_rel!(nlopt, 1e-7)
	NLopt.ftol_rel!(nlopt, 1e-10)
	NLopt.maxeval!(nlopt, 2000)
	f̂, θ̂, nlopt_status = NLopt.optimize(nlopt, θ0)
	solved = nlopt_status in (:SUCCESS, :STOPVAL_REACHED, :FTOL_REACHED, :XTOL_REACHED)
	solved || @warn "NLopt did not converge" nlopt_status
	return (
		  θ̂=θ̂
		, f̂=f̂
		, W=W_gmm
		, raw_status=String(nlopt_status)
		, termination_status=solved ? "LOCALLY_SOLVED" : "OTHER_ERROR"
		, optimizer="NLopt Nelder-Mead"
		, objective_evaluations=nlopt.numevals
		, maximum_iterations=missing
		, maximum_evaluations=2000
	)
end

# Run a two-step GMM with unchanged parameter dimension and flexible αγ block length len_ag
function two_step_gmm!(
      IV
    , αγ_init
    , αr_init::AbstractVector
    , αη_init;
      b_dc_s
    , w_b_dc_s::AbstractVector
    , lower
    , upper
    , len_ag::Int=(αγ_init isa AbstractVector ? length(αγ_init) : 1)
    , verbose::Bool=true
	, recenter_gbar = nothing
	, obs_weights=nothing
	, W_first = nothing
	, sb_mat_local=sb_mat
	, dc_mat_local=dc_mat
)
    # Helper to unpack step 1 results back into structural parameters
    function _split(θ)
        αγv = (len_ag == 1) ? θ[1] : θ[1:len_ag]
        αr_len = length(αr_init)
        αrv = θ[(len_ag+1):(len_ag+αr_len)]
        αηs = θ[end]
        return αγv, αrv, αηs
    end

    # Step 1
    res1 = estimate_gmm!(
          IV
        , αγ_init
        , αr_init
        , αη_init;
          b_dc_s=b_dc_s
        , w_b_dc_s=w_b_dc_s
        , lower=lower
        , upper=upper
		, W=W_first
		, recenter_gbar=recenter_gbar
		, obs_weights=obs_weights
		, sb_mat_local=sb_mat_local
		, dc_mat_local=dc_mat_local
	)

    if verbose
        println("Step 1: f̂ = ", res1.f̂, "  θ̂ = ", res1.θ̂)
    end

    # Step 2 starts from the Step-1 estimates and uses their optimal weight matrix.
    αγ1, αr1, αη1 = _split(res1.θ̂)
    W_step2 = build_W(
          IV
        , αγ1
        , αr1
        , αη1;
          b_dc_s=b_dc_s
		, w_b_dc_s=w_b_dc_s
		, obs_weights=obs_weights
		, sb_mat_local=sb_mat_local
        , dc_mat_local=dc_mat_local
    )

    res2 = estimate_gmm!(
          IV
        , αγ1
        , αr1
        , αη1;
          b_dc_s=b_dc_s
        , w_b_dc_s=w_b_dc_s
        , lower=lower
        , upper=upper
		, W=W_step2
		, recenter_gbar=recenter_gbar
		, obs_weights=obs_weights
		, sb_mat_local=sb_mat_local
		, dc_mat_local=dc_mat_local
	)

    if verbose
        println("Step 2: f̂ = ", res2.f̂, "  θ̂ = ", res2.θ̂)
    end

    return (
          θ̂=res2.θ̂
        , f̂=res2.f̂
        , W=res2.W
        , raw_status=res2.raw_status
        , termination_status=res2.termination_status
		, step1=res1
		, step2=res2
    )
end

# Calculates predicted bid shares for Goodness-of-Fit checks.
# θ̂: Vector of estimated parameters.
# αγ_len: 1 (homogeneous log gamma) or 2 (heterogeneous log gamma).
# b_dc_s: (obs_dc × S) matrix of price nodes.
# w_b_dc_s: Vector of length S with weights for b_dc_s.
function predict_bid_shares(
	  θ̂::AbstractVector{<:Real}
	, αγ_len::Int;
	  b_dc_s::AbstractMatrix
	, w_b_dc_s::AbstractVector
)
	# 0) Unpack θ̂
	αr_len = length(θ̂) - αγ_len - 1
	αγ = (αγ_len == 1) ? θ̂[1] : θ̂[1:αγ_len]
	αr = θ̂[(αγ_len+1):(αγ_len+αr_len)]
	ση = exp(θ̂[end])
	standard_normal = Normal()
	censored_normal_mean(q_tilde, lower, upper) = begin
		z_lower = (lower - q_tilde) / ση
		z_upper = (upper - q_tilde) / ση
		p_lower = cdf(standard_normal, z_lower)
		p_upper = ccdf(standard_normal, z_upper)
		p_interior = max(1.0 - p_lower - p_upper, 0.0)
		interior_first_moment =
			q_tilde * p_interior +
			ση * (pdf(standard_normal, z_lower) - pdf(standard_normal, z_upper))
		lower * p_lower + interior_first_moment + upper * p_upper
	end

	# 1) Pull data/vectors from memory
	q_sb = sb_mat[:, 1]
	b_sb = sb_mat[:, 2]
	K_log_sb = sb_mat[:, 3]
	qlb_sb_loc = sb_mat[:, 4]
	qub_sb_loc = sb_mat[:, 5]

	b_dc = dc_mat[:, 2]
	K_log_dc = dc_mat[:, 3]
	qlb_dc_loc = dc_mat[:, 4]
	qub_dc_loc = dc_mat[:, 5]

	# 2) γ per obs (αγ_len ∈ {1,2})
	if αγ_len == 1
		γ_sb = fill(exp(αγ), obs_sb)
		γ_dc = fill(exp(αγ), obs_dc)
	else
		γ_sb = [exp(αγ[1] + αγ[2] * K_log_sb[i]) for i = 1:obs_sb]
		γ_dc = [exp(αγ[1] + αγ[2] * K_log_dc[i]) for i = 1:obs_dc]
	end

	# 3) μ_r per obs
	if αr_len == 1
		μr_sb = fill(αr[1], obs_sb)
		μr_dc = fill(αr[1], obs_dc)
	elseif αr_len == 2
		μr_sb = [αr[1] + αr[2] * m_sb[i] for i = 1:obs_sb]
		μr_dc = [αr[1] + αr[2] * m_dc[i] for i = 1:obs_dc]
	else
		throw(ArgumentError("predict_bid_shares requires one or two αr parameters"))
	end

	# 4) Expected observed bid shares in SB auctions, integrating over the
	# censored-normal bid-share shock.
	q_latent_sb = similar(q_sb)
	@inbounds for i in 1:obs_sb
		q_latent_sb[i] = foc_sb[a_sb[i]](b_sb[i], γ_sb[i], μr_sb[i])
	end
	q_pred_sb = [
		censored_normal_mean(q_latent_sb[i], qlb_sb_loc[i], qub_sb_loc[i])
		for i in 1:obs_sb
	]

	# 5) Expected observed bid shares in DC auctions, integrating over both the
	# latent equilibrium bid price and the censored-normal bid-share shock.
	S_local = size(b_dc_s, 2)
	q_pred_dc = zeros(obs_dc)

	@inbounds for i in 1:obs_dc
		expected_q = 0.0
		for s in 1:S_local
			q_tilde = solve(ZeroProblem(x -> foc_dc[a_dc[i]-A_sb](x, b_dc_s[i, s], γ_dc[i], μr_dc[i]), x0), Order0())
			expected_q += w_b_dc_s[s] * censored_normal_mean(
				q_tilde, qlb_dc_loc[i], qub_dc_loc[i])
		end
		q_pred_dc[i] = expected_q
	end

	return (
		  q_pred_sb=q_pred_sb
		, q_pred_dc=q_pred_dc
		, b_sb=b_sb
		, b_dc=b_dc
		, K_log_sb=K_log_sb
		, K_log_dc=K_log_dc
		, a_sb=a_sb_id
		, a_dc=a_dc_id
	)
end


# descending-clock auction bid price integral by Gauss–Legendre quadrature
# winner bid prices and weights
U_b_dc_s, w_b_dc_s = uniform_quadrature_legendre(Sgmm)
b_dc_s = Array{Float64}(undef, obs_dc, Sgmm)

function result_by_step(result, step::Symbol)
	θ = step == :first ? result.step1.θ̂ : result.θ̂
	f = step == :first ? result.step1.f̂ : result.f̂
	return θ, f
end

function params1_values_homo_gamma_mur(result, step::Symbol)
	θ, f = result_by_step(result, step)
	return vcat(θ[1], missing, θ[2], missing, θ[3], f)
end

function params1_values_homo_gamma(result, step::Symbol)
	θ, f = result_by_step(result, step)
	return vcat(θ[1], missing, θ[2:3], θ[4], f)
end

function params1_values_hetero(result, step::Symbol)
	θ, f = result_by_step(result, step)
	return vcat(θ[1:5], f)
end

function params1_values_homo_gamma_mur(result)
	return vcat(result.θ̂[1], missing, result.θ̂[2], missing, result.θ̂[3], result.f̂)
end

function params1_values_homo_gamma(result)
	return vcat(result.θ̂[1], missing, result.θ̂[2:3], result.θ̂[4], result.f̂)
end

function params1_values_hetero(result)
	return vcat(result.θ̂[1:5], result.f̂)
end

function params1_status(result, model::Symbol)
	return DataFrame(
		  model=fill(String(model), 2)
		, step=["first", "second"]
		, raw_status=[result.step1.raw_status, result.step2.raw_status]
			, termination_status=[result.step1.termination_status, result.step2.termination_status]
			, optimizer=[result.step1.optimizer, result.step2.optimizer]
			, objective_evaluations=[result.step1.objective_evaluations, result.step2.objective_evaluations]
			, objective_value=[result.step1.f̂, result.step2.f̂]
			, moment_count=fill(size(result.W, 1), 2)
			, maximum_iterations=[result.step1.maximum_iterations, result.step2.maximum_iterations]
			, maximum_evaluations=[result.step1.maximum_evaluations, result.step2.maximum_evaluations]
		)
end

function params1_status_fixed_weight(result, model::Symbol)
	return DataFrame(
		  model=String(model)
		, step="fixed_weight"
		, raw_status=result.raw_status
		, termination_status=result.termination_status
		, optimizer=result.optimizer
		, objective_evaluations=result.objective_evaluations
		, objective_value=result.f̂
		, moment_count=size(result.W, 1)
		, maximum_iterations=result.maximum_iterations
		, maximum_evaluations=result.maximum_evaluations
	)
end

function fill_b_dc_s!(b_dc_s, inv_cdf_b_dc_local)
	for i = 1:obs_dc, s = 1:Sgmm
		b_dc_s[i, s] = inv_cdf_b_dc_local[i](U_b_dc_s[s])
	end
	return b_dc_s
end

function build_iv_homo_gamma_mur(sb_mat_local, dc_mat_local)
	return hcat(
		  vcat(fill(1.0, obs_sb), fill(0.0, obs_dc))
		, vcat(fill(0.0, obs_sb), fill(1.0, obs_dc))
		, vcat(sb_mat_local[:, 2], fill(0.0, obs_dc))
		, vcat(fill(0.0, obs_sb), dc_mat_local[:, 2])
	)
end

function build_iv_homo(sb_mat_local, dc_mat_local)
	return hcat(
		  build_iv_homo_gamma_mur(sb_mat_local, dc_mat_local)
		, m_vec
	)
end

function build_iv_hetero(sb_mat_local, dc_mat_local)
	price = vcat(sb_mat_local[:, 2], dc_mat_local[:, 2])
	K_log = vcat(sb_mat_local[:, 3], dc_mat_local[:, 3])
	return hcat(
		  build_iv_homo(sb_mat_local, dc_mat_local)
		, price .* K_log
	)
end

function identity_weight(IV)
	moment_count = 3 * size(IV, 2)
	return Matrix{Float64}(I, moment_count, moment_count)
end

function bounds_homo_gamma_mur()
	# The restricted bootstrap occasionally implies sigma_eta above 0.30.
	return [lower_homo[1], lower_homo[2], lower_homo[4]], [upper_homo[1], upper_homo[2], log(1.0)]
end

function bounds_homo_gamma()
	return lower_homo, upper_homo
end

function bounds_hetero()
	return lower_hetero, upper_hetero
end

if !definitions_only
if b == 0
	# ==============================================================================
	# Main Estimation (b = 0)
	# ==============================================================================
	println("--- Running Main Estimation ---")

    # winner bid prices for each auction and each simulation (for numerical integration in GMM)
	fill_b_dc_s!(b_dc_s, inv_cdf_b_dc)

	# IV
	IV_homo_gamma_mur = build_iv_homo_gamma_mur(sb_mat, dc_mat)
    IV_homo = build_iv_homo(sb_mat, dc_mat)
    IV_hetero = build_iv_hetero(sb_mat, dc_mat)

	# homogeneous risk aversion + homogeneous wholesale market price
	# starting values for the restricted model; subsequent models start from
	# their nested-model estimates
	αγ_init = log(0.5)
	αr_init = [25.0]
	αη_init = log(0.1)
	lower_homo_gamma_mur, upper_homo_gamma_mur = bounds_homo_gamma_mur()
	result_homo_gamma_mur = two_step_gmm!(IV_homo_gamma_mur, αγ_init, αr_init, αη_init;
        b_dc_s=b_dc_s
      , w_b_dc_s=w_b_dc_s
		, lower=lower_homo_gamma_mur
		, upper=upper_homo_gamma_mur
		, W_first=identity_weight(IV_homo_gamma_mur))

	# homogeneous risk aversion + heterogeneous wholesale market price
	αγ_init = result_homo_gamma_mur.θ̂[1]
	αr_init = [result_homo_gamma_mur.θ̂[2], 0.0]
	αη_init = result_homo_gamma_mur.θ̂[3]
	lower_homo_gamma, upper_homo_gamma = bounds_homo_gamma()
	W_homo_gamma = build_W(IV_homo, αγ_init, αr_init, αη_init;
		b_dc_s=b_dc_s, w_b_dc_s=w_b_dc_s)
	result_homo_gamma = two_step_gmm!(IV_homo, αγ_init, αr_init, αη_init;
		  b_dc_s=b_dc_s
		, w_b_dc_s=w_b_dc_s
		, lower=lower_homo_gamma
		, upper=upper_homo_gamma
		, W_first=W_homo_gamma)

	# prior-started-wind-capacity-dependent risk aversion
	αγ_init = [result_homo_gamma.θ̂[1], 0.0]
	αr_init = result_homo_gamma.θ̂[2:3]
	αη_init = result_homo_gamma.θ̂[4]
	lower_hetero_model, upper_hetero_model = bounds_hetero()
	W_hetero = build_W(IV_hetero, αγ_init, αr_init, αη_init;
		b_dc_s=b_dc_s, w_b_dc_s=w_b_dc_s)
    result_hetero = two_step_gmm!(IV_hetero, αγ_init, αr_init, αη_init;
        b_dc_s=b_dc_s
      , w_b_dc_s=w_b_dc_s
		, lower=lower_hetero_model
		, upper=upper_hetero_model
		, W_first=W_hetero)

	# Get GOF predictions for heterogeneous model
	preds = predict_bid_shares(result_hetero.θ̂, 2; b_dc_s=b_dc_s, w_b_dc_s=w_b_dc_s)

	# save for goodness-of-fit plot
	bid_share_predicted = (
		  auction_id=[preds.a_sb; preds.a_dc]
		, bid_price=[preds.b_sb; preds.b_dc]
		, prior_started_wind_capacity=[
			auction_data_sb.prior_started_wind_capacity
			; auction_data_dc.prior_started_wind_capacity
		]
		, bid_share=[sb_mat[:, 1]; dc_mat[:, 1]]
		, bid_share_predicted=[preds.q_pred_sb; preds.q_pred_dc]
	)
	CSV.write(intermediate_path("bid_share_predicted"), bid_share_predicted)

	# save results
	variable = vcat(
		# parameters
		  ["alpha_gamma0", "alpha_gammaK"]
		, [string("alpha_r", k) for k in 0:1]
		, ["alpha_eta"]
		# objective function value
		, ["objective_value"]
	)
	# estimates
	params1_hat = (
		  variable = variable
		, homo_gamma_mur_step1=params1_values_homo_gamma_mur(result_homo_gamma_mur, :first)
		, homo_gamma_mur=params1_values_homo_gamma_mur(result_homo_gamma_mur, :second)
		, homo_gamma_step1=params1_values_homo_gamma(result_homo_gamma, :first)
		, homo_gamma=params1_values_homo_gamma(result_homo_gamma, :second)
		, hetero_step1=params1_values_hetero(result_hetero, :first)
		, hetero=params1_values_hetero(result_hetero, :second)
	)
	CSV.write(result_path("params1_hat"), params1_hat)
	CSV.write(
		  result_path("params1_status")
		, vcat(
			  params1_status(result_homo_gamma_mur, :homo_gamma_mur)
			, params1_status(result_homo_gamma, :homo_gamma)
			, params1_status(result_hetero, :hetero)
		)
	)

else
	# ==============================================================================
	# Bootstrap Estimation (b > 0)
	# ==============================================================================
	println("--- Running Bootstrap Iteration $b ---")

	# load first stage parameter estimates
	params1_hat = CSV.read(result_path("params1_hat"), DataFrame)

	# unpack estimates for recentered moments
	# homogeneous risk aversion + homogeneous wholesale market price
	αγ_homo_gamma_mur_step1 = Float64.(filter(
		row -> row.variable == "alpha_gamma0", params1_hat
	)[:, :homo_gamma_mur_step1])[1]
	αr_homo_gamma_mur_step1 = Float64.(filter(
		row -> row.variable == "alpha_r0", params1_hat
	)[:, :homo_gamma_mur_step1])
	αη_homo_gamma_mur_step1 = Float64.(filter(
		row -> row.variable == "alpha_eta", params1_hat
	)[:, :homo_gamma_mur_step1])[1]
	αγ_homo_gamma_mur = Float64.(filter(row -> row.variable == "alpha_gamma0", params1_hat)[:, :homo_gamma_mur])[1]
	αr_homo_gamma_mur = Float64.(filter(row -> row.variable == "alpha_r0", params1_hat)[:, :homo_gamma_mur])
	αη_homo_gamma_mur = Float64.(filter(row -> row.variable == "alpha_eta", params1_hat)[:, :homo_gamma_mur])[1]
	# homogeneous risk aversion + heterogeneous wholesale market price
	αγ_homo_gamma_step1 = Float64.(filter(
		row -> row.variable == "alpha_gamma0", params1_hat
	)[:, :homo_gamma_step1])[1]
	αr_homo_gamma_step1 = Float64.(filter(
		row -> row.variable in ("alpha_r0", "alpha_r1"), params1_hat
	)[:, :homo_gamma_step1])
	αη_homo_gamma_step1 = Float64.(filter(
		row -> row.variable == "alpha_eta", params1_hat
	)[:, :homo_gamma_step1])[1]
	αγ_homo_gamma = Float64.(filter(row -> row.variable == "alpha_gamma0", params1_hat)[:, :homo_gamma])[1]
	αr_homo_gamma = Float64.(filter(
		row -> row.variable in ("alpha_r0", "alpha_r1"), params1_hat
	)[:, :homo_gamma])
	αη_homo_gamma = Float64.(filter(row -> row.variable == "alpha_eta", params1_hat)[:, :homo_gamma])[1]
	# prior-started-wind-capacity-dependent risk aversion
	αγ_hetero_step1 = Float64.(filter(
		row -> occursin(r"alpha_gamma", row.variable), params1_hat
	)[:, :hetero_step1])
	αr_hetero_step1 = Float64.(filter(
		row -> row.variable in ("alpha_r0", "alpha_r1"), params1_hat
	)[:, :hetero_step1])
	αη_hetero_step1 = Float64.(filter(
		row -> row.variable == "alpha_eta", params1_hat
	)[:, :hetero_step1])[1]
	αγ_hetero = Float64.(filter(row -> occursin(r"alpha_gamma", row.variable), params1_hat)[:, :hetero])
	αr_hetero = Float64.(filter(
		row -> row.variable in ("alpha_r0", "alpha_r1"), params1_hat
	)[:, :hetero])
	αη_hetero = Float64.(filter(row -> row.variable == "alpha_eta", params1_hat)[:, :hetero])[1]

	# recenter bootstrap moments by the baseline-sample moments at the original estimates
	b_dc_s_full = fill_b_dc_s!(Array{Float64}(undef, obs_dc, Sgmm), inv_cdf_b_dc)
	IV_homo_gamma_mur = build_iv_homo_gamma_mur(sb_mat, dc_mat)
	IV_homo = build_iv_homo(sb_mat, dc_mat)
	IV_hetero = build_iv_hetero(sb_mat, dc_mat)
	recenter_gbar_homo_gamma_mur = mean_g_obs(
		  IV_homo_gamma_mur
		, αγ_homo_gamma_mur
		, αr_homo_gamma_mur
		, αη_homo_gamma_mur;
		  b_dc_s=b_dc_s_full
		, w_b_dc_s=w_b_dc_s
	)
	recenter_gbar_homo_gamma = mean_g_obs(
		  IV_homo
		, αγ_homo_gamma
		, αr_homo_gamma
		, αη_homo_gamma;
		  b_dc_s=b_dc_s_full
		, w_b_dc_s=w_b_dc_s
	)
	recenter_gbar_hetero = mean_g_obs(
		  IV_hetero
		, αγ_hetero
		, αr_hetero
		, αη_hetero;
		  b_dc_s=b_dc_s_full
		, w_b_dc_s=w_b_dc_s
	)

	# Fix each bootstrap criterion at the original-sample second-step weighting
	# matrix, evaluated at that model's original first-step estimate.
	W_homo_gamma_mur_fixed = build_W(
		  IV_homo_gamma_mur
		, αγ_homo_gamma_mur_step1
		, αr_homo_gamma_mur_step1
		, αη_homo_gamma_mur_step1;
		  b_dc_s=b_dc_s_full
		, w_b_dc_s=w_b_dc_s
	)
	W_homo_gamma_fixed = build_W(
		  IV_homo
		, αγ_homo_gamma_step1
		, αr_homo_gamma_step1
		, αη_homo_gamma_step1;
		  b_dc_s=b_dc_s_full
		, w_b_dc_s=w_b_dc_s
	)
	W_hetero_fixed = build_W(
		  IV_hetero
		, αγ_hetero_step1
		, αr_hetero_step1
		, αη_hetero_step1;
		  b_dc_s=b_dc_s_full
		, w_b_dc_s=w_b_dc_s
	)

	# output directory
	bootstrap_dir = "params1_bs"
	bootstrap_status_dir = "params1_bs_status"
	mkpath(joinpath(pathtodata, "intermediate", bootstrap_dir))
	mkpath(joinpath(pathtodata, "intermediate", bootstrap_status_dir))

	# bootstrap data
	seed_bs = seed_bs_list[b]
	Random.seed!(seed_bs)
	obs_weights_bs, bootstrap_diagnostics = draw_global_firm_weights()
	println("Global firm-pairs bootstrap diagnostics: ", bootstrap_diagnostics)

	# Firm resampling changes observation multiplicities, not covariates or
	# auction-specific numerical-integration objects.
	fill_b_dc_s!(b_dc_s, inv_cdf_b_dc)
	IV_homo_gamma_mur_bs = build_iv_homo_gamma_mur(sb_mat, dc_mat)
    IV_homo_bs = build_iv_homo(sb_mat, dc_mat)
	IV_hetero_bs = build_iv_hetero(sb_mat, dc_mat)

	# Fixed-weight bootstrap estimates start from the corresponding
	# original-sample second-step estimates.
	lower_homo_gamma_mur, upper_homo_gamma_mur = bounds_homo_gamma_mur()
	result_homo_gamma_mur = estimate_gmm!(
		  IV_homo_gamma_mur_bs
		, αγ_homo_gamma_mur
		, αr_homo_gamma_mur
		, αη_homo_gamma_mur;
          b_dc_s=b_dc_s
        , w_b_dc_s=w_b_dc_s
		, lower=lower_homo_gamma_mur
		, upper=upper_homo_gamma_mur
		, W=W_homo_gamma_mur_fixed
		, obs_weights=obs_weights_bs
		, recenter_gbar=recenter_gbar_homo_gamma_mur)

	# homogeneous risk aversion + heterogeneous wholesale market price
	lower_homo_gamma, upper_homo_gamma = bounds_homo_gamma()
	result_homo_gamma = estimate_gmm!(
		  IV_homo_bs
		, αγ_homo_gamma
		, αr_homo_gamma
		, αη_homo_gamma;
		  b_dc_s=b_dc_s
		, w_b_dc_s=w_b_dc_s
		, lower=lower_homo_gamma
		, upper=upper_homo_gamma
		, W=W_homo_gamma_fixed
		, obs_weights=obs_weights_bs
		, recenter_gbar=recenter_gbar_homo_gamma
	)

	# prior-started-wind-capacity-dependent risk aversion
	lower_hetero_model, upper_hetero_model = bounds_hetero()
	result_hetero = estimate_gmm!(
		  IV_hetero_bs
		, αγ_hetero
		, αr_hetero
		, αη_hetero;
          b_dc_s=b_dc_s
        , w_b_dc_s=w_b_dc_s
		, lower=lower_hetero_model
		, upper=upper_hetero_model
		, W=W_hetero_fixed
		, obs_weights=obs_weights_bs
        , recenter_gbar=recenter_gbar_hetero)

	# save results
	output = (
		  variable=params1_hat.variable
		, homo_gamma_mur_step1=fill(missing, nrow(params1_hat))
		, homo_gamma_mur=params1_values_homo_gamma_mur(result_homo_gamma_mur)
		, homo_gamma_step1=fill(missing, nrow(params1_hat))
		, homo_gamma=params1_values_homo_gamma(result_homo_gamma)
		, hetero_step1=fill(missing, nrow(params1_hat))
		, hetero=params1_values_hetero(result_hetero)
	)
	CSV.write(
		joinpath(pathtodata, "intermediate", bootstrap_dir, string("bs", b, ".csv"))
		, output
	)
	status_output = vcat(
		  params1_status_fixed_weight(result_homo_gamma_mur, :homo_gamma_mur)
		, params1_status_fixed_weight(result_homo_gamma, :homo_gamma)
		, params1_status_fixed_weight(result_hetero, :hetero)
	)
	status_output[!, :bootstrap_draw] = fill(b, nrow(status_output))
	for (name, value) in pairs(bootstrap_diagnostics)
		status_output[!, name] = fill(value, nrow(status_output))
	end
	CSV.write(
		  joinpath(
			pathtodata
			, "intermediate"
			, bootstrap_status_dir
			, string("bs", b, ".csv")
		  )
		, status_output
	)
end
end
