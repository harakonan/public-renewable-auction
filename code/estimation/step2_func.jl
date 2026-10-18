#=
Brazil new energy auction structural model estimation
sealed-bid auction winning probability function and cost distribution
functions

@author: Konan Hara

@date: September 13, 2026

@software: Julia version 1.12

=#

# sealed-bid auction unconstrained optimal bid share function
# b: bid price
# (γ,μr): structural parameter
func_q_tilde_sb = Array{Function}(undef, A_sb)
for a = 1:A_sb
    function func_q_tilde_sb_a(b, γi, αr)
        q_tilde_sb = (1 - ((δ^l_sb_a[a]) * (αr[1] + αr[2] * m_sb_a[a]) - δtilde_sb_a[a] * b) / (γi * σ2r_sb_a[a]))
        return q_tilde_sb
    end
    func_q_tilde_sb[a] = func_q_tilde_sb_a
end

# function that generates auction parameters
# a: auction ID
# (βC,σ2C,βK,σ2K,αγ,αr,βb,σ2b): structural parameters
function func_aparam(a, βC, σ2C, βK, σ2K, αγ, αr, βb, σ2b)
    func_q_tilde(b, γi) = func_q_tilde_sb[a](b, γi, αr)
    aparam = (
        # f(log D)
        μD=βDhat[1] + βDhat[2] * (t_sb_a[a] + l_sb_a[a]) + βDhat[3] * m_sb_a[a] + βDhat[4] * N_sb_a[a]
        , σ2D=σ2Dhat
        # f(r)
        , μr=αr[1] + αr[2] * m_sb_a[a]
        , σ2r=σ2r_sb_a[a]
		# f(log C | b*,a)
		, βC_b=βC[1]
		, μC_b0=βC[1 + a]
		, σ2C=σ2C
		# censored-normal f(K_log | C,b*,a), K_log=log(1+K/100)
		, βK_b=βK[1]
		, βK_logC=βK[2]
		, μK_bC0=βK[2 + a]
		, σ2K=σ2K
		, K_degenerate=structural_zero_prior_sb_a[a] || !isfinite(βK[2 + a])
		# log gamma = alpha_gamma0 + alpha_gammaK * K_log
		, αγ_K=αγ[2]
		, μγ_prior0=αγ[1]
        # f(q|b,γ)
        , func_q_tilde=func_q_tilde
        # f(b)
        , μb=βb[a]
        , σ2b=σ2b
        # number of participants
        , N=N_sb_a[a]
        # bid share lower bound
        , qlb_a=qlb_sb_a[a]
        # quantity minimum unit
        , Qmu=Qmu_sb_a[a]
    )
    return aparam
end

function draw_K_log(b, C, aparam, ZK)
	aparam.K_degenerate && return zeros(size(b))
	μK = aparam.μK_bC0 .+ aparam.βK_b .* b .+
		aparam.βK_logC .* log.(C)
	y_latent = μK .+ sqrt(aparam.σ2K) .* ZK
	return max.(0.0, y_latent)
end

function markup_from_winprob(b, q, γ, a, W, dW)
	slope = dW[a](b)
	winprob = W[a](b)
	if !isfinite(slope) || !isfinite(winprob) || slope >= 0
		return missing
	end
	log_argument = 1 - (winprob * γ * q * δtilde_sb_a[a]) / slope
	if !isfinite(log_argument) || log_argument <= 0
		return missing
	end
	return log(log_argument) / γ
end

function safe_cor(x, y)
	length(x) > 1 || return missing
	sx = std(x)
	sy = std(y)
	tol_x = sqrt(eps(Float64)) * max(abs(mean(x)), 1.0)
	tol_y = sqrt(eps(Float64)) * max(abs(mean(y)), 1.0)
	(sx > tol_x && sy > tol_y) || return missing
	return cor(x, y)
end

const projection_terms = ["intercept", "t", "t2", "t3", "l", "t_l", "t2_l", "t3_l"]
const joint_projection_parameters = vcat(
	["loggamma_" * term for term in projection_terms],
	["cost_" * term for term in projection_terms],
	["var_loggamma", "cov_loggamma_cost", "var_cost"],
)
const joint_summary_parameters = [
	"loggamma_mean", "loggamma_var", "gamma_mean", "gamma_var",
	"cost_mean", "cost_var", "loggamma_cost_cor", "gamma_cost_cor",
]
const distribution_weightings = ["equal_auction", "participant", "capacity"]

function joint_projection_inputs(γ_sim, c_sim, Capacity_sim, a_list)
	X_parts = Matrix{Float64}[]
	logγ_parts = Vector{Float64}[]
	γ_parts = Vector{Float64}[]
	c_parts = Vector{Float64}[]
	weight_parts = Dict(name => Vector{Float64}[] for name in distribution_weightings)

	for a in a_list
		valid = .!ismissing.(γ_sim[:, a]) .& .!ismissing.(c_sim[:, a]) .&
			.!ismissing.(Capacity_sim[:, a])
		any(valid) || continue
		γ_a = Float64.(γ_sim[valid, a])
		c_a = Float64.(c_sim[valid, a])
		Capacity_a = Float64.(Capacity_sim[valid, a])
		n_a = length(γ_a)
		t = Float64(t_sb_a[a])
		l = Float64(l_sb_a[a])
		x_a = [1.0, t, t^2, t^3, l, t * l, t^2 * l, t^3 * l]

		push!(X_parts, repeat(reshape(x_a, 1, :), n_a, 1))
		push!(logγ_parts, log.(γ_a))
		push!(γ_parts, γ_a)
		push!(c_parts, c_a)
		push!(weight_parts["equal_auction"], fill(1.0 / n_a, n_a))
		push!(weight_parts["participant"], fill(Float64(N_sb_a[a]) / n_a, n_a))
		total_capacity = Float64(auction_sb_a[a, :total_capacity])
		push!(weight_parts["capacity"], total_capacity .* Capacity_a ./ sum(Capacity_a))
	end

	isempty(X_parts) && error("No simulated participant costs were recovered")
	return (
		  X=reduce(vcat, X_parts)
		, logγ=reduce(vcat, logγ_parts)
		, γ=reduce(vcat, γ_parts)
		, c=reduce(vcat, c_parts)
		, weights=Dict(name => reduce(vcat, weight_parts[name]) for name in distribution_weightings)
	)
end

function weighted_joint_projection(data, weights)
	X = data.X
	length(weights) == size(X, 1) || error("Projection weights do not align")
	XtWX = X' * (weights .* X)
	rank(XtWX) == size(X, 2) || return fill(NaN, length(joint_projection_parameters))
	βlogγ = XtWX \ (X' * (weights .* data.logγ))
	βc = XtWX \ (X' * (weights .* data.c))
	e_logγ = data.logγ - X * βlogγ
	e_c = data.c - X * βc
	w = weights ./ sum(weights)
	σ2logγ = sum(w .* e_logγ .^ 2)
	σlogγc = sum(w .* e_logγ .* e_c)
	σ2c = sum(w .* e_c .^ 2)
	return [βlogγ; βc; σ2logγ; σlogγc; σ2c]
end

function weighted_joint_summary(data, weights)
	w = weights ./ sum(weights)
	μlogγ = sum(w .* data.logγ)
	μγ = sum(w .* data.γ)
	μc = sum(w .* data.c)
	logγ_dev = data.logγ .- μlogγ
	γ_dev = data.γ .- μγ
	c_dev = data.c .- μc
	var_logγ = sum(w .* logγ_dev .^ 2)
	var_γ = sum(w .* γ_dev .^ 2)
	var_c = sum(w .* c_dev .^ 2)
	cor_logγ_c = sum(w .* logγ_dev .* c_dev) / sqrt(var_logγ * var_c)
	cor_γ_c = sum(w .* γ_dev .* c_dev) / sqrt(var_γ * var_c)
	return [μlogγ, var_logγ, μγ, var_γ, μc, var_c, cor_logγ_c, cor_γ_c]
end

function joint_distribution_outputs(γ_sim, c_sim, Capacity_sim, a_list)
	data = joint_projection_inputs(γ_sim, c_sim, Capacity_sim, a_list)
	projection = hcat([
		weighted_joint_projection(data, data.weights[name])
		for name in distribution_weightings
	]...)
	summary = hcat([
		weighted_joint_summary(data, data.weights[name])
		for name in distribution_weightings
	]...)
	return projection, summary
end

# smoothing function Φ(x/h)
func_I(x,h) = cdf(Normal(), x / h)

# draw auction-level inputs used throughout the winning-probability grid
# aparam: auction parameters
# Umax: standard normal draws for competitors
# UD: standard normal draws for procurement quantity
function draw_W_inputs(aparam, Umax, UD, UK)
    # unpack parameters
    μD = aparam.μD
    σ2D = aparam.σ2D
	βC_b = aparam.βC_b
	μC_b0 = aparam.μC_b0
	σ2C = aparam.σ2C
	αγ_K = aparam.αγ_K
    μγ_prior0 = aparam.μγ_prior0
    func_q_tilde = aparam.func_q_tilde
    μb = aparam.μb
    σ2b = aparam.σ2b
    N = aparam.N
    qlb_a = aparam.qlb_a
    Qmu = aparam.Qmu

    # draw procurement quantity
    D_draws = exp.(μD .+ UD * sqrt(σ2D))

    # draw competitors
    b_draws = Array{Float64}(undef, SW, N - 1)
    Q_draws = Array{Float64}(undef, SW, N - 1)
	for sw = 1:SW, j = 1:N-1
		b_draws[sw, j] = μb + Umax[1, sw, j] * sqrt(σ2b)
	end
	Capacity_draws = exp.(βC_b .* b_draws .+ μC_b0 .+
		sqrt(σ2C) .* @view(Umax[2, :, 1:N-1]))
	K_log_draws = draw_K_log(
		b_draws, Capacity_draws, aparam, @view(UK[:, 1:N-1]))
    for sw = 1:SW
        for j = 1:N-1
			# bidder variables
            b_sj = b_draws[sw, j]
            Capacity_sj = Capacity_draws[sw, j]
			# γ
			γ_sj = exp(μγ_prior0 + K_log_draws[sw, j] * αγ_K)
			# unconstrained optimal bid share
            q_tilde_sj = func_q_tilde(b_sj, γ_sj)
            # optimal bid share adjusted for quantity minimum unit
            qlb_sj = (ceil(Capacity_sj * qlb_a / Qmu - eps(1.0)) * Qmu) / Capacity_sj
            qub_sj = (floor(Capacity_sj / Qmu + eps(1.0)) * Qmu) / Capacity_sj
            q_ast_sj = min(max(q_tilde_sj, qlb_sj), qub_sj)
            Q_sj = q_ast_sj * Capacity_sj
            b_draws[sw, j] = b_sj
            Q_draws[sw, j] = Q_sj
        end
    end

	return D_draws, b_draws, Q_draws
end

# function that calculates a winning probability nonparametrically
function func_W_nonparametric(b, D_draws, b_draws, Q_draws)
    mean([func_I(D_draws[sw] - sum(Q_draws[sw, :] .* func_I.(b .- b_draws[sw, :], hb)), hD) for sw = 1:SW])
end

# function that generates winning probability functions with linear interpolation
# aparam: auction parameters
# Umax: standard normal draws
function func_W_functions(aparam, Umax, UD, UK)
	D_draws, b_draws, Q_draws = draw_W_inputs(aparam, Umax, UD, UK)
	W_at(b) = func_W_nonparametric(b, D_draws, b_draws, Q_draws)

    # refine bid price lower/upper bound by the winning probability upper/lower bound
	index = 0
	distance_W_ub = 0.0
	b_lb = first(b_range)
	while distance_W_ub >= 0 && index < length(b_range)
		index += 1
		distance_W_ub = W_at(b_range[index]) - W_ub
		b_lb = index == 1 ? first(b_range) : b_range[index - 1]
	end
	distance_W_ub < 0 || error(
		"Winning probability remains above W_ub at the upper bid-price bound")

	distance_W_lb = 0.0
	b_ub = last(b_range)
	while distance_W_lb >= 0 && index < length(b_range)
		index += 1
		distance_W_lb = W_at(b_range[index]) - W_lb
		b_ub = b_range[index]
	end
	use_asymptote_fallback = distance_W_lb >= 0
	W_inf = missing
	if use_asymptote_fallback
		W_inf = mean([
			func_I(D_draws[sw] - sum(Q_draws[sw, :]), hD)
			for sw in 1:SW
		])
		distance_W_inf = W_at(last(b_range)) - W_inf
		distance_W_inf < W_asymptote_tol || error(
			"Winning probability neither falls below W_lb nor converges " *
			"to its lower asymptote at the upper bid-price bound")
	end

    # calculate W
    W_b_grid = b_lb:grid_interval:b_ub
    W_grid = [W_at(b) for b = W_b_grid]

    return (
		W_b_grid=W_b_grid,
		W_grid=W_grid,
		W_inf=W_inf,
		use_asymptote_fallback=use_asymptote_fallback,
    )
end

# function for step2 estimation
# sb_mat: sealed-bid bidder matrix
# a_sb_local: auction ID for each row of sb_mat
# a_sb_list_local: sealed-bid auctions to use in step2
# (βC,σ2C,βK,σ2K,γ,αr,βb,σ2b): structural parameters
# winprob_grid: winning probability functions for sealed-bid auctions
# Uγc: draws for bid and capacity in the cost-distribution simulation
# prior_assurance_100_factual: observed prior-started wind capacity / 100 MW
# UK_dist: standard-normal draws for the shifted-lognormal Tobit
function func_step2(sb_mat, a_sb_local, a_sb_list_local, βC, σ2C, βK, σ2K,
	γ, αr, βb, σ2b, winprob_grid, Uγc, prior_assurance_100_factual, UK_dist)

    @assert size(sb_mat, 1) == length(a_sb_local)
    @assert size(sb_mat, 1) == length(prior_assurance_100_factual)

    # generate W and dW/db function from the grid information
    W = Array{Any}(undef, A_sb)
    dW = Array{Any}(undef, A_sb)
    for a in a_sb_list_local
        b_grid = winprob_grid[winprob_grid.auction_id.==a, :bid_price]
        W_grid = winprob_grid[winprob_grid.auction_id.==a, :winprob]
		fallback_a = unique(winprob_grid[
			winprob_grid.auction_id.==a, :use_asymptote_fallback])
		length(fallback_a) == 1 || error(
			"Winning-probability fallback indicator is not unique within auction $a")

        # truncate winning probability function
		if fallback_a[1]
			W_inf_a = unique(skipmissing(winprob_grid[
				winprob_grid.auction_id.==a, :winprob_inf]))
			length(W_inf_a) == 1 || error(
				"Winning-probability asymptote is not unique within auction $a")
			valid_lower_grid = W_grid .- W_inf_a[1] .> W_asymptote_tol
		else
			valid_lower_grid = W_grid .> W_lb
		end
		valid_grid = (W_grid .< W_ub) .& valid_lower_grid
		b_grid = b_grid[valid_grid]
		W_grid = W_grid[valid_grid]

        # calculate W'(=dW/db) with approximate derivative
        dW_grid = [(W_grid[k+1] - W_grid[k]) / grid_interval for k = 1:length(b_grid)-1]
        # grid points for the derivative
        db_grid = b_grid[1:end-1] .+ grid_interval / 2

        # approximate W and dW/db using linear interpolation
        # (extrapolate outside flat)
        W[a] = linear_interpolation(b_grid, W_grid, extrapolation_bc=Interpolations.Flat())
        dW[a] = linear_interpolation(db_grid, dW_grid, extrapolation_bc=Interpolations.Flat())
    end

    # bidder variables
    b_sb = sb_mat[:, 2]
    Capacity_sb = sb_mat[:, 3]
    qlb_sb = sb_mat[:, 4]
    qub_sb = sb_mat[:, 5]
    obs_sb_local = size(sb_mat, 1)

    # recover (γ,c) for factual observations
    γ_f = Array{Union{Missing,Float64}}(missing, obs_sb_local)
    q_ast_f = Array{Union{Missing,Float64}}(missing, obs_sb_local)
    rp_f = Array{Union{Missing,Float64}}(missing, obs_sb_local)
    markup_f = Array{Union{Missing,Float64}}(missing, obs_sb_local)
    c_f = Array{Union{Missing,Float64}}(missing, obs_sb_local)
    # simulate the entire (γ,c) distribution
    γ_sim = Array{Union{Missing,Float64}}(missing, (Sγc, A_sb))
	Capacity_sim = Array{Union{Missing,Float64}}(missing, (Sγc, A_sb))
    c_sim = Array{Union{Missing,Float64}}(missing, (Sγc, A_sb))
    # undiscounted for the lead time
    c_sim_undiscount = Array{Union{Missing,Float64}}(missing, (Sγc, A_sb))
    for a in a_sb_list_local
        # auction parameters
        aparam = func_aparam(a, βC, σ2C, βK, σ2K, γ, αr, βb, σ2b)

        # unpack parameters
        μr = aparam.μr
        σ2r = aparam.σ2r
		βC_b = aparam.βC_b
		μC_b0 = aparam.μC_b0
		σ2C_a = aparam.σ2C
		αγ_K = aparam.αγ_K
        μγ_prior0 = aparam.μγ_prior0
        func_q_tilde = aparam.func_q_tilde
        μb = aparam.μb
        σ2b = aparam.σ2b
        qlb_a = aparam.qlb_a
        Qmu = aparam.Qmu

        # factual
        a_mask = a_sb_local .== a
        NW = sum(a_mask)
        γ_a = Array{Union{Missing,Float64}}(missing, NW)
        q_ast_a = Array{Union{Missing,Float64}}(missing, NW)
        rp_a = Array{Union{Missing,Float64}}(missing, NW)
        markup_a = Array{Union{Missing,Float64}}(missing, NW)
        c_a = Array{Union{Missing,Float64}}(missing, NW)
        for i in 1:NW
            # bidder variables
            b_i = b_sb[a_mask][i]
            Capacity_i = Capacity_sb[a_mask][i]
			# γ
			K_log_i = log1p(prior_assurance_100_factual[a_mask][i])
			γ_i = exp(μγ_prior0 + K_log_i * αγ_K)
			# unconstrained optimal bid share
            q_tilde_i = func_q_tilde(b_i, γ_i)
			# adjust for quantity minimum unit
            q_ast_i = min(max(q_tilde_i, qlb_sb[a_mask][i]), qub_sb[a_mask][i])
			# risk premium
            rp_i = (γ_i / 2) * σ2r
			# markup
			markup_i = markup_from_winprob(b_i, q_ast_i, γ_i, a, W, dW)
			γ_a[i] = γ_i
			q_ast_a[i] = q_ast_i
			rp_a[i] = rp_i
			ismissing(markup_i) && continue
			# cost
            c_i = q_ast_i * δtilde_sb_a[a] * b_i + (1 - q_ast_i) * μr - (γ_i / 2) * ((1 - q_ast_i)^2) * σ2r - markup_i
			# store results
            markup_a[i] = markup_i
            c_a[i] = c_i
        end
		# store results
        γ_f[a_mask] = γ_a
        q_ast_f[a_mask] = q_ast_a
        rp_f[a_mask] = rp_a
        markup_f[a_mask] = markup_a
        c_f[a_mask] = c_a

		# Simulate f(gamma,C,b*) = f(gamma,C | b*) f(b*)
		b_sim_a = μb .+ Uγc[1, :] .* sqrt(σ2b)
		Capacity_sim_a = exp.(βC_b .* b_sim_a .+ μC_b0 .+
			Uγc[2, :] .* sqrt(σ2C_a))
		K_log_sim_a = draw_K_log(
			b_sim_a, Capacity_sim_a, aparam, UK_dist)
        for s in 1:Sγc
			# bidder variables
            b_s = b_sim_a[s]
            Capacity_s = Capacity_sim_a[s]
            # γ
			γ_s = exp(μγ_prior0 + K_log_sim_a[s] * αγ_K)
            # unconstrained optimal bid share
            q_tilde_s = func_q_tilde(b_s, γ_s)
            # optimal bid share adjusted for quantity minimum unit
            qlb_s = (ceil(Capacity_s * qlb_a / Qmu - eps(1.0)) * Qmu) / Capacity_s
            qub_s = (floor(Capacity_s / Qmu + eps(1.0)) * Qmu) / Capacity_s
            q_ast_s = min(max(q_tilde_s, qlb_s), qub_s)
            # markup
			markup_s = markup_from_winprob(b_s, q_ast_s, γ_s, a, W, dW)
			γ_sim[s, a] = γ_s
			ismissing(markup_s) && continue
            # cost
            c_s = q_ast_s * δtilde_sb_a[a] * b_s + (1 - q_ast_s) * μr - (γ_s / 2) * ((1 - q_ast_s)^2) * σ2r - markup_s
            # undiscounted cost
            c_s_undiscount = (1 / (δ^l_sb_a[a])) * (q_ast_s * δtilde_sb_a[a] * b_s + (1 - q_ast_s) * μr - (γ_s / 2) * ((1 - q_ast_s)^2) * σ2r - markup_s)
			# store results
			Capacity_sim[s, a] = Capacity_s
            c_sim[s, a] = c_s
            c_sim_undiscount[s, a] = c_s_undiscount
        end
    end

	# Joint linear projections under equal-auction, participant, and observed
	# participant-total-capacity weighting.
	joint_lp, joint_summary = joint_distribution_outputs(
		γ_sim, c_sim, Capacity_sim, a_sb_list_local)
	joint_lp_undiscount, joint_summary_undiscount = joint_distribution_outputs(
		γ_sim, c_sim_undiscount, Capacity_sim, a_sb_list_local)
	# Preserve the existing equal-auction cost-projection output.
	βc_rows = (length(projection_terms) + 1):(2 * length(projection_terms))
	var_c_row = length(joint_projection_parameters)
	c_lp = [joint_lp[βc_rows, 1]; joint_lp[var_c_row, 1]]
	c_lp_undiscount = [joint_lp_undiscount[βc_rows, 1];
		joint_lp_undiscount[var_c_row, 1]]

    # summary statistics of simulated (γ,c) distribution for each auction
    logγ_mean = [mean(skipmissing(log.(γ_sim[:, a]))) for a in a_sb_list_local]
    logγ_var = [var(skipmissing(log.(γ_sim[:, a]))) for a in a_sb_list_local]
    γ_mean = [mean(skipmissing(γ_sim[:, a])) for a in a_sb_list_local]
    γ_var = [var(skipmissing(γ_sim[:, a])) for a in a_sb_list_local]
    c_mean = [mean(skipmissing(c_sim[:, a])) for a in a_sb_list_local]
    c_var = [var(skipmissing(c_sim[:, a])) for a in a_sb_list_local]
    valid_γc = [.!ismissing.(γ_sim[:, a]) .& .!ismissing.(c_sim[:, a]) for a in a_sb_list_local]
    logγ_c_cor = [safe_cor(log.(Float64.(γ_sim[:, a][valid_γc[k]])),
		Float64.(c_sim[:, a][valid_γc[k]])) for (k, a) in enumerate(a_sb_list_local)]
    γ_c_cor = [safe_cor(Float64.(γ_sim[:, a][valid_γc[k]]),
		Float64.(c_sim[:, a][valid_γc[k]])) for (k, a) in enumerate(a_sb_list_local)]

	# output
	return(
		  γ_f = γ_f
		, q_ast_f = q_ast_f
		, rp_f = rp_f
		, markup_f = markup_f
		, c_f = c_f
		, c_lp = c_lp
		, c_lp_undiscount = c_lp_undiscount
		, joint_lp = joint_lp
		, joint_lp_undiscount = joint_lp_undiscount
		, joint_summary = joint_summary
		, joint_summary_undiscount = joint_summary_undiscount
		, logγ_mean = logγ_mean
		, logγ_var = logγ_var
        , γ_mean = γ_mean
		, γ_var = γ_var
		, c_mean = c_mean
		, c_var = c_var
		, logγ_c_cor = logγ_c_cor
		, γ_c_cor = γ_c_cor
		, γ_sim = γ_sim
		, Capacity_sim = Capacity_sim
		, c_sim = c_sim
		, c_sim_undiscount = c_sim_undiscount
	)
end
