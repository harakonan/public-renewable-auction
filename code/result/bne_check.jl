#=
Numerical diagnostics for the sufficient monotone-BNE conditions in
Online Appendix, Section "Equilibrium Strategy".

The script evaluates d^2 EU(b | gamma, c) / db dgamma. Sealed-bid
auctions are evaluated at func_step2-style simulated draws of
(q*, b*, gamma, c). Descending-clock auctions are evaluated on a
theoretical (gamma, c) grid at the optimal cutoff bid price implied by
CE(b | gamma, c) = 0, keeping only types whose implied winning probability
is in [0.01, 0.99].

Outputs:
  <pathtodata>/intermediate/bne_check_grid_sb.csv
  <pathtodata>/intermediate/bne_check_grid_dc.csv
  output/result/bne_check.csv

@author: Konan Hara

@date: September 14, 2026

@software: Julia version 1.12

=#

include(joinpath(@__DIR__, "..", "estimation", "prep.jl"))

# =============================================================================
# Setup
# =============================================================================

mkpath(joinpath("output", "result"))
mkpath(joinpath(pathtodata, "intermediate"))

const delta = δ
const delta_tilde_sb_a = δtilde_sb_a
const delta_tilde_dc_a = δtilde_dc_a
const sigma2r_sb_a = σ2r_sb_a
const sigma2r_dc_a = σ2r_dc_a

const DC_GAMMA_GRID = collect(range(0.05, 1.0; length=39))
const DC_C_GRID = collect(range(10.0, 50.0; length=33))
const FD_EPS_SCALE = 1e-5
const DC_W_LB = 0.01
const DC_W_UB = 0.99

params_specified = CSV.read(
	joinpath("output", "result", "params_specified.csv"),
	DataFrame,
)
params1_hat = CSV.read(
	joinpath("output", "result", "params1_hat.csv"),
	DataFrame,
)
params2_dist_result = CSV.read(
	joinpath("output", "result", "params2_dist_result.csv"),
	DataFrame,
)
gamma_c_dist = CSV.read(
	joinpath("output", "result", "gamma_c_dist.csv"),
	DataFrame,
)
winprob_grid_raw = CSV.read(
	joinpath(pathtodata, "intermediate", "winprob_grid.csv"),
	DataFrame,
)

function truncate_winprob_grid(
	winprob_grid_raw::DataFrame,
)
	rows = DataFrame[]
	for sub in groupby(winprob_grid_raw, :auction_id)
		subdf = DataFrame(sub)
		fallback_values = unique(subdf.use_asymptote_fallback)
		@assert length(fallback_values) == 1
		use_asymptote_fallback = only(fallback_values)
		lower_keep = if use_asymptote_fallback
			winprob_inf_values = unique(subdf.winprob_inf)
			@assert length(winprob_inf_values) == 1
			subdf.winprob .- only(winprob_inf_values) .> W_asymptote_tol
		else
			subdf.winprob .> W_lb
		end
		keep = (subdf.winprob .< W_ub) .& lower_keep
		push!(rows, subdf[keep, :])
	end
	return vcat(rows...)
end

winprob_grid = truncate_winprob_grid(winprob_grid_raw)

alpha_r_hat = Float64.(
	filter(row -> occursin(r"alpha_r", row.variable), params1_hat)[:, :hetero],
)
alpha_gamma_hat = Float64.(
	filter(
		row -> occursin(r"alpha_gamma", row.variable),
		params1_hat,
	)[:, :hetero],
)
mu_r_sb(a::Int) = alpha_r_hat[1] + alpha_r_hat[2] * m_sb_a[a]
mu_r_dc(a::Int) = alpha_r_hat[1] + alpha_r_hat[2] * m_dc_a[a]

params2_dist_hat = DataFrame(
	variable=params2_dist_result[:, :variable],
	estimate=Float64.(params2_dist_result[:, :estimate]),
)
beta_C_hat = Float64.(
	filter(
		row -> occursin(r"beta_C", row.variable),
		params2_dist_hat,
	)[:, :estimate],
)
sigma2_C_hat = Float64(
	filter(row -> row.variable == "sigma2_C", params2_dist_hat)[:, :estimate][1],
)
beta_b_hat = Array{Union{Missing,Float64}}(missing, A_sb)
beta_b_hat[.!ismissing.(N_sb_a)] = Float64.(
	filter(
		row -> occursin(r"beta_b", row.variable),
		params2_dist_hat,
	)[:, :estimate],
)
sigma2_b_hat = Float64(
	filter(row -> row.variable == "sigma2_b", params2_dist_hat)[:, :estimate][1],
)
beta_K_hat = Float64.([
	filter(row -> row.variable == "beta_K_b", params2_dist_hat)[:, :estimate][1];
	filter(row -> row.variable == "beta_K_logC", params2_dist_hat)[:, :estimate][1];
	filter(
		row -> occursin(r"beta_K_a", row.variable),
		params2_dist_hat,
	)[:, :estimate]
])
sigma2_K_hat = Float64(
	filter(row -> row.variable == "sigma2_K", params2_dist_hat)[:, :estimate][1],
)

mu_p_hat = Float64.(
	filter(row -> occursin(r"mu_p_", row.variable), params_specified)[:, :value],
)
sigma2_p_hat = Float64(
	filter(row -> row.variable == "sigma2_p", params_specified)[:, :value][1],
)
sigma_p_hat = sqrt(sigma2_p_hat)

function c_np_mean_sb(a::Int)
	return Float64(
		filter(row -> row.auction_id == a, gamma_c_dist)[:, :c_mean][1],
	)
end

function c_np_sd_sb(a::Int)
	return sqrt(Float64(
		filter(row -> row.auction_id == a, gamma_c_dist)[:, :c_var][1],
	))
end

# =============================================================================
# Common Helpers
# =============================================================================

# CARA utility derivatives used in both auction formats.
u_x_cara(x, gamma) = exp(-gamma * x)
u_xx_cara(x, gamma) = -gamma * exp(-gamma * x)
u_gamma_cara(x, gamma) =
	(exp(-gamma * x) * (1.0 + gamma * x) - 1.0) / gamma^2
u_x_gamma_cara(x, gamma) = -x * exp(-gamma * x)

normal = Normal()
small_prob = 1e-14

function inv_mills(z)
	return pdf(normal, z) / max(cdf(normal, z), small_prob)
end

function inv_mills_prime(z)
	lambda = inv_mills(z)
	return -lambda * (z + lambda)
end

function safe_linear_itp(x, y)
	return linear_interpolation(x, y, extrapolation_bc=Interpolations.Flat())
end

function find_center_on_grid(f, b_grid)
	vals = [f(b) for b in b_grid]
	finite = isfinite.(vals)
	if !any(finite)
		return NaN
	end

	for k in 1:(length(b_grid) - 1)
		if finite[k] && finite[k + 1]
			if vals[k] == 0.0
				return b_grid[k]
			elseif vals[k] * vals[k + 1] < 0.0
				try
					return find_zero(f, (b_grid[k], b_grid[k + 1]))
				catch
					# Fall through to nearest-grid fallback.
				end
			end
		end
	end

	idx_finite = findall(finite)
	idx = idx_finite[argmin(abs.(vals[idx_finite]))]
	return b_grid[idx]
end

# =============================================================================
# Sealed-Bid Diagnostic
# =============================================================================

function build_sb_winprob_functions(a::Int)
	sub = winprob_grid[winprob_grid.auction_id .== a, :]
	sort!(sub, :bid_price)

	b_grid = Float64.(sub[:, :bid_price])
	w_grid = Float64.(sub[:, :winprob])

	dw_grid = diff(w_grid) ./ diff(b_grid)
	db_grid = b_grid[1:(end - 1)] .+ diff(b_grid) ./ 2.0

	return (
		W=safe_linear_itp(b_grid, w_grid),
		dW=safe_linear_itp(db_grid, dw_grid),
		b_min=minimum(b_grid),
		b_max=maximum(b_grid),
	)
end

function qstar_sb_bounds(b, gamma, a, q_lb, q_ub)
	mu = mu_r_sb(a)
	sigma2_r = sigma2r_sb_a[a]
	delta_tilde = delta_tilde_sb_a[a]
	delta_l = delta^l_sb_a[a]

	q_tilde = 1.0 - (delta_l * mu - delta_tilde * b) / (gamma * sigma2_r)
	q = clamp(q_tilde, q_lb, q_ub)
	status = q_tilde < q_lb ? "lower" :
		(q_tilde > q_ub ? "upper" : "interior")
	dq_dgamma = status == "interior" ?
		(delta_l * mu - delta_tilde * b) / (gamma^2 * sigma2_r) : 0.0
	return q, dq_dgamma, status, q_tilde
end

function ce_sb_bounds(b, gamma, c, a, q_lb, q_ub)
	q, _, _, _ = qstar_sb_bounds(b, gamma, a, q_lb, q_ub)
	mu = mu_r_sb(a)
	sigma2_r = sigma2r_sb_a[a]
	delta_tilde = delta_tilde_sb_a[a]
	rp = 0.5 * gamma * (1.0 - q)^2 * sigma2_r
	return q * delta_tilde * b + (1.0 - q) * mu - c - rp
end

function markup_sb_bounds(b, gamma, a, winprob, q_lb, q_ub)
	q, _, _, _ = qstar_sb_bounds(b, gamma, a, q_lb, q_ub)
	w = winprob.W(b)
	dw = winprob.dW(b)
	inside = 1.0 - (w * gamma * q * delta_tilde_sb_a[a]) / dw
	return inside > 0.0 ? log(inside) / gamma : NaN
end

function rhs_sb_bounds(b, gamma, c, a, winprob, q_lb, q_ub)
	q, dq_dgamma, status, q_tilde = qstar_sb_bounds(b, gamma, a, q_lb, q_ub)
	ce = ce_sb_bounds(b, gamma, c, a, q_lb, q_ub)
	w = winprob.W(b)
	dw = winprob.dW(b)
	delta_tilde = delta_tilde_sb_a[a]
	sigma2_r = sigma2r_sb_a[a]

	rp_gamma = 0.5 * (1.0 - q)^2 * sigma2_r

	ux = u_x_cara(ce, gamma)
	uxx = u_xx_cara(ce, gamma)
	ugamma = u_gamma_cara(ce, gamma)
	uxgamma = u_x_gamma_cara(ce, gamma)

	term_wprime = dw * (-ux * rp_gamma + ugamma)
	term_rp = -w * uxx * rp_gamma * q * delta_tilde
	term_q = w * (uxgamma * q * delta_tilde + ux * dq_dgamma * delta_tilde)
	rhs = term_wprime + term_rp + term_q

	return (
		rhs=rhs,
		q=q,
		q_tilde=q_tilde,
		dq_dgamma=dq_dgamma,
		q_status=status,
		ce=ce,
		W=w,
		dW=dw,
		term_wprime=term_wprime,
		term_rp=term_rp,
		term_q=term_q,
	)
end

function simulated_K_log_sb(b, capacity, a::Int, Z_K)
	if structural_zero_prior_sb_a[a]
		return 0.0
	end
	mu_K = beta_K_hat[2 + a] +
		beta_K_hat[1] * b +
		beta_K_hat[2] * log(capacity)
	return max(0.0, mu_K + sqrt(sigma2_K_hat) * Z_K)
end

function simulated_sb_draw(a::Int, s::Int, winprob, Ugamma_c, U_K_dist)
	b = Float64(beta_b_hat[a]) + Ugamma_c[1, s] * sqrt(sigma2_b_hat)
	capacity = exp(
		beta_C_hat[1] * b +
		beta_C_hat[1 + a] +
		Ugamma_c[2, s] * sqrt(sigma2_C_hat)
	)
	K_log = simulated_K_log_sb(b, capacity, a, U_K_dist[s])
	gamma = exp(alpha_gamma_hat[1] + alpha_gamma_hat[2] * K_log)

	q_lb = (
		ceil(capacity * qlb_sb_a[a] / Qmu_sb_a[a] - eps(1.0)) *
		Qmu_sb_a[a]
	) / capacity
	q_ub = (
		floor(capacity / Qmu_sb_a[a] + eps(1.0)) *
		Qmu_sb_a[a]
	) / capacity
	q, dq_dgamma, q_status, q_tilde = qstar_sb_bounds(b, gamma, a, q_lb, q_ub)

	markup = markup_sb_bounds(b, gamma, a, winprob, q_lb, q_ub)
	mu = mu_r_sb(a)
	sigma2_r = sigma2r_sb_a[a]
	c = q * delta_tilde_sb_a[a] * b +
		(1.0 - q) * mu -
		0.5 * gamma * (1.0 - q)^2 * sigma2_r -
		markup

	r = rhs_sb_bounds(b, gamma, c, a, winprob, q_lb, q_ub)

	return (
		b=b,
		capacity=capacity,
		K_log=K_log,
		prior_started_wind_capacity_100=exp(K_log) - 1.0,
		gamma=gamma,
		loggamma=log(gamma),
		c=c,
		q=q,
		q_lb=q_lb,
		q_ub=q_ub,
		q_tilde=q_tilde,
		dq_dgamma=dq_dgamma,
		q_status=q_status,
		markup=markup,
		rhs=r.rhs,
		ce=r.ce,
		W=r.W,
		dW=r.dW,
		term_wprime=r.term_wprime,
		term_rp=r.term_rp,
		term_q=r.term_q,
	)
end

# =============================================================================
# Descending-Clock Diagnostic
# =============================================================================

function k_normal(mu, sigma2, tau)
	return mu * tau + 0.5 * sigma2 * tau^2
end

function k_trunc_lower_normal(mu, sigma, b, tau)
	z0 = (mu - b) / sigma
	z1 = z0 + sigma * tau
	return (
		mu * tau +
		0.5 * sigma^2 * tau^2 +
		log(max(cdf(normal, z1), small_prob)) -
		log(max(cdf(normal, z0), small_prob))
	)
end

function mu_tilde_dc(q, b, gamma, a)
	z = (
		mu_p_hat[a] -
		q * gamma * delta_tilde_dc_a[a] * sigma2_p_hat -
		b
	) / sigma_p_hat
	return mu_p_hat[a] + sigma_p_hat * inv_mills(z)
end

function q_foc_residual_dc(q, b, gamma, a)
	mu = mu_r_dc(a)
	sigma2_r = sigma2r_dc_a[a]
	delta_tilde = delta_tilde_dc_a[a]
	delta_l = delta^l_dc_a[a]
	mu_tilde = mu_tilde_dc(q, b, gamma, a)
	rhs = (1.0 / (1.0 + delta_tilde^2 * sigma2_p_hat / sigma2_r)) *
		(1.0 - (delta_l * mu - delta_tilde * mu_tilde) / (gamma * sigma2_r))
	return rhs - q
end

function qstar_dc(b, gamma, a)
	q_unclamped = try
		solve(ZeroProblem(q -> q_foc_residual_dc(q, b, gamma, a), x0), Order0())
	catch
		root_grid = collect(range(-0.5, 1.5; length=401))
		find_center_on_grid(q -> q_foc_residual_dc(q, b, gamma, a), root_grid)
	end

	q = clamp(q_unclamped, 0.0, 1.0)
	status = q_unclamped < 0.0 ? "lower" :
		(q_unclamped > 1.0 ? "upper" : "interior")

	mu = mu_r_dc(a)
	sigma2_r = sigma2r_dc_a[a]
	delta_tilde = delta_tilde_dc_a[a]
	delta_l = delta^l_dc_a[a]
	z = (mu_p_hat[a] - b) / sigma_p_hat - q * gamma * delta_tilde * sigma_p_hat
	lambda_prime = inv_mills_prime(z)
	mu_tilde = mu_p_hat[a] + sigma_p_hat * inv_mills(z)
	denom = gamma^2 * (
		sigma2_r +
		(1.0 + lambda_prime) * delta_tilde^2 * sigma2_p_hat
	)
	dq_dgamma_unclamped = (
		delta_l * mu -
		delta_tilde * mu_tilde -
		q * gamma * delta_tilde^2 * sigma2_p_hat * lambda_prime
	) / denom
	dq_dgamma = status == "interior" ? dq_dgamma_unclamped : 0.0

	return q, dq_dgamma, status, q_unclamped
end

function ce_dc(b, gamma, c, a)
	q, _, _, _ = qstar_dc(b, gamma, a)
	tau_p = -gamma * q * delta_tilde_dc_a[a]
	tau_r = -gamma * (1.0 - q)
	mu = delta^l_dc_a[a] * mu_r_dc(a)
	kp = k_trunc_lower_normal(mu_p_hat[a], sigma_p_hat, b, tau_p)
	kr = k_normal(mu, sigma2r_dc_a[a], tau_r)
	return -(kp + kr) / gamma - c
end

function ce_gamma_fd_dc(b, gamma, c, a)
	eps = max(FD_EPS_SCALE, FD_EPS_SCALE * gamma)
	lo = max(gamma - eps, gamma / 2.0)
	hi = gamma + eps
	return (ce_dc(b, hi, c, a) - ce_dc(b, lo, c, a)) / (hi - lo)
end

function rhs_dc(b, gamma, c, a)
	q, dq_dgamma, status, q_unclamped = qstar_dc(b, gamma, a)
	ce = ce_dc(b, gamma, c, a)
	ce_gamma = ce_gamma_fd_dc(b, gamma, c, a)

	mu_p = mu_p_hat[a]
	sigma_p = sigma_p_hat
	delta_tilde = delta_tilde_dc_a[a]
	z0 = (mu_p - b) / sigma_p
	z1 = z0 - q * gamma * delta_tilde * sigma_p
	delta_lambda = (inv_mills(z0) - inv_mills(z1)) / sigma_p
	lambda_prime = inv_mills_prime(z1)

	ux = u_x_cara(ce, gamma)
	uxx = u_xx_cara(ce, gamma)
	ugamma = u_gamma_cara(ce, gamma)
	uxgamma = u_x_gamma_cara(ce, gamma)

	fp = pdf(Normal(mu_p, sigma_p), b)
	one_minus_F = 1.0 - cdf(Normal(mu_p, sigma_p), b)

	first = -fp * (ux * ce_gamma + ugamma)
	curly = gamma^(-2) * ux - gamma^(-1) * (uxx * ce_gamma + uxgamma)
	not_too_risk_averse = 1.0 + gamma * ce + gamma^2 * ce_gamma
	last = -gamma^(-1) * ux * (1.0 / sigma_p) * lambda_prime *
		(gamma * dq_dgamma + q) * delta_tilde * sigma_p
	bracket = delta_lambda * curly + last
	rhs = first + one_minus_F * bracket

	return (
		rhs=rhs,
		q=q,
		q_tilde=q_unclamped,
		dq_dgamma=dq_dgamma,
		q_status=status,
		ce=ce,
		ce_gamma=ce_gamma,
		first=first,
		curly=curly,
		not_too_risk_averse=not_too_risk_averse,
		delta_lambda=delta_lambda,
		lambda_prime=lambda_prime,
		last=last,
		bracket=bracket,
		fp=fp,
		one_minus_F=one_minus_F,
	)
end

function optimal_b_dc(gamma, c, a)
	# Descending-clock winning probability is 1 - F_p(b). The theoretical
	# optimal cutoff bid price solves CE(b | gamma, c) = 0.
	p_dist = Normal(mu_p_hat[a], sigma_p_hat)
	b_min = max(0.0, quantile(p_dist, 1.0 - DC_W_UB))
	b_max = quantile(p_dist, 1.0 - DC_W_LB)
	f(b) = ce_dc(b, gamma, c, a)

	f_min = f(b_min)
	f_max = f(b_max)
	if !isfinite(f_min) || !isfinite(f_max)
		return missing, missing, "no_finite_bracket"
	elseif f_min == 0.0
		w = 1.0 - cdf(p_dist, b_min)
		return b_min, w, "root_lower"
	elseif f_max == 0.0
		w = 1.0 - cdf(p_dist, b_max)
		return b_max, w, "root_upper"
	elseif f_min * f_max > 0.0
		return missing, missing, "outside_winprob_bounds"
	end

	try
		b = find_zero(f, (b_min, b_max))
		w = 1.0 - cdf(p_dist, b)
		return b, w, "root"
	catch
		return missing, missing, "root_failed"
	end
end

# =============================================================================
# Output
# =============================================================================

function summarize_by_auction(df::DataFrame, format::String)
	rows = NamedTuple[]
	for subdf in groupby(df, :auction_id)
		neg = subdf[subdf.rhs .< 0.0, :]
		pos = subdf[subdf.rhs .> 0.0, :]
		violation = format == "sealed_bid" ? pos : neg
		worst_idx = argmin(subdf.rhs)
		best_idx = argmax(subdf.rhs)
		expected_sign = format == "sealed_bid" ? "nonpositive" : "nonnegative"
		violation_sign = format == "sealed_bid" ? "positive" : "negative"
		c_mean_np = format == "sealed_bid" ?
			c_np_mean_sb(subdf.auction_id[1]) : missing
		c_sd_np = format == "sealed_bid" ?
			c_np_sd_sb(subdf.auction_id[1]) : missing
		min_c_pos = nrow(pos) == 0 ? missing : minimum(pos.c)
		max_c_pos = nrow(pos) == 0 ? missing : maximum(pos.c)
		c_gap_min_pos = (
			ismissing(c_mean_np) || ismissing(min_c_pos) ?
			missing : c_mean_np - min_c_pos
		)
		c_gap_max_pos = (
			ismissing(c_mean_np) || ismissing(max_c_pos) ?
			missing : c_mean_np - max_c_pos
		)
		c_gap_min_pos_sd = (
			ismissing(c_gap_min_pos) || ismissing(c_sd_np) ?
			missing : c_gap_min_pos / c_sd_np
		)
		c_gap_max_pos_sd = (
			ismissing(c_gap_max_pos) || ismissing(c_sd_np) ?
			missing : c_gap_max_pos / c_sd_np
		)
		push!(rows, (
			format=format,
			auction_id=subdf.auction_id[1],
			n_grid=nrow(subdf),
			expected_sign=expected_sign,
			violation_sign=violation_sign,
			n_violation=nrow(violation),
			share_violation=nrow(violation) / nrow(subdf),
			n_negative=nrow(neg),
			share_negative=nrow(neg) / nrow(subdf),
			n_positive=nrow(pos),
			share_positive=nrow(pos) / nrow(subdf),
			min_rhs=subdf.rhs[worst_idx],
			max_rhs=subdf.rhs[best_idx],
			min_gamma_positive=nrow(pos) == 0 ? missing : minimum(pos.gamma),
			max_gamma_positive=nrow(pos) == 0 ? missing : maximum(pos.gamma),
			min_c_positive=min_c_pos,
			max_c_positive=max_c_pos,
			min_b_positive=nrow(pos) == 0 ? missing : minimum(pos.b),
			max_b_positive=nrow(pos) == 0 ? missing : maximum(pos.b),
			min_q_positive=nrow(pos) == 0 ? missing : minimum(pos.q),
			max_q_positive=nrow(pos) == 0 ? missing : maximum(pos.q),
			c_np_mean=c_mean_np,
			c_np_sd=c_sd_np,
			c_np_mean_minus_min_c_positive=c_gap_min_pos,
			c_np_mean_minus_max_c_positive=c_gap_max_pos,
			c_np_mean_minus_min_c_positive_sd=c_gap_min_pos_sd,
			c_np_mean_minus_max_c_positive_sd=c_gap_max_pos_sd,
		))
	end
	return DataFrame(rows)
end

# Sealed-bid draws follow the same simulated type recovery as func_step2.
println("Evaluating sealed-bid auctions...")
Random.seed!(seed_γc)
U_gamma_c = rand(Normal(), (2, Sγc))
Random.seed!(seed_γc + 10_000)
U_K_dist = rand(Normal(), Sγc)
sb_rows = NamedTuple[]
for a in a_sb_list[.!ismissing.(N_sb_a)]
	winprob = build_sb_winprob_functions(a)
	for s in 1:Sγc
		draw = simulated_sb_draw(a, s, winprob, U_gamma_c, U_K_dist)
		push!(sb_rows, (
			format="sealed_bid",
			auction_id=a,
			draw_id=s,
			gamma=draw.gamma,
			loggamma=draw.loggamma,
			c=draw.c,
			b=draw.b,
			capacity=draw.capacity,
			K_log=draw.K_log,
			prior_started_wind_capacity_100=
				draw.prior_started_wind_capacity_100,
			q=draw.q,
			q_lb=draw.q_lb,
			q_ub=draw.q_ub,
			q_tilde=draw.q_tilde,
			dq_dgamma=draw.dq_dgamma,
			q_status=draw.q_status,
			markup=draw.markup,
			ce=draw.ce,
			rhs=draw.rhs,
			negative=draw.rhs < 0.0,
			W=draw.W,
			dW=draw.dW,
			term_wprime=draw.term_wprime,
			term_rp=draw.term_rp,
			term_q=draw.term_q,
		))
	end
end
sb_detail = DataFrame(sb_rows)

# Type grid points outside the win-probability band are omitted.
println("Evaluating descending-clock auctions...")
dc_rows = NamedTuple[]
for a in 1:A_dc
	global_auction_id = A_sb + a
	for gamma in DC_GAMMA_GRID, c in DC_C_GRID
		b, W, b_source = optimal_b_dc(gamma, c, a)
		if !ismissing(b) && isfinite(b)
			r = rhs_dc(b, gamma, c, a)
			push!(dc_rows, (
				format="descending_clock",
				auction_id=global_auction_id,
				auction_index=a,
				gamma=gamma,
				c=c,
				b=b,
				W=W,
				b_center=b,
				b_source=b_source,
				q=r.q,
				q_tilde=r.q_tilde,
				dq_dgamma=r.dq_dgamma,
				q_status=r.q_status,
				ce=r.ce,
				ce_gamma=r.ce_gamma,
				rhs=r.rhs,
				negative=r.rhs < 0.0,
				first=r.first,
				curly=r.curly,
				not_too_risk_averse=r.not_too_risk_averse,
				delta_lambda=r.delta_lambda,
				lambda_prime=r.lambda_prime,
				last=r.last,
				bracket=r.bracket,
				fp=r.fp,
				one_minus_F=r.one_minus_F,
			))
		end
	end
end
dc_detail = DataFrame(dc_rows)

auction_summary = vcat(
	summarize_by_auction(sb_detail, "sealed_bid"),
	summarize_by_auction(dc_detail, "descending_clock");
	cols=:union,
)

CSV.write(joinpath("output", "result", "bne_check.csv"), auction_summary)

CSV.write(
	joinpath(pathtodata, "intermediate", "bne_check_grid_sb.csv"),
	sb_detail,
)
CSV.write(
	joinpath(pathtodata, "intermediate", "bne_check_grid_dc.csv"),
	dc_detail,
)

println(auction_summary)
