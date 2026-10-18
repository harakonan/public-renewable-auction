#=
Brazil new energy auction structural model estimation
sealed-bid auction concomitant and bid-price distributions

@author: Konan Hara

@date: September 13, 2026

@software: Julia version 1.12

=#

include("prep.jl")
using JuMP
using Ipopt

tol = 1e-6
max_iter = 10000
σ2b_init = 10.0
σ2b_lb = 1.0
σ2b_ub = 100.0
σ2K_init = 1.0
σ2K_lb = 1e-4
σ2K_ub = 100.0

include("step2_dist_func.jl")

function dist_estimates(result)
	return vcat(
		result.βC_est, result.σ2C_est,
		result.βK_est, result.σ2K_est,
		result.βb_est, result.σ2b_est,
	)
end

result = func_dist(sb_mat, a_sb_list, N_sb, obs_sb, missing, σ2b_init)
v = dist_estimates(result)

params2_dist_bs = Array{Float64}(undef, bs_num, length(v))
params2_dist_bs_status = DataFrame(
	  bootstrap_index=Int[]
	, raw_status_K=String[]
	, termination_status_K=String[]
	, primal_status_K=String[]
	, objective_value_K=Float64[]
	, sigma2_K=Float64[]
	, sigma2_K_bound_hit=Bool[]
	, retried_K=Bool[]
	, raw_status_b=String[]
	, termination_status_b=String[]
	, primal_status_b=String[]
	, objective_value_b=Float64[]
	, sigma2_b=Float64[]
	, sigma2_b_bound_hit=Bool[]
	, capacity_active_auctions=Int[]
	, prior_active_auctions=Int[]
	, prior_estimable_auctions=Int[]
	, bid_active_auctions=Int[]
	, bootstrap_weight_sum=Float64[]
	, bootstrap_effective_sample_size=Float64[]
	, bootstrap_active_firms=Int[]
)

for index_bs in 1:bs_num
	seed_bs = seed_bs_list[index_bs]
	Random.seed!(seed_bs)
	weights_all, bootstrap_diagnostics = draw_global_firm_weights()
	weights_sb = weights_all[1:obs_sb]
	result_bs = func_dist(
		sb_mat, a_sb_list, N_sb, obs_sb, result.βb_est, result.σ2b_est;
		weights=weights_sb, βK_init=result.βK_est, σ2K_start=result.σ2K_est,
	)
	params2_dist_bs[index_bs, :] = dist_estimates(result_bs)
	push!(params2_dist_bs_status, (
		  bootstrap_index=index_bs
		, raw_status_K=result_bs.raw_status_K
		, termination_status_K=result_bs.termination_status_K
		, primal_status_K=result_bs.primal_status_K
		, objective_value_K=result_bs.objective_value_K
		, sigma2_K=result_bs.σ2K_est
		, sigma2_K_bound_hit=(result_bs.σ2K_est <= σ2K_lb + 1e-6 || result_bs.σ2K_est >= σ2K_ub - 1e-6)
		, retried_K=result_bs.retried_K
		, raw_status_b=result_bs.raw_status_b
		, termination_status_b=result_bs.termination_status_b
		, primal_status_b=result_bs.primal_status_b
		, objective_value_b=result_bs.objective_value_b
		, sigma2_b=result_bs.σ2b_est
		, sigma2_b_bound_hit=(result_bs.σ2b_est <= σ2b_lb + 1e-6 || result_bs.σ2b_est >= σ2b_ub - 1e-6)
		, capacity_active_auctions=result_bs.capacity_active_auctions
		, prior_active_auctions=result_bs.prior_active_auctions
		, prior_estimable_auctions=result_bs.prior_estimable_auctions
		, bid_active_auctions=result_bs.bid_active_auctions
		, bootstrap_weight_sum=bootstrap_diagnostics.bootstrap_weight_sum
		, bootstrap_effective_sample_size=bootstrap_diagnostics.bootstrap_effective_sample_size
		, bootstrap_active_firms=bootstrap_diagnostics.bootstrap_active_firms
	))
end

solved_statuses = ["LOCALLY_SOLVED", "ALMOST_LOCALLY_SOLVED"]
bootstrap_solved =
	all(in.(params2_dist_bs_status.termination_status_K, Ref(solved_statuses))) &&
	all(params2_dist_bs_status.termination_status_b .== "LOCALLY_SOLVED") &&
	!any(params2_dist_bs_status.sigma2_K_bound_hit) &&
	!any(params2_dist_bs_status.sigma2_b_bound_hit)
raw_status_bs = bootstrap_solved ? "LOCALLY_SOLVED" : "CHECK_BOOTSTRAP_STATUS"

variable = vcat(
	["beta_C_b"], [string("beta_C_a", a) for a in 1:A_sb], ["sigma2_C"],
	["beta_K_b", "beta_K_logC"], [string("beta_K_a", a) for a in 1:A_sb], ["sigma2_K"],
	[string("beta_b_a", a) for a in a_sb_list[.!ismissing.(N_sb_a)]], ["sigma2_b"],
)
CSV.write(
	joinpath(pathtodata, "intermediate", "params2_dist_bs.csv"),
	Tables.table(params2_dist_bs);
	header=variable,
)
CSV.write(
	joinpath(pathtodata, "intermediate", "params2_dist_bs_status.csv"),
	params2_dist_bs_status,
)

function finite_sd(x)
	finite_x = x[isfinite.(x)]
	return length(finite_x) > 1 ? std(finite_x) : missing
end
params2_dist_se = [finite_sd(params2_dist_bs[:, j]) for j in axes(params2_dist_bs, 2)]
params2_dist_result = (
	  variable=variable
	, estimate=v
	, se=params2_dist_se
)
CSV.write(joinpath("output", "result", "params2_dist_result.csv"), params2_dist_result)

params2_dist_status = (
	  sample=["main", "bootstrap"]
	, raw_status_K=[result.raw_status_K, raw_status_bs]
	, termination_status_K=[result.termination_status_K, raw_status_bs]
	, objective_value_K=[result.objective_value_K, missing]
	, raw_status_b=[result.raw_status_b, raw_status_bs]
	, termination_status_b=[result.termination_status_b, raw_status_bs]
	, objective_value_b=[result.objective_value_b, missing]
)
CSV.write(joinpath("output", "result", "params2_dist_status.csv"), params2_dist_status)
