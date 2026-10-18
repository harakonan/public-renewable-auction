#=
Brazil new energy auction structural model estimation
first stage
estimation results

@author: Konan Hara

@date: September 13, 2026

@software: Julia version 1.12

=#

# Load the preferred Step 1 data construction and function definitions without
# rerunning the estimator.
ENV["STEP1_DEFINITIONS_ONLY"] = "true"
include("step1.jl")

# load first stage parameter estimates
params1_hat = CSV.read(joinpath("output", "result", "params1_hat.csv"), DataFrame)
# variable names
variable = filter(row -> occursin(r"alpha_gamma|alpha_r|alpha_eta", row.variable), params1_hat)[:, :variable]

# K_log=log(1+K/100), where K is prior-started wind capacity in MW.
K_log_sb = sb_mat[:, 3]
K_log_dc = dc_mat[:, 3]

# risk aversion distribution
# [average, standard deviation] of γ
# homogeneous risk aversion
αγ_homo = filter(row -> row.variable == "alpha_gamma0", params1_hat)[:, :homo_gamma][1]
γ_dist_homo = [exp(αγ_homo), missing]
# heterogeneous risk aversion
αγ = filter(row -> occursin(r"alpha_gamma", row.variable), params1_hat)[:, :hetero]
# γ for each observation
γ_sb = [exp(αγ[1] + αγ[2] * K_log_sb[i]) for i = 1:obs_sb]
γ_dc = [exp(αγ[1] + αγ[2] * K_log_dc[i]) for i = 1:obs_dc]
γ_dist_hetero = [(sum(γ_sb) + sum(γ_dc)) / obs, std(vcat(γ_sb, γ_dc))]

# collect bootstrap estimates
params1_homo_bs = Array{Float64}(undef, bs_num, length(variable) - 1)
params1_hetero_bs = Array{Float64}(undef, bs_num, length(variable))
γ_dist_homo_bs = Array{Float64}(undef, bs_num)
γ_dist_hetero_bs = Array{Float64}(undef, bs_num, 2)
for index_bs in 1:bs_num
    status_path = joinpath(
		pathtodata, "intermediate", "params1_bs_status", string("bs", index_bs, ".csv")
	)
	isfile(status_path) || error("Missing bootstrap status file: $status_path")
	status_b = CSV.read(status_path, DataFrame)
	all(status_b.termination_status .== "LOCALLY_SOLVED") ||
		error("Bootstrap draw $index_bs contains an unsuccessful optimization")

    # load bootstrap estimates
	result_path_b = joinpath(
		pathtodata, "intermediate", "params1_bs", string("bs", index_bs, ".csv")
	)
	isfile(result_path_b) || error("Missing bootstrap result file: $result_path_b")
    params1_bs_b = CSV.read(result_path_b, DataFrame)
    params1_homo_bs[index_bs, :] = collect(skipmissing(params1_bs_b[1:end-1, :homo_gamma]))
    params1_hetero_bs[index_bs, :] = params1_bs_b[1:end-1, :hetero]

    # risk aversion distribution
    # [average, standard deviation] of γ
    # homogeneous risk aversion
	αγ_homo_b = filter(row -> row.variable == "alpha_gamma0", params1_bs_b)[:, :homo_gamma][1]
	γ_dist_homo_bs[index_bs] = exp(αγ_homo_b)
    # heterogeneous risk aversion
    αγ_b = filter(row -> occursin(r"alpha_gamma", row.variable), params1_bs_b)[:, :hetero]
    # γ for each observation
	γ_sb_b = [exp(αγ_b[1] + αγ_b[2] * K_log_sb[i]) for i = 1:obs_sb]
	γ_dc_b = [exp(αγ_b[1] + αγ_b[2] * K_log_dc[i]) for i = 1:obs_dc]
    γ_dist_hetero_bs[index_bs, :] = [
		(sum(γ_sb_b) + sum(γ_dc_b)) / obs
		, std(vcat(γ_sb_b, γ_dc_b))
	]
end

# save results
# bootstrap standard errors
vcm = cov(params1_homo_bs)
params1_homo_se = sqrt.([vcm[i, i] for i in 1:size(vcm)[1]])
params1_homo_se = [params1_homo_se[1]; missing; params1_homo_se[2:end]]
γ_dist_homo_se = [sqrt(var(γ_dist_homo_bs)), missing]
vcm = cov(params1_hetero_bs)
params1_hetero_se = sqrt.([vcm[i, i] for i in 1:size(vcm)[1]])
vcm = cov(γ_dist_hetero_bs)
γ_dist_hetero_se = sqrt.([vcm[i, i] for i in 1:size(vcm)[1]])
# result table
params1_homo = filter(row -> occursin(r"alpha_gamma|alpha_r|alpha_eta", row.variable), params1_hat)[:, :homo_gamma]
params1_hetero = filter(row -> occursin(r"alpha_gamma|alpha_r|alpha_eta", row.variable), params1_hat)[:, :hetero]
params1_result = (
	  variable = [variable; "gamma_mean"; "gamma_sd"]  
    , homo= [params1_homo; γ_dist_homo]
    , homo_se=[params1_homo_se; γ_dist_homo_se]
    , hetero= [params1_hetero; γ_dist_hetero]
    , hetero_se=[params1_hetero_se; γ_dist_hetero_se]
)
CSV.write(joinpath("output", "result", "params1_result.csv")
	, params1_result)
