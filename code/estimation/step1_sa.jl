#=
Brazil new energy auction structural model estimation
first stage
sensitivity analysis

@author: Konan Hara

@date: September 13, 2026

@software: Julia version 1.12
=#

# Reuse the production moment and optimization implementation without running it.
ENV["STEP1_DEFINITIONS_ONLY"] = "true"
include("step1.jl")

# Sensitivity analysis index
# SLURM_ARRAY_TASK_ID maps directly to the sensitivity-analysis index.
index_sa = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))

# Inferred latent DC winning-bid distribution, constructed from mu_pD and
# sigma2_pD. Keep it fixed when varying mu_p and sigma2_p below: those
# sensitivities isolate clearing-price beliefs in the portfolio-choice FOC.
inv_cdf_b_dc = deserialize(joinpath(pathtodata, "intermediate", "inv_cdf_b_dc.dat"))

# load specified parameters
params_specified = CSV.read(joinpath("output", "result", "params_specified.csv"), DataFrame)

# load first stage parameter estimates
params1_hat = CSV.read(joinpath("output", "result", "params1_hat.csv"), DataFrame)

# unpack original specified parameters
μphat_original = filter(row -> occursin(r"mu_p_", row.variable), params_specified)[:, :value]
σ2phat_original = filter(row -> row.variable == "sigma2_p", params_specified)[:, :value][1]
spot_variance_original = copy(spot_variance)
δ_original = δ
spot_result = CSV.read(joinpath("output", "result", "spot_result.csv"), DataFrame)
ρhat = filter(row -> row.variable == "rho", spot_result)[:, :est][1]
σ2ξhat_original = filter(row -> row.variable == "var_xi", spot_result)[:, :est][1]
auction_lead_time = unique(auction_data[:, [:auction_id, :lead_time]])

# sensitivity analysis grids
sa_grid_sd = [0.50, 0.75, 1.0, 1.25, 1.50]
sa_grid_mean = [-2.0, -1.0, 0.0, 1.0, 2.0]
sa_grid_delta = [0.90, 0.95, 0.99]
sa_num = 2 * length(sa_grid_sd) + length(sa_grid_mean) + length(sa_grid_delta)

if !(1 <= index_sa <= sa_num)
	error("Sensitivity analysis index must be between 1 and $sa_num. Got $index_sa.")
end

# calculate discounted spot market variance for a given discount factor
function func_spot_variance(ρ, σ2ξ, δ_local, T_local, l)
	l_fc = [floor(l), ceil(l)]
	l_plus_T = [Int(l_fc[k] + T_local) for k in 1:2]
	var_prep = zeros(2)
	for k in 1:2
		for t in 1:(l_plus_T[k] - 1)
			if t <= l_fc[k]
				var_prep[k] += (((δ_local^(l_fc[k])) * (ρ^(l_fc[k] - t)) * ((1 - (δ_local^T_local) * (ρ^T_local)) / (1 - δ_local * ρ)))^2)
			else
				var_prep[k] += (((δ_local^t) * (1 - (δ_local^(l_plus_T[k] - t)) * (ρ^(l_plus_T[k] - t))) / (1 - δ_local * ρ))^2)
			end
		end
	end
	spot_variance_local = l == l_fc[1] ? (σ2ξ / (T_local^2)) * var_prep[1] :
		(σ2ξ / (T_local^2)) * ((l_fc[2] - l) * var_prep[1] + (l - l_fc[1]) * var_prep[2])
	return spot_variance_local
end

# default sensitivity-analysis parameter values
sa_type = ""
sa_value = NaN
μphat = copy(μphat_original)
σ2phat = σ2phat_original
spot_variance = copy(spot_variance_original)

# map sensitivity-analysis index to sensitivity type and grid index
function parse_sa_index(index_sa)
	if index_sa <= length(sa_grid_sd)
		return "sigma2_xi", index_sa
	elseif index_sa <= length(sa_grid_sd) + length(sa_grid_mean)
		return "mu_p", index_sa - length(sa_grid_sd)
	elseif index_sa <= 2 * length(sa_grid_sd) + length(sa_grid_mean)
		return "sigma2_p", index_sa - (length(sa_grid_sd) + length(sa_grid_mean))
	else
		return "delta", index_sa - (2 * length(sa_grid_sd) + length(sa_grid_mean))
	end
end

# baseline value and internal increment for continuation paths
function baseline_scenario_value(sa_type)
	if sa_type == "delta"
		return δ_original
	elseif sa_type == "mu_p"
		return 0.0
	elseif sa_type == "sigma2_xi" || sa_type == "sigma2_p"
		return 1.0
	end
	error("Unknown sensitivity-analysis type: $sa_type")
end

function continuation_increment(sa_type)
	if sa_type == "delta"
		return 0.01
	elseif sa_type == "mu_p"
		return 0.20
	elseif sa_type == "sigma2_xi" || sa_type == "sigma2_p"
		return 0.05
	end
	error("Unknown sensitivity-analysis type: $sa_type")
end

# Fine continuation path from the baseline to a reported sensitivity value.
# The baseline is excluded because the original estimate is already available.
function continuation_values(sa_type, target_value)
	base_value = baseline_scenario_value(sa_type)
	increment = continuation_increment(sa_type)
	difference = target_value - base_value
	isapprox(difference, 0.0; atol=1e-10) && return Float64[]
	n_steps = round(Int, abs(difference) / increment)
	isapprox(n_steps * increment, abs(difference); atol=1e-10) ||
		error("Target $target_value is not reachable from $base_value in increments of $increment")
	direction = sign(difference)
	values = [round(base_value + direction * increment * k; digits=10) for k in 1:n_steps]
	values[end] = Float64(target_value)
	return values
end

# value associated with a sensitivity-analysis grid point
function scenario_value(sa_type, grid_index)
	if sa_type == "sigma2_xi" || sa_type == "sigma2_p"
		return sa_grid_sd[grid_index]
	elseif sa_type == "mu_p"
		return sa_grid_mean[grid_index]
	elseif sa_type == "delta"
		return sa_grid_delta[grid_index]
	end
	error("Unknown sensitivity-analysis type: $sa_type")
end

# apply a one-at-a-time perturbation at a supplied value
function apply_scenario_value!(sa_type, sa_value_local)
	global μphat = copy(μphat_original)
	global σ2phat = σ2phat_original
	global δ = δ_original
	global δtilde_sb_a = [sum([δ^t for t = l_sb_a[a]:1.0:(l_sb_a[a] + T - 1)]) / T for a in 1:A_sb]
	global δtilde_dc_a = [sum([δ^t for t = l_dc_a[a]:1.0:(l_dc_a[a] + T - 1)]) / T for a in 1:A_dc]

	spot_variance_local = copy(spot_variance_original)
	if sa_type == "sigma2_xi"
		σ2ξhat = (sa_value_local^2) * σ2ξhat_original
		spot_variance_local = [func_spot_variance(ρhat, σ2ξhat, δ, T, auction_lead_time[a, :lead_time]) for a in 1:nrow(auction_lead_time)]
	elseif sa_type == "mu_p"
		global μphat = μphat_original .+ sa_value_local .* sqrt(σ2phat_original)
	elseif sa_type == "sigma2_p"
		global σ2phat = (sa_value_local^2) * σ2phat_original
	elseif sa_type == "delta"
		global δ = sa_value_local
		global δtilde_sb_a = [sum([δ^t for t = l_sb_a[a]:1.0:(l_sb_a[a] + T - 1)]) / T for a in 1:A_sb]
		global δtilde_dc_a = [sum([δ^t for t = l_dc_a[a]:1.0:(l_dc_a[a] + T - 1)]) / T for a in 1:A_dc]
		spot_variance_local = [func_spot_variance(ρhat, σ2ξhat_original, δ, T, auction_lead_time[a, :lead_time]) for a in 1:nrow(auction_lead_time)]
	else
		error("Unknown sensitivity-analysis type: $sa_type")
	end

	global spot_variance = spot_variance_local
	global σ2r_sb_a = spot_variance_local[a_sb_a]
	global σ2r_dc_a = spot_variance_local[a_dc_a]

	return sa_value_local
end

# apply one-at-a-time perturbation for a given sensitivity type and grid index
apply_scenario!(sa_type, grid_index) =
	apply_scenario_value!(sa_type, scenario_value(sa_type, grid_index))

sa_type, target_grid_index = parse_sa_index(index_sa)
sa_value = scenario_value(sa_type, target_grid_index)
path_values = continuation_values(sa_type, sa_value)
apply_scenario_value!(sa_type, sa_value)

# descending-clock auction bid price integral by Gauss–Legendre quadrature
# winner bid prices and weights
U_b_dc_s, w_b_dc_s = uniform_quadrature_legendre(Sgmm)
b_dc_s = Array{Float64}(undef, obs_dc, Sgmm)

println("--- Running Sensitivity Analysis $index_sa ---")
println("Sensitivity type: $sa_type")
println("Sensitivity value: $sa_value")
println("Continuation path: ", path_values)

# winner bid prices for each auction and each simulation (for numerical integration in GMM)
for i = 1:obs_dc, s = 1:Sgmm
	b_dc_s[i, s] = inv_cdf_b_dc[i](U_b_dc_s[s])
end

# load first-stage estimates for initial values
# homogeneous risk aversion
αγ_homo = Float64.(filter(row -> row.variable == "alpha_gamma0", params1_hat)[:, :homo_gamma])[1]
αr_homo = Float64.(filter(row -> occursin(r"alpha_r", row.variable), params1_hat)[:, :homo_gamma])
αη_homo = Float64.(filter(row -> row.variable == "alpha_eta", params1_hat)[:, :homo_gamma])[1]
f_homo = Float64.(filter(row -> row.variable == "objective_value", params1_hat)[:, :homo_gamma])[1]
αγ_homo_step1 = Float64.(filter(row -> row.variable == "alpha_gamma0", params1_hat)[:, :homo_gamma_step1])[1]
αr_homo_step1 = Float64.(filter(row -> occursin(r"alpha_r", row.variable), params1_hat)[:, :homo_gamma_step1])
αη_homo_step1 = Float64.(filter(row -> row.variable == "alpha_eta", params1_hat)[:, :homo_gamma_step1])[1]
αγ_hetero = Float64.(filter(row -> occursin(r"alpha_gamma", row.variable), params1_hat)[:, :hetero])
αr_hetero = Float64.(filter(row -> occursin(r"alpha_r", row.variable), params1_hat)[:, :hetero])
αη_hetero = Float64.(filter(row -> row.variable == "alpha_eta", params1_hat)[:, :hetero])[1]
f_hetero = Float64.(filter(row -> row.variable == "objective_value", params1_hat)[:, :hetero])[1]
αγ_hetero_step1 = Float64.(filter(row -> occursin(r"alpha_gamma", row.variable), params1_hat)[:, :hetero_step1])
αr_hetero_step1 = Float64.(filter(row -> occursin(r"alpha_r", row.variable), params1_hat)[:, :hetero_step1])
αη_hetero_step1 = Float64.(filter(row -> row.variable == "alpha_eta", params1_hat)[:, :hetero_step1])[1]

# IV
# homogeneous risk aversion
IV_homo = build_iv_homo(sb_mat, dc_mat)
# heterogeneous risk aversion
IV_hetero = build_iv_hetero(sb_mat, dc_mat)

# Compute the original-sample second-step weight matrices from the recorded
# original first-step estimates, then restore the target scenario.
μphat = copy(μphat_original)
σ2phat = σ2phat_original
δ = δ_original
δtilde_sb_a = [sum([δ^t for t = l_sb_a[a]:1.0:(l_sb_a[a] + T - 1)]) / T for a in 1:A_sb]
δtilde_dc_a = [sum([δ^t for t = l_dc_a[a]:1.0:(l_dc_a[a] + T - 1)]) / T for a in 1:A_dc]
spot_variance = copy(spot_variance_original)
σ2r_sb_a = spot_variance[a_sb_a]
σ2r_dc_a = spot_variance[a_dc_a]
W_opt_homo = build_W(
	  IV_homo
	, αγ_homo_step1
	, αr_homo_step1
	, αη_homo_step1;
	  b_dc_s=b_dc_s
	, w_b_dc_s=w_b_dc_s
)
W_opt_hetero = build_W(
	  IV_hetero
	, αγ_hetero_step1
	, αr_hetero_step1
	, αη_hetero_step1;
	  b_dc_s=b_dc_s
	, w_b_dc_s=w_b_dc_s
)
apply_scenario_value!(sa_type, sa_value)

function empty_sa_status()
	return DataFrame(
		  model=String[]
		, step=String[]
		, sa_type=String[]
		, sa_value=Float64[]
		, continuation_index=Union{Missing, Int}[]
		, reported_target=Bool[]
		, selected=Bool[]
		, optimizer=String[]
		, raw_status=String[]
		, termination_status=String[]
		, objective_value=Float64[]
		, objective_evaluations=Union{Missing, Int}[]
		, maximum_evaluations=Union{Missing, Int}[]
		, parameter_bound_hit=Bool[]
	)
end

function parameter_bound_hit(θ, lower, upper; tolerance=1e-6)
	return any(min.(abs.(θ .- lower), abs.(upper .- θ)) .<= tolerance)
end

function push_result_status!(
	  status
	, result;
	  model
	, step
	, sa_type
	, sa_value
	, continuation_index=missing
	, reported_target=false
	, selected=false
	, lower
	, upper
)
	push!(status, (
		  model=model
		, step=step
		, sa_type=sa_type
		, sa_value=sa_value
		, continuation_index=continuation_index
		, reported_target=reported_target
		, selected=selected
		, optimizer=result.optimizer
		, raw_status=result.raw_status
		, termination_status=result.termination_status
		, objective_value=result.f̂
		, objective_evaluations=result.objective_evaluations
		, maximum_evaluations=result.maximum_evaluations
		, parameter_bound_hit=parameter_bound_hit(result.θ̂, lower, upper)
	))
end

function check_solved(result, model, path_value)
	result.termination_status == "LOCALLY_SOLVED" ||
		error("$model sensitivity optimization failed at $sa_type = $path_value: $(result.raw_status)")
end

# Continue each model independently from its original-sample estimate.
function run_homo_continuation(path_values)
	status = empty_sa_status()
	result = (
		  θ̂=vcat(αγ_homo, αr_homo, αη_homo)
		, f̂=f_homo
		, optimizer="none"
		, raw_status="BASELINE_REUSED"
		, termination_status="LOCALLY_SOLVED"
		, objective_evaluations=missing
		, maximum_evaluations=missing
	)

	if isempty(path_values)
		push_result_status!(status, result;
			model="homo_gamma", step="baseline_reused", sa_type=sa_type,
			sa_value=sa_value, reported_target=true, selected=true,
			lower=lower_homo, upper=upper_homo)
		return result, status
	end

	for (path_index, path_value) in enumerate(path_values)
		apply_scenario_value!(sa_type, path_value)
		println("Running homogeneous continuation: index = $path_index, value = $path_value")
		result = estimate_gmm!(IV_homo, result.θ̂[1], result.θ̂[2:3], result.θ̂[4];
			b_dc_s=b_dc_s, w_b_dc_s=w_b_dc_s, lower=lower_homo,
			upper=upper_homo, W=W_opt_homo)
		check_solved(result, "homo_gamma", path_value)
		is_target = path_index == length(path_values)
		push_result_status!(status, result;
			model="homo_gamma", step="continuation", sa_type=sa_type,
			sa_value=path_value, continuation_index=path_index,
			reported_target=is_target, selected=is_target,
			lower=lower_homo, upper=upper_homo)
	end
	return result, status
end

function run_hetero_continuation(path_values, result_homo)
	status = empty_sa_status()
	result_continuation = (
		  θ̂=vcat(αγ_hetero, αr_hetero, αη_hetero)
		, f̂=f_hetero
		, optimizer="none"
		, raw_status="BASELINE_REUSED"
		, termination_status="LOCALLY_SOLVED"
		, objective_evaluations=missing
		, maximum_evaluations=missing
	)

	if isempty(path_values)
		push_result_status!(status, result_continuation;
			model="hetero", step="baseline_reused", sa_type=sa_type,
			sa_value=sa_value, reported_target=true, selected=true,
			lower=lower_hetero, upper=upper_hetero)
		return result_continuation, status
	end

	for (path_index, path_value) in enumerate(path_values)
		apply_scenario_value!(sa_type, path_value)
		println("Running heterogeneous continuation: index = $path_index, value = $path_value")
		result_continuation = estimate_gmm!(
			IV_hetero, result_continuation.θ̂[1:2],
			result_continuation.θ̂[3:4], result_continuation.θ̂[5];
			b_dc_s=b_dc_s, w_b_dc_s=w_b_dc_s, lower=lower_hetero,
			upper=upper_hetero, W=W_opt_hetero)
		check_solved(result_continuation, "hetero", path_value)
		is_target = path_index == length(path_values)
		push_result_status!(status, result_continuation;
			model="hetero", step="continuation", sa_type=sa_type,
			sa_value=path_value, continuation_index=path_index,
			reported_target=is_target, lower=lower_hetero, upper=upper_hetero)
	end

	# At the reported target, compare the continuation solution with the usual
	# nested-model start and retain the lower GMM objective.
	apply_scenario_value!(sa_type, sa_value)
	result_nested = estimate_gmm!(
		IV_hetero, [result_homo.θ̂[1], 0.0], result_homo.θ̂[2:3], result_homo.θ̂[4];
		b_dc_s=b_dc_s, w_b_dc_s=w_b_dc_s, lower=lower_hetero,
		upper=upper_hetero, W=W_opt_hetero)
	push_result_status!(status, result_nested;
		model="hetero", step="target_nested", sa_type=sa_type,
		sa_value=sa_value, reported_target=true,
		lower=lower_hetero, upper=upper_hetero)

	continuation_row = findlast(status.step .== "continuation")
	nested_row = nrow(status)
	nested_solved = result_nested.termination_status == "LOCALLY_SOLVED"
	if nested_solved && result_nested.f̂ < result_continuation.f̂
		status[nested_row, :selected] = true
		return result_nested, status
	end
	status[continuation_row, :selected] = true
	return result_continuation, status
end

result_homo, status_homo_path = run_homo_continuation(path_values)
result_hetero, status_hetero_path = run_hetero_continuation(path_values, result_homo)
apply_scenario_value!(sa_type, sa_value)

# output directories
mkpath(joinpath(pathtodata, "intermediate", "params1_sa"))
mkpath(joinpath(pathtodata, "intermediate", "params1_sa_status"))
# initialize result file if exist
rm(joinpath(pathtodata, "intermediate", "params1_sa", string("sa", index_sa, ".csv")), force=true)
rm(joinpath(pathtodata, "intermediate", "params1_sa_status", string("sa", index_sa, ".csv")), force=true)

# save results
variable = vcat(
	  ["alpha_gamma0", "alpha_gammaK"]
	, [string("alpha_r", k) for k in 0:1]
	, ["alpha_eta"]
	, ["objective_value"]
	, filter(row -> occursin(r"var_r_a", row.variable), spot_result)[:, :variable]
)
spot_variance_output = Any[spot_variance...]
output = DataFrame(
	  sa_type = fill(sa_type, length(variable))
	, sa_value = fill(sa_value, length(variable))
	, variable = variable
	, homo = vcat(Any[
		  result_homo.θ̂[1]
		, missing
		, result_homo.θ̂[2]
		, result_homo.θ̂[3]
		, result_homo.θ̂[4]
		, result_homo.f̂
	], spot_variance_output)
	, hetero = vcat(Any[
		  result_hetero.θ̂[1]
		, result_hetero.θ̂[2]
		, result_hetero.θ̂[3]
		, result_hetero.θ̂[4]
		, result_hetero.θ̂[5]
		, result_hetero.f̂
	], spot_variance_output)
)
CSV.write(joinpath(pathtodata, "intermediate", "params1_sa", string("sa", index_sa, ".csv")), output)

CSV.write(
	  joinpath(pathtodata, "intermediate", "params1_sa_status", string("sa", index_sa, ".csv"))
	, vcat(status_homo_path, status_hetero_path)
)
