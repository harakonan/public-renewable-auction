#=
Brazil new energy auction structural model estimation
cost distribution

@author: Konan Hara

@date: September 13, 2026

@software: Julia version 1.12

=#

# preparation code
include("prep.jl")

# Simple bootstrap control
# SLURM_ARRAY_TASK_ID maps directly to `b`. If empty/missing, b = 0 (main estimation).
b = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "0"))

# Remove this draw's prior outputs before loading or calculating anything. A
# failed rerun must not leave stale files that appear to be current results.
if b > 0
	bootstrap_output_paths = [
		joinpath(pathtodata, "intermediate", "params2_bs", subdir, "bs$(b).csv")
		for subdir in (
			"gamma_c_data", "c_dist_lp", "gamma_c_lp", "gamma_c_summary", "gamma_c_dist")
	]
	push!(bootstrap_output_paths,
		joinpath(pathtodata, "intermediate", "params2_bs_status", "bs$(b).csv"))
	for output_path in bootstrap_output_paths
		rm(output_path, force=true)
	end
end

# use auctions with N_sb in step2, while keeping a_sb_list itself unchanged
# so bootstrap draws stay aligned with step1 and step2_dist.
a_sb_list_step2 = a_sb_list[.!ismissing.(N_sb_a)]

# load specified parameters
params_specified = CSV.read(joinpath("output", "result", "params_specified.csv"), DataFrame)

# load first stage / step2_dist estimates
if b == 0
    params1_hat = CSV.read(joinpath("output", "result", "params1_hat.csv"), DataFrame)

    params2_dist_result = CSV.read(joinpath("output", "result", "params2_dist_result.csv"), DataFrame)
    params2_dist_hat = DataFrame(
          variable = params2_dist_result[:, :variable]
        , estimate = Float64.(params2_dist_result[:, :estimate])
    )
else
    params1_hat = CSV.read(joinpath(pathtodata, "intermediate", "params1_bs", string("bs", b, ".csv")), DataFrame)

    params2_dist_bs = CSV.read(joinpath(pathtodata, "intermediate", "params2_dist_bs.csv"), DataFrame)
    params2_dist_hat = DataFrame(
          variable = names(params2_dist_bs)
        , estimate = Float64.(collect(params2_dist_bs[b, :]))
    )
end

# unpack parameters
βDhat = filter(row -> occursin(r"beta_D._sb", row.variable), params_specified)[:, :value]
σ2Dhat = filter(row -> row.variable == "sigma2_D_sb", params_specified)[:, :value][1]

# unpack estimates
αγhat = filter(row -> occursin(r"alpha_gamma", row.variable), params1_hat)[:, :hetero]
αrhat = filter(row -> occursin(r"alpha_r", row.variable), params1_hat)[:, :hetero]
βChat = filter(row -> occursin(r"beta_C", row.variable), params2_dist_hat)[:, :estimate]
σ2Chat = filter(row -> row.variable == "sigma2_C", params2_dist_hat)[:, :estimate][1]
βKhat = filter(row -> occursin(r"beta_K", row.variable), params2_dist_hat)[:, :estimate]
σ2Khat = filter(row -> row.variable == "sigma2_K", params2_dist_hat)[:, :estimate][1]
βbhat = Array{Union{Missing,Float64}}(missing, A_sb)
βbhat[.!ismissing.(N_sb_a)] = filter(row -> occursin(r"beta_b", row.variable), params2_dist_hat)[:, :estimate]
σ2bhat = filter(row -> row.variable == "sigma2_b", params2_dist_hat)[:, :estimate][1]

# draws for winning probability function calculation
# set seed for auction draws
Random.seed!(seed_winprob)
# draws are fixed across auctions
# draw standard normals for competitor bid prices
# maximum number of competitors is maximum number of participants - 1
Umax = rand(Normal(), (2, SW, maximum(N_sb_a[.!ismissing.(N_sb_a)]) - 1))
# draw standard normals for procurement quantity
UD = rand(Normal(), SW)

# draws for (γ,c) distribution simulation
# set seed for bidder draws
Random.seed!(seed_γc)
# draws are fixed across auctions
Uγc = rand(Normal(), (2, Sγc))

# Separate normal draws generate prior-started wind capacity from its
# shifted-lognormal Tobit without changing the existing bid/capacity draws.
Random.seed!(seed_winprob + 10_000)
UK_win = rand(Normal(), SW, maximum(N_sb_a[.!ismissing.(N_sb_a)]) - 1)
Random.seed!(seed_γc + 10_000)
UK_dist = rand(Normal(), Sγc)

# step2 functions
include("step2_func.jl")

function step2_status(sample, output, observation_weights, active_auctions, W_functions)
	valid = .!ismissing.(output.c_f)
	active_rows = in.(a_sb, Ref(Set(active_auctions)))
	fallback_auctions = [
		a for a in active_auctions if W_functions[a].use_asymptote_fallback]
	return DataFrame(
		  sample=[sample]
		, observations=[length(valid)]
		, recovered_observations=[count(valid)]
		, recovered_weight=[sum(observation_weights[valid])]
		, total_weight=[sum(observation_weights)]
		, missing_gamma=[count(ismissing, output.γ_f)]
		, missing_risk_premium=[count(ismissing, output.rp_f)]
		, missing_markup=[count(ismissing, output.markup_f)]
		, missing_cost=[count(ismissing, output.c_f)]
		, active_auctions=[length(active_auctions)]
		, active_auction_ids=[join(a_sb_a[active_auctions], ";")]
		, active_observations=[count(active_rows)]
		, active_participants=[sum(Int.(N_sb_a[active_auctions]))]
		, active_total_capacity=[sum(Float64.(auction_sb_a[active_auctions, :total_capacity]))]
		, asymptote_fallback_auctions=[length(fallback_auctions)]
		, asymptote_fallback_auction_ids=[join(a_sb_a[fallback_auctions], ";")]
	)
end

function joint_output_table(output)
	table = DataFrame(variable=joint_projection_parameters)
	for (j, name) in enumerate(distribution_weightings)
		table[!, Symbol(name)] = output.joint_lp[:, j]
		table[!, Symbol(name * "_undiscount")] = output.joint_lp_undiscount[:, j]
	end
	return table
end

function joint_summary_table(output)
	table = DataFrame(variable=joint_summary_parameters)
	for (j, name) in enumerate(distribution_weightings)
		table[!, Symbol(name)] = output.joint_summary[:, j]
		table[!, Symbol(name * "_undiscount")] = output.joint_summary_undiscount[:, j]
	end
	return table
end

if b == 0
    println("--- Running Main Estimation ---")

    sb_mat_step2 = sb_mat
    a_sb_step2 = a_sb
	prior_assurance_100_step2 = prior_assurance_100_sb
	observation_weights = ones(obs_sb)
	active_auctions = a_sb_list_step2

    # generate winning probability functions
    W_functions = Array{Any}(undef, A_sb)
    for a in active_auctions
        aparam = func_aparam(
			a, βChat, σ2Chat, βKhat, σ2Khat, αγhat, αrhat, βbhat, σ2bhat)
        W_functions[a] = func_W_functions(aparam, Umax, UD, UK_win)
    end

    # save winning probability functions
    winprob_grid = DataFrame(
          auction_id=reduce(vcat, [repeat([a], length(W_functions[a].W_b_grid)) for a in active_auctions])
        , bid_price=reduce(vcat, [collect(W_functions[a].W_b_grid) for a in active_auctions])
        , winprob=reduce(vcat, [W_functions[a].W_grid for a in active_auctions])
		, winprob_inf=reduce(vcat, [fill(W_functions[a].W_inf,
			length(W_functions[a].W_b_grid)) for a in active_auctions])
		, use_asymptote_fallback=reduce(vcat, [fill(
			W_functions[a].use_asymptote_fallback,
			length(W_functions[a].W_b_grid)) for a in active_auctions])
    )
    CSV.write(joinpath(pathtodata, "intermediate", "winprob_grid.csv"), winprob_grid)

	# (γ,c) distribution estimation
    output = func_step2(
		sb_mat_step2, a_sb_step2, active_auctions, βChat, σ2Chat, βKhat,
		σ2Khat, αγhat, αrhat, βbhat, σ2bhat, winprob_grid, Uγc,
		prior_assurance_100_step2, UK_dist)

    # concise original-sample diagnostics
	println("Observed bid-share mean: ", mean(sb_mat_step2[.!ismissing.(N_sb), 1]))
	println("Model bid-share mean: ", mean(skipmissing(output.q_ast_f)))
	println("Winner gamma mean: ", mean(skipmissing(output.γ_f)))
	println("Winner risk-premium mean: ", mean(skipmissing(output.rp_f)))
	println("Winner markup mean: ", mean(skipmissing(output.markup_f)))
	println("Winner cost mean: ", mean(skipmissing(output.c_f)))

    # save parameter estimates
    # winner costs
    gamma_c_data = (
          auction_id = a_sb
        , bid_share = sb_mat_step2[:, 1]
        , bid_price = sb_mat_step2[:, 2]
		, prior_started_wind_capacity_100 = prior_assurance_100_step2
        , gamma = output.γ_f
        , bid_share_opt = output.q_ast_f
        , rp = output.rp_f
        , markup = output.markup_f
        , cost = output.c_f
        , capacity = sb_mat_step2[:, 3]
		, bootstrap_weight = observation_weights
    )
    CSV.write(joinpath(pathtodata, "intermediate", "gamma_c_data.csv"), gamma_c_data)

    # cost distribution linear projection across auctions
    c_dist_lp = (
        variable=["intercept", "t", "t2", "t3", "l", "t*l", "t2*l", "t3*l", "sigma2"]
        , main = output.c_lp
        , undiscounted = output.c_lp_undiscount
    )
    CSV.write(joinpath("output", "result", "c_dist_lp.csv"), c_dist_lp)
	CSV.write(joinpath("output", "result", "gamma_c_lp.csv"), joint_output_table(output))
	CSV.write(joinpath("output", "result", "gamma_c_summary.csv"), joint_summary_table(output))

	# Joint participant draws used as counterfactual primitives.
	cf_primitives = DataFrame(
		  auction_id=Int[]
		, draw=Int[]
		, gamma=Float64[]
		, cost=Float64[]
		, capacity=Float64[]
	)
	for a in active_auctions
		valid = .!ismissing.(output.γ_sim[:, a]) .&
			.!ismissing.(output.c_sim[:, a]) .&
			.!ismissing.(output.Capacity_sim[:, a])
		draw_index = findall(valid)
		append!(cf_primitives, DataFrame(
			  auction_id=fill(a_sb_a[a], length(draw_index))
			, draw=draw_index
			, gamma=Float64.(output.γ_sim[valid, a])
			, cost=Float64.(output.c_sim[valid, a])
			, capacity=Float64.(output.Capacity_sim[valid, a])
		))
	end
	counterfactual_directory = joinpath(
		pathtodata, "intermediate", "counterfactual_task")
	mkpath(counterfactual_directory)
	CSV.write(joinpath(counterfactual_directory, "participant_primitives.csv")
		, cf_primitives)

    # (γ,c) distribution summary statistics for each auction
    gamma_c_dist = (
	          auction_id = active_auctions
        , loggamma_mean = output.logγ_mean
        , loggamma_var = output.logγ_var
        , gamma_mean = output.γ_mean
        , gamma_var = output.γ_var
        , c_mean = output.c_mean
        , c_var = output.c_var
        , loggamma_c_cor = output.logγ_c_cor
        , gamma_c_cor = output.γ_c_cor
    )
    CSV.write(joinpath("output", "result", "gamma_c_dist.csv"), gamma_c_dist)
	CSV.write(joinpath("output", "result", "step2_status.csv"),
		step2_status("main", output, observation_weights, active_auctions, W_functions))
else
    println("--- Running Bootstrap Iteration $b ---")

	# Reproduce the Step 1 global firm-pairs draw. Keep the observed auctions and
	# rows fixed and carry multiplicities as frequency weights.
    seed_bs = seed_bs_list[b]
    Random.seed!(seed_bs)
	weights_all, bootstrap_diagnostics = draw_global_firm_weights()
	observation_weights = weights_all[1:obs_sb]
	sb_mat_step2 = sb_mat
	a_sb_step2 = a_sb
	prior_assurance_100_step2 = prior_assurance_100_sb
	active_auctions = [a for a in a_sb_list_step2 if
		sum(observation_weights[a_sb .== a]) > 0]

    # generate winning probability functions
    W_functions = Array{Any}(undef, A_sb)
    for a in active_auctions
        aparam = func_aparam(
			a, βChat, σ2Chat, βKhat, σ2Khat, αγhat, αrhat, βbhat, σ2bhat)
        W_functions[a] = func_W_functions(aparam, Umax, UD, UK_win)
    end

    winprob_grid = DataFrame(
          auction_id=reduce(vcat, [repeat([a], length(W_functions[a].W_b_grid)) for a in active_auctions])
        , bid_price=reduce(vcat, [collect(W_functions[a].W_b_grid) for a in active_auctions])
        , winprob=reduce(vcat, [W_functions[a].W_grid for a in active_auctions])
		, winprob_inf=reduce(vcat, [fill(W_functions[a].W_inf,
			length(W_functions[a].W_b_grid)) for a in active_auctions])
		, use_asymptote_fallback=reduce(vcat, [fill(
			W_functions[a].use_asymptote_fallback,
			length(W_functions[a].W_b_grid)) for a in active_auctions])
    )

    # (γ,c) distribution estimation
    output = func_step2(
		sb_mat_step2, a_sb_step2, active_auctions, βChat, σ2Chat, βKhat,
		σ2Khat, αγhat, αrhat, βbhat, σ2bhat, winprob_grid, Uγc,
		prior_assurance_100_step2, UK_dist)

    # save bootstrap outputs
    mkpath(joinpath(pathtodata, "intermediate", "params2_bs", "gamma_c_data"))
    mkpath(joinpath(pathtodata, "intermediate", "params2_bs", "c_dist_lp"))
    mkpath(joinpath(pathtodata, "intermediate", "params2_bs", "gamma_c_dist"))
	mkpath(joinpath(pathtodata, "intermediate", "params2_bs", "gamma_c_lp"))
	mkpath(joinpath(pathtodata, "intermediate", "params2_bs", "gamma_c_summary"))
	mkpath(joinpath(pathtodata, "intermediate", "params2_bs_status"))

    gamma_c_data = (
          auction_id = a_sb_step2
        , bid_share = sb_mat_step2[:, 1]
        , bid_price = sb_mat_step2[:, 2]
		, prior_started_wind_capacity_100 = prior_assurance_100_step2
        , gamma = output.γ_f
        , bid_share_opt = output.q_ast_f
        , rp = output.rp_f
        , markup = output.markup_f
        , cost = output.c_f
        , capacity = sb_mat_step2[:, 3]
		, bootstrap_weight = observation_weights
    )
    CSV.write(joinpath(pathtodata, "intermediate", "params2_bs", "gamma_c_data", string("bs", b, ".csv")), gamma_c_data)

    # cost distribution linear projection across auctions
    c_dist_lp = (
        variable=["intercept", "t", "t2", "t3", "l", "t*l", "t2*l", "t3*l", "sigma2"]
        , main = output.c_lp
        , undiscounted = output.c_lp_undiscount
    )
    CSV.write(joinpath(pathtodata, "intermediate", "params2_bs", "c_dist_lp", string("bs", b, ".csv")), c_dist_lp)
	CSV.write(joinpath(pathtodata, "intermediate", "params2_bs", "gamma_c_lp", string("bs", b, ".csv")), joint_output_table(output))
	CSV.write(joinpath(pathtodata, "intermediate", "params2_bs", "gamma_c_summary", string("bs", b, ".csv")), joint_summary_table(output))

    gamma_c_dist = (
	          auction_id = active_auctions
        , loggamma_mean = output.logγ_mean
        , loggamma_var = output.logγ_var
        , gamma_mean = output.γ_mean
        , gamma_var = output.γ_var
        , c_mean = output.c_mean
        , c_var = output.c_var
        , loggamma_c_cor = output.logγ_c_cor
        , gamma_c_cor = output.γ_c_cor
    )
    CSV.write(joinpath(pathtodata, "intermediate", "params2_bs", "gamma_c_dist", string("bs", b, ".csv")), gamma_c_dist)

	status = step2_status(
		"bootstrap", output, observation_weights, active_auctions, W_functions)
	status[!, :bootstrap_draw] = [b]
	status[!, :bootstrap_effective_sample_size] =
		[bootstrap_diagnostics.bootstrap_effective_sample_size]
	status[!, :bootstrap_active_firms] = [bootstrap_diagnostics.bootstrap_active_firms]
	CSV.write(joinpath(pathtodata, "intermediate", "params2_bs_status", string("bs", b, ".csv")), status)
end
