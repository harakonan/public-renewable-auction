#=
Brazil new energy auction structural model estimation
cost distribution
sensitivity analysis

@author: Konan Hara

@date: September 13, 2026

@software: Julia version 1.12

=#

# preparation code
include("prep.jl")

# Sensitivity analysis index
# SLURM_ARRAY_TASK_ID maps directly to the sensitivity-analysis index.
index_sa = parse(Int, get(ENV, "SLURM_ARRAY_TASK_ID", "1"))

# use auctions with N_sb in step2, while keeping a_sb_list itself unchanged
# so outputs match the bootstrap step2 sample.
a_sb_list_step2 = a_sb_list[.!ismissing.(N_sb_a)]
id_sb_step2 = collect(id_sb)[.!ismissing.(N_sb)]
sb_mat_step2 = sb_mat[id_sb_step2, :]
a_sb_step2 = a_sb[id_sb_step2]

# load specified parameters
params_specified = CSV.read(joinpath("output", "result", "params_specified.csv"), DataFrame)

# load step2_dist estimates
params2_dist_result = CSV.read(joinpath("output", "result", "params2_dist_result.csv"), DataFrame)
params2_dist_hat = DataFrame(
      variable = params2_dist_result[:, :variable]
    , estimate = Float64.(params2_dist_result[:, :estimate])
)

# load first-stage sensitivity estimates
params1_hat = CSV.read(joinpath(pathtodata, "intermediate", "params1_sa", string("sa", index_sa, ".csv")), DataFrame)

# sensitivity-analysis scenario recorded by step1_sa.jl
sa_type = String(params1_hat[1, :sa_type])
sa_value = Float64(params1_hat[1, :sa_value])
if !(sa_type in ["sigma2_xi", "mu_p", "sigma2_p", "delta"])
    error("Unknown sensitivity-analysis type: $sa_type")
end

# sensitivity-analysis spot variance recorded by step1_sa.jl
spot_variance = Float64.(filter(row -> occursin(r"var_r_a", row.variable), params1_hat)[:, :hetero])
σ2r_sb_a = spot_variance[a_sb_a]
σ2r_dc_a = spot_variance[a_dc_a]

# sensitivity-analysis discount factor
if sa_type == "delta"
    δ = sa_value
end
δtilde_sb_a = [sum([δ^t for t = l_sb_a[a]:1.0:(l_sb_a[a] + T - 1)]) / T for a in 1:A_sb]
δtilde_dc_a = [sum([δ^t for t = l_dc_a[a]:1.0:(l_dc_a[a] + T - 1)]) / T for a in 1:A_dc]

println("--- Running Sensitivity Analysis $index_sa ---")
println("Sensitivity type: $sa_type")
println("Sensitivity value: $sa_value")

# unpack parameters
βDhat = filter(row -> occursin(r"beta_D._sb", row.variable), params_specified)[:, :value]
σ2Dhat = filter(row -> row.variable == "sigma2_D_sb", params_specified)[:, :value][1]

# unpack estimates
αγhat = Float64.(filter(row -> occursin(r"alpha_gamma", row.variable), params1_hat)[:, :hetero])
αrhat = Float64.(filter(row -> occursin(r"alpha_r", row.variable), params1_hat)[:, :hetero])
βChat = Float64.(filter(row -> occursin(r"beta_C", row.variable), params2_dist_hat)[:, :estimate])
σ2Chat = filter(row -> row.variable == "sigma2_C", params2_dist_hat)[:, :estimate][1]
βKhat = Float64.(filter(row -> occursin(r"beta_K", row.variable), params2_dist_hat)[:, :estimate])
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

Random.seed!(seed_winprob + 10_000)
UK_win = rand(Normal(), SW, maximum(N_sb_a[.!ismissing.(N_sb_a)]) - 1)
Random.seed!(seed_γc + 10_000)
UK_dist = rand(Normal(), Sγc)

# step2 functions
include("step2_func.jl")

# generate winning probability functions
W_functions = Array{Any}(undef, A_sb)
for a in a_sb_list_step2
    aparam = func_aparam(
		a, βChat, σ2Chat, βKhat, σ2Khat, αγhat, αrhat, βbhat, σ2bhat)
    W_functions[a] = func_W_functions(aparam, Umax, UD, UK_win)
end

winprob_grid = DataFrame(
      auction_id=reduce(vcat, [repeat([a], length(W_functions[a].W_b_grid)) for a in a_sb_list_step2])
    , bid_price=reduce(vcat, [collect(W_functions[a].W_b_grid) for a in a_sb_list_step2])
    , winprob=reduce(vcat, [W_functions[a].W_grid for a in a_sb_list_step2])
	, winprob_inf=reduce(vcat, [fill(W_functions[a].W_inf,
		length(W_functions[a].W_b_grid)) for a in a_sb_list_step2])
	, use_asymptote_fallback=reduce(vcat, [fill(
		W_functions[a].use_asymptote_fallback,
		length(W_functions[a].W_b_grid)) for a in a_sb_list_step2])
)

# (γ,c) distribution estimation
output = func_step2(
	sb_mat_step2, a_sb_step2, a_sb_list_step2, βChat, σ2Chat, βKhat,
	σ2Khat, αγhat, αrhat, βbhat, σ2bhat, winprob_grid, Uγc,
	prior_assurance_100_sb[id_sb_step2], UK_dist)

# save sensitivity-analysis outputs
mkpath(joinpath(pathtodata, "intermediate", "params2_sa", "gamma_c_data"))

gamma_c_data = (
      auction_id = a_sb_step2
    , bid_share = sb_mat_step2[:, 1]
    , bid_price = sb_mat_step2[:, 2]
    , prior_started_wind_capacity_100 = prior_assurance_100_sb[id_sb_step2]
    , gamma = output.γ_f
    , bid_share_opt = output.q_ast_f
    , rp = output.rp_f
    , markup = output.markup_f
    , cost = output.c_f
    , capacity = sb_mat_step2[:, 3]
)
CSV.write(joinpath(pathtodata, "intermediate", "params2_sa", "gamma_c_data", string("sa", index_sa, ".csv")), gamma_c_data)

fallback_auctions = [
	a for a in a_sb_list_step2 if W_functions[a].use_asymptote_fallback]
status = DataFrame(
      sensitivity_index=[index_sa]
    , sa_type=[sa_type]
    , sa_value=[sa_value]
    , observations=[length(output.c_f)]
    , recovered_observations=[count(x -> !ismissing(x), output.c_f)]
    , missing_cost=[count(ismissing, output.c_f)]
	, asymptote_fallback_auctions=[length(fallback_auctions)]
	, asymptote_fallback_auction_ids=[join(a_sb_a[fallback_auctions], ";")]
)
mkpath(joinpath(pathtodata, "intermediate", "params2_sa_status"))
CSV.write(joinpath(pathtodata, "intermediate", "params2_sa_status", string("sa", index_sa, ".csv")), status)
