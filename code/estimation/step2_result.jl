#=
Brazil new energy auction structural model estimation
cost distribution
estimation results

@author: Konan Hara

@date: September 13, 2026

@software: Julia version 1.12

=#

# preparation code
include("prep.jl")

# output directory
mkpath(joinpath("output", "result"))

main_status_path = joinpath("output", "result", "step2_status.csv")
isfile(main_status_path) || error("Missing original-sample Step 2 status file: $main_status_path")
main_status = CSV.read(main_status_path, DataFrame)
main_status.recovered_observations[1] == main_status.active_observations[1] ||
	error("At least one active original-sample winner cost was not recovered")

dist_bs_status_path = joinpath(pathtodata, "intermediate", "params2_dist_bs_status.csv")
isfile(dist_bs_status_path) || error("Missing Step 2 distribution bootstrap status file")
dist_bs_status = CSV.read(dist_bs_status_path, DataFrame)
nrow(dist_bs_status) == bs_num ||
	error("Expected $bs_num Step 2 distribution bootstrap statuses")
all(dist_bs_status.termination_status_b .== "LOCALLY_SOLVED") ||
	error("At least one Step 2 distribution bootstrap optimization failed")
all(in.(dist_bs_status.termination_status_K,
	Ref(["LOCALLY_SOLVED", "ALMOST_LOCALLY_SOLVED"]))) ||
	error("At least one prior-capacity Tobit bootstrap optimization failed")
any(dist_bs_status.sigma2_b_bound_hit) &&
	error("At least one Step 2 bid-price variance estimate hit its bound")
any(dist_bs_status.sigma2_K_bound_hit) &&
	error("At least one prior-capacity Tobit variance estimate hit its bound")

# winners' capacity-weighted average
# load parameter estimates
gamma_c_data = CSV.read(joinpath(pathtodata, "intermediate", "gamma_c_data.csv"), DataFrame)
# use auctions with number of participants data
dropmissing!(gamma_c_data)
# variable names
variable = ["bid_price", "gamma", "bid_share_opt", "rp", "markup", "cost"]
# variable columns
col = Symbol.(variable)
# estimates
estimate = map(col) do c
	weights = gamma_c_data.capacity .* gamma_c_data.bootstrap_weight
	sum(gamma_c_data[!, c] .* weights) / sum(weights)
end
# collect bootstrap estimates
gamma_c_winner_bs = Array{Float64}(undef, bs_num, length(variable))
for index_bs in 1:bs_num
	status_path = joinpath(pathtodata, "intermediate", "params2_bs_status",
		string("bs", index_bs, ".csv"))
	isfile(status_path) || error("Missing Step 2 bootstrap status file: $status_path")
	status_b = CSV.read(status_path, DataFrame)
	status_b.recovered_observations[1] == status_b.active_observations[1] ||
		error("At least one active winner cost was not recovered in draw $index_bs")

	# load bootstrap estimates
	result_path_b = joinpath(pathtodata, "intermediate", "params2_bs", "gamma_c_data",
		string("bs", index_bs, ".csv"))
	isfile(result_path_b) || error("Missing Step 2 bootstrap result file: $result_path_b")
	gamma_c_data_index_bs = CSV.read(result_path_b, DataFrame)
	# winners' capacity-weighted average
	dropmissing!(gamma_c_data_index_bs)
	weights = gamma_c_data_index_bs.capacity .* gamma_c_data_index_bs.bootstrap_weight
	sum(weights) > 0 || error("Step 2 bootstrap draw $index_bs has zero winner weight")
	gamma_c_winner_bs[index_bs, :] = map(col) do c
		sum(gamma_c_data_index_bs[!, c] .* weights) / sum(weights)
	end
end
# save results
# bootstrap standard errors
vcm = cov(gamma_c_winner_bs)
se = sqrt.([vcm[i,i] for i in 1:size(vcm)[1]])
# central 95% bootstrap percentile intervals
ci_lower = [quantile(gamma_c_winner_bs[:, i], 0.025)
	for i in 1:size(gamma_c_winner_bs, 2)]
ci_upper = [quantile(gamma_c_winner_bs[:, i], 0.975)
	for i in 1:size(gamma_c_winner_bs, 2)]
# results table
gamma_c_winner_result = (
	  variable = variable
	, estimate = estimate
	, se = se
	, ci_lower = ci_lower
	, ci_upper = ci_upper
)
CSV.write(joinpath("output", "result", "gamma_c_winner_result.csv")
	, gamma_c_winner_result)

# cost distribution linear projection across auctions
# load parameter estimates
c_dist_lp = CSV.read(joinpath("output", "result", "c_dist_lp.csv"), DataFrame)
# variable names
variable = c_dist_lp[:, :variable]
# collect bootstrap estimates
c_dist_lp_bs = Array{Float64}(undef, bs_num, length(variable))
c_dist_lp_undiscount_bs = Array{Float64}(undef, bs_num, length(variable))
for index_bs in 1:bs_num
	# load bootstrap estimates
	c_dist_lp_index_bs = CSV.read(joinpath(pathtodata, "intermediate", "params2_bs", "c_dist_lp"
			, string("bs", index_bs, ".csv"))
		, DataFrame)
	# main
	c_dist_lp_bs[index_bs,:] = c_dist_lp_index_bs[:,:main]
	# undiscounted
	c_dist_lp_undiscount_bs[index_bs,:] = c_dist_lp_index_bs[:,:undiscounted]
end
# save results
# bootstrap standard errors
vcm = cov(c_dist_lp_bs)
main_se = sqrt.([vcm[i,i] for i in 1:size(vcm)[1]])
vcm = cov(c_dist_lp_undiscount_bs)
undiscounted_se = sqrt.([vcm[i,i] for i in 1:size(vcm)[1]])
# results table
c_dist_lp_result = (
	  variable = variable
	, main_estimate = c_dist_lp[:, :main]
	, main_se = main_se
	, undiscounted_estimate = c_dist_lp[:, :undiscounted]
	, undiscounted_se = undiscounted_se
)
CSV.write(joinpath("output", "result", "c_dist_lp_result.csv")
	, c_dist_lp_result)

# Joint (log gamma,cost) projections and aggregate summaries under all three
# weighting conventions. Draws with a rank-deficient projection are omitted
# parameter by parameter and their valid count is recorded.
function aggregate_bootstrap_table(main_path, bootstrap_dir, output_path)
	main = CSV.read(main_path, DataFrame)
	value_columns = names(main)[2:end]
	result = DataFrame(
		  variable=String[]
		, weighting=String[]
		, cost_type=String[]
		, estimate=Float64[]
		, se=Union{Missing,Float64}[]
		, valid_draws=Int[]
	)

	bootstrap = Dict(
		name => fill(NaN, bs_num, nrow(main)) for name in value_columns)
	for index_bs in 1:bs_num
		path = joinpath(bootstrap_dir, string("bs", index_bs, ".csv"))
		isfile(path) || error("Missing Step 2 bootstrap aggregate file: $path")
		draw = CSV.read(path, DataFrame)
		draw.variable == main.variable || error("Aggregate variables do not align in draw $index_bs")
		for name in value_columns
			bootstrap[name][index_bs, :] = Float64.(draw[!, name])
		end
	end

	for name in value_columns
		undiscount = endswith(name, "_undiscount")
		weighting = replace(name, "_undiscount" => "")
		for j in 1:nrow(main)
			values = bootstrap[name][:, j]
			finite_values = values[isfinite.(values)]
			se = length(finite_values) > 1 ? std(finite_values) : missing
			push!(result, (
				  variable=String(main.variable[j])
				, weighting=weighting
				, cost_type=undiscount ? "undiscounted" : "main"
				, estimate=Float64(main[j, name])
				, se=se
				, valid_draws=length(finite_values)
			))
		end
	end
	CSV.write(output_path, result)
end

aggregate_bootstrap_table(
	joinpath("output", "result", "gamma_c_lp.csv"),
	joinpath(pathtodata, "intermediate", "params2_bs", "gamma_c_lp"),
	joinpath("output", "result", "gamma_c_lp_result.csv"),
)
aggregate_bootstrap_table(
	joinpath("output", "result", "gamma_c_summary.csv"),
	joinpath(pathtodata, "intermediate", "params2_bs", "gamma_c_summary"),
	joinpath("output", "result", "gamma_c_summary_result.csv"),
)

# (γ,c) distribution summary statistics
# total capacity-weighted correlation
function weighted_mean_complete(x, w)
	keep = .!ismissing.(x) .& isfinite.(coalesce.(x, NaN))
	any(keep) || return NaN
	x_complete = Float64.(x[keep])
	w_complete = Float64.(w[keep])
	return sum(x_complete .* w_complete) / sum(w_complete)
end

function weighted_cor(x, y, w)
	w_sum = sum(w)
	x_mean = sum(x .* w) / w_sum
	y_mean = sum(y .* w) / w_sum
	x_dev = x .- x_mean
	y_dev = y .- y_mean
	sum(w .* x_dev .* y_dev) / sqrt(sum(w .* x_dev .^ 2) * sum(w .* y_dev .^ 2))
end

# total capacity-weighted average across auctions
# load parameter estimates
gamma_c_dist = CSV.read(joinpath("output", "result", "gamma_c_dist.csv"), DataFrame)
# use auctions with total capacity and participant-count data
auction_dist_weights = auction_sb_a[.!ismissing.(auction_sb_a.number_participant),
	[:auction_id, :total_capacity]]
gamma_c_dist = leftjoin(gamma_c_dist, auction_dist_weights; on=:auction_id)
any(ismissing, gamma_c_dist.total_capacity) &&
	error("Missing auction capacity weight in gamma_c_dist.csv")
totcap_sb_a = Float64.(gamma_c_dist.total_capacity)
# variable names
dist_variable = names(gamma_c_dist)[2:(end - 1)]
# variable columns
col = Symbol.(dist_variable)
# estimates
estimate = map(col) do c
	weighted_mean_complete(gamma_c_dist[!, c], totcap_sb_a)
end
variable = vcat(dist_variable, "gamma_c_cor_across")
estimate = vcat(estimate
	, weighted_cor(gamma_c_dist[:, :gamma_mean], gamma_c_dist[:, :c_mean], totcap_sb_a))
# collect bootstrap estimates
gamma_c_dist_bs = Array{Float64}(undef, bs_num, length(variable))
for index_bs in 1:bs_num
	# load bootstrap estimates
	gamma_c_dist_index_bs = CSV.read(joinpath(pathtodata, "intermediate", "params2_bs", "gamma_c_dist"
			, string("bs", index_bs, ".csv"))
		, DataFrame)
	gamma_c_dist_index_bs = leftjoin(gamma_c_dist_index_bs, auction_dist_weights; on=:auction_id)
	any(ismissing, gamma_c_dist_index_bs.total_capacity) &&
		error("Missing auction capacity weight in gamma_c_dist bootstrap draw $index_bs")
	totcap_bs = Float64.(gamma_c_dist_index_bs.total_capacity)
	# total capacity-weighted average across auctions
	gamma_c_dist_bs[index_bs, 1:length(col)] = map(col) do c
		weighted_mean_complete(gamma_c_dist_index_bs[!, c], totcap_bs)
	end
	gamma_c_dist_bs[index_bs, end] = weighted_cor(
		gamma_c_dist_index_bs[:, :gamma_mean]
		, gamma_c_dist_index_bs[:, :c_mean]
		, totcap_bs)
end
# save results
# bootstrap standard errors
vcm = cov(gamma_c_dist_bs)
se = sqrt.([vcm[i,i] for i in 1:size(vcm)[1]])
# results table
gamma_c_dist_result = (
	  variable = variable
	, estimate = estimate
	, se = se
)
CSV.write(joinpath("output", "result", "gamma_c_dist_result.csv")
	, gamma_c_dist_result)
