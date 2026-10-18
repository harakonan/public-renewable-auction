#=
Brazil new energy auction structural model estimation
first-stage GMM goodness of fit

@author: Konan Hara

@date: September 13, 2026

@software: Julia version 1.12
=#

ENV["STEP1_DEFINITIONS_ONLY"] = "true"
include("../../code/estimation/step1.jl")
params1_hat = CSV.read(joinpath("output", "result", "params1_hat.csv"), DataFrame)

inv_cdf_b_dc = deserialize(joinpath(pathtodata, "intermediate", "inv_cdf_b_dc.dat"))

U_b_dc_s, w_b_dc_s = uniform_quadrature_legendre(Sgmm)
b_dc_s = [
	inv_cdf_b_dc[i](U_b_dc_s[s])
	for i in 1:obs_dc, s in 1:Sgmm
]

function parameter_value(model::String, variable::String)
	row = filter(r -> r.variable == variable, params1_hat)
	return Float64(row[1, Symbol(model)])
end

function model_parameters(model::String)
	if model == "hetero"
		αγ = [
			  parameter_value(model, "alpha_gamma0")
			, parameter_value(model, "alpha_gammaK")
		]
	else
		throw(ArgumentError("Unsupported GMM goodness-of-fit model: $model"))
	end
	αr = [
		  parameter_value(model, "alpha_r0")
		, parameter_value(model, "alpha_r1")
	]
	αη = parameter_value(model, "alpha_eta")
	return αγ, αr, αη
end

function censored_normal_cdf(x, q_tilde, ση, lower, upper)
	if x < lower
		return 0.0
	elseif x >= upper
		return 1.0
	end
	return cdf(Normal(), (x - q_tilde) / ση)
end

function censored_normal_mean(q_tilde, ση, lower, upper)
	standard_normal = Normal()
	z_lower = (lower - q_tilde) / ση
	z_upper = (upper - q_tilde) / ση
	p_lower = cdf(standard_normal, z_lower)
	p_upper = ccdf(standard_normal, z_upper)
	p_interior = max(1.0 - p_lower - p_upper, 0.0)
	interior_first_moment =
		q_tilde * p_interior +
		ση * (pdf(standard_normal, z_lower) - pdf(standard_normal, z_upper))
	return lower * p_lower + interior_first_moment + upper * p_upper
end

function latent_bid_shares(αγ, αr)
	K_log_sb = sb_mat[:, 3]
	K_log_dc = dc_mat[:, 3]
	γ_sb = exp.(αγ[1] .+ αγ[2] .* K_log_sb)
	γ_dc = exp.(αγ[1] .+ αγ[2] .* K_log_dc)
	if length(αr) == 1
		μr_sb = fill(αr[1], obs_sb)
		μr_dc = fill(αr[1], obs_dc)
	else
		μr_sb = αr[1] .+ αr[2] .* m_sb
		μr_dc = αr[1] .+ αr[2] .* m_dc
	end

	q_tilde_sb = [
		foc_sb[a_sb[i]](sb_mat[i, 2], γ_sb[i], μr_sb[i])
		for i in 1:obs_sb
	]
	q_tilde_dc = Matrix{Float64}(undef, obs_dc, Sgmm)
	for i in 1:obs_dc, s in 1:Sgmm
		q_tilde_dc[i, s] = solve(
			  ZeroProblem(
				q -> foc_dc[a_dc[i] - A_sb](
					q, b_dc_s[i, s], γ_dc[i], μr_dc[i]
				)
				, x0
			  )
			, Order0()
		)
	end
	return q_tilde_sb, q_tilde_dc
end

function model_diagnostics(model::String; grid=collect(0.0:0.01:1.0))
	αγ, αr, αη = model_parameters(model)
	ση = exp(αη)
	q_tilde_sb, q_tilde_dc = latent_bid_shares(αγ, αr)

	q = [sb_mat[:, 1]; dc_mat[:, 1]]
	q_lower = [sb_mat[:, 4]; dc_mat[:, 4]]
	q_upper = [sb_mat[:, 5]; dc_mat[:, 5]]
	p_lower = Vector{Float64}(undef, obs)
	p_upper = Vector{Float64}(undef, obs)
	predicted_mean = Vector{Float64}(undef, obs)
	predicted_cdf = Matrix{Float64}(undef, obs, length(grid))

	for i in 1:obs_sb
		p_lower[i] = cdf(Normal(), (q_lower[i] - q_tilde_sb[i]) / ση)
		p_upper[i] = ccdf(Normal(), (q_upper[i] - q_tilde_sb[i]) / ση)
		predicted_mean[i] = censored_normal_mean(
			q_tilde_sb[i], ση, q_lower[i], q_upper[i]
		)
		for (j, x) in enumerate(grid)
			predicted_cdf[i, j] = censored_normal_cdf(
				x, q_tilde_sb[i], ση, q_lower[i], q_upper[i]
			)
		end
	end

	for i in 1:obs_dc
		row = obs_sb + i
		p_lower[row] = 0.0
		p_upper[row] = 0.0
		predicted_mean[row] = 0.0
		predicted_cdf[row, :] .= 0.0
		for s in 1:Sgmm
			q_tilde = q_tilde_dc[i, s]
			p_lower[row] +=
				w_b_dc_s[s] * cdf(Normal(), (q_lower[row] - q_tilde) / ση)
			p_upper[row] +=
				w_b_dc_s[s] * ccdf(Normal(), (q_upper[row] - q_tilde) / ση)
			predicted_mean[row] += w_b_dc_s[s] * censored_normal_mean(
				q_tilde, ση, q_lower[row], q_upper[row]
			)
			for (j, x) in enumerate(grid)
				predicted_cdf[row, j] += w_b_dc_s[s] * censored_normal_cdf(
					x, q_tilde, ση, q_lower[row], q_upper[row]
				)
			end
		end
	end

	observed_lower = q .< q_lower .+ small_num
	observed_upper = q .> q_upper .- small_num
	observed_interior = .!(observed_lower .| observed_upper)
	diagnostics = DataFrame(
		  model=fill(model, obs)
		, auction_id=[a_sb_id; a_dc_id]
		, auction_type=[
			String.(auction_data_sb.auction_type)
			; String.(auction_data_dc.auction_type)
		  ]
		, auction_format=[
			fill("sealed_bid", obs_sb)
			; fill("descending_clock", obs_dc)
		  ]
		, contract_type=[
			String.(auction_data_sb.contract_type)
			; String.(auction_data_dc.contract_type)
		  ]
		, bid_share=q
		, bid_share_lower=q_lower
		, bid_share_upper=q_upper
		, observed_lower=observed_lower
		, observed_interior=observed_interior
		, observed_upper=observed_upper
		, predicted_lower=p_lower
		, predicted_interior=max.(1.0 .- p_lower .- p_upper, 0.0)
		, predicted_upper=p_upper
		, predicted_mean=predicted_mean
	)

	group_masks = [
		  ("all", "all", trues(obs))
		, ("auction_format", "sealed_bid", [trues(obs_sb); falses(obs_dc)])
		, (
			"auction_format"
			, "descending_clock"
			, [falses(obs_sb); trues(obs_dc)]
		  )
	]
	for auction_type in sort(unique(diagnostics.auction_type))
		push!(
			group_masks
			, ("auction_type", auction_type, diagnostics.auction_type .== auction_type)
		)
	end
	predictive_summary = DataFrame(
		  model=String[]
		, group_type=String[]
		, group_value=String[]
		, bid_share=Float64[]
		, observed_cdf=Float64[]
		, predicted_cdf=Float64[]
	)
	for (group_type, group_value, mask) in group_masks
		for (j, x) in enumerate(grid)
			push!(predictive_summary, (
				  model=model
				, group_type=group_type
				, group_value=group_value
				, bid_share=x
				, observed_cdf=mean(q[mask] .<= x)
				, predicted_cdf=mean(predicted_cdf[mask, j])
			))
		end
	end
	return diagnostics, predictive_summary
end

models = ["hetero"]
diagnostics = DataFrame[]
predictive_cdfs = DataFrame[]
for model in models
	diagnostic, predictive_cdf = model_diagnostics(model)
	push!(diagnostics, diagnostic)
	push!(predictive_cdfs, predictive_cdf)
end

CSV.write(
	  joinpath(
		pathtodata
		, "intermediate"
		, "bid_share_gmm_diagnostics.csv"
	  )
	, vcat(diagnostics...)
)
CSV.write(
	  joinpath(
		"notes"
		, "gmm"
		, "step1_gmm_predictive_cdf.csv"
	  )
	, vcat(predictive_cdfs...)
)
