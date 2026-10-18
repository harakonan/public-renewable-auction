# Calibrate the primary two-risk-type common-support model.

include("counterfactual_model.jl")

using Dates

function main()
	cache_path = get(ENV, "CF_BASELINE_FILE", joinpath(
		intermediate_directory, "counterfactual_baseline.jls"))
	loggamma_sd = parse(Float64, get(ENV, "CF_PRIMARY_LOGGAMMA_SD",
		"0.0383370135"))
	rho = representative_primitive.factual_rho
	draws, achieved_rho, outside = two_type_draws(loggamma_sd, rho, 0.0)
	map = allocation_map()
	anchor_result = marginal_participant_bid_anchor()
	upper_bid_anchor = anchor_result.anchor
	anchor_path = joinpath(intermediate_directory,
		"counterfactual_marginal_participant_anchor.csv")
	CSV.write(anchor_path, anchor_result.diagnostic)
	println("Calibrating the primary two-type common-support model.")
	println("Baseline marginal-participant bid anchor: ",
		round(upper_bid_anchor, digits=8))
	println("Marginal-participant anchor diagnostics written to $anchor_path")
	flush(stdout)
	calibration = calibrate_registered_pool(draws, map, upper_bid_anchor)
	payload = (
		cache_kind="two_type_common_support_marginal_participant_bid_anchor",
		generated_at=string(now()),
		model="two discrete gamma types; common upper bid support",
		entry_rule="two type-specific cost cutoffs",
		upper_bid_anchor_method="expected maximum participant bid",
		upper_bid_anchor_draws=parse(Int, get(
			ENV, "CF_MARGINAL_BID_DRAWS", "5000")),
		factual_rho=rho,
		achieved_rho=achieved_rho,
		primary_loggamma_sd=loggamma_sd,
		gamma_mean=mean(draws.gamma),
		K_anchor=calibration.boundary.K,
		upper_bid_anchor=upper_bid_anchor,
		kappa=calibration.kappa,
		upper_bid=calibration.boundary.bid,
		registered_tilt_alpha=calibration.alpha,
		registered_moment_residual=calibration.moment_residual,
		registered_moment_rmse=calibration.moment_rmse,
		registered_moment_maximum=calibration.moment_maximum,
		registered_entry_rate=calibration.entry_rate,
		registered_entry_rate_residual=calibration.entry_residual,
		primitive_floor_share=mean(outside.active),
		n_aux=controls.n_aux,
		n_winprob=controls.n_winprob,
		n_alloc=controls.n_alloc,
		seed=controls.seed)
	mkpath(dirname(cache_path))
	temporary_path = tempname(dirname(cache_path))
	open(temporary_path, "w") do stream
		serialize(stream, payload)
	end
	mv(temporary_path, cache_path; force=true)
	println("Baseline cache written to $cache_path")
	println("Baseline entry rate: $(round(calibration.entry_rate, digits=6))")
	println("Baseline kappa per MWh: ",
		@sprintf("%.8e", calibration.kappa))
	println("Registered moment RMSE: $(round(calibration.moment_rmse, digits=6))")
end

main()
