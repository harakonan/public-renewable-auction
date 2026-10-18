#=
Brazil new energy auction structural model estimation
prespecified parameters

@author: Konan Hara

@date: April 1, 2026

@software: Julia version 1.12

=#

# preparation code
include("prep.jl")

# MLE for f(D) (use OLS to get the same result)
# sealed-bid auction
cons = repeat([1],A_sb)
# remove missing N_sb_a auctions
X_D_sb = hcat(cons, (t_sb_a + l_sb_a), m_sb_a, N_sb_a)[.!ismissing.(N_sb_a), :]
Y_D_sb = log.(D_sb_a[.!ismissing.(N_sb_a)])
βDhat_sb = (X_D_sb'*X_D_sb)\X_D_sb'*Y_D_sb
e = Y_D_sb-X_D_sb*βDhat_sb
σ2Dhat_sb = e'e/size(X_D_sb, 1)
# descending-clock auction
cons = repeat([1],A_dc)
X_D_dc = hcat(cons, (t_dc_a + l_dc_a), m_dc_a, N_dc_a)
Y_D_dc = log.(D_dc_a)
βDhat_dc = (X_D_dc'*X_D_dc)\X_D_dc'*Y_D_dc
e = Y_D_dc-X_D_dc*βDhat_dc
σ2Dhat_dc = e'e/A_dc

# MLE for f(p|D) (use OLS to get the same result)
cons = repeat([1],A_dc)
X_pD = hcat(cons, log.(D_dc_a), (t_dc_a + l_dc_a), m_dc_a, N_dc_a)
Y_pD = p_dc_a
βpDhat = (X_pD'*X_pD)\X_pD'*Y_pD
μpDhat = X_pD*βpDhat
e = Y_pD-X_pD*βpDhat
σ2pDhat = e'e/A_dc

# predicted μp
X_p = hcat(cons, X_D_dc * βDhat_dc, (t_dc_a + l_dc_a), m_dc_a, N_dc_a)
μphat = X_p*βpDhat
# predicted σ2p
σ2phat = (βpDhat[2]^2)*σ2Dhat_dc + σ2pDhat

# save results
# specified parameters
variable = vcat(
	# procurement capacity distribution
	  [string("beta_D", k, "_sb") for k in 0:3], ["sigma2_D_sb"]
	, [string("beta_D", k, "_dc") for k in 0:3], ["sigma2_D_dc"]
	# equilibrium clearing price distribution
	, [string("beta_pD", k) for k in 0:4], ["sigma2_pD"]
	# expected clearing price conditional on observed procurement capacity
	# for each descending-clock auction
	, [string("mu_pD_a", a) for a in (A_sb+1):(A_sb+A_dc)]
	# expected clearing price for each descending-clock auction
	, [string("mu_p_a", a) for a in (A_sb+1):(A_sb+A_dc)]
	# clearing price variance
	, ["sigma2_p"]
)
# values
v = vcat(
	  βDhat_sb
	, σ2Dhat_sb
	, βDhat_dc
	, σ2Dhat_dc
	, βpDhat
	, σ2pDhat
	, μpDhat
	, μphat
	, σ2phat
)
params_specified = (
	  variable = variable  
	, value = v
)
CSV.write(joinpath("output", "result", "params_specified.csv")
	, params_specified)
