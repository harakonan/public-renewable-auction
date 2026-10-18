#=
Brazil new energy auction structural model estimation
wholesale market variance evaluation using historical spot prices

@author: Konan Hara

@date: April 1, 2026

@software: Julia version 1.12

=#

# package loading
using DataFrames
using CSV
using Distributions
using FiniteDifferences

# set working directory to project folder
cd(dirname(dirname((@__DIR__))))

# path to data
pathtodata = readlines(open("data.txt", "r"))[1]

# load auction data
auction_data = CSV.read(joinpath(pathtodata, "processed", "generation_auction_wind20.csv")
	, DataFrame)

# load annual spot price data
spot_price_data = CSV.read(joinpath(pathtodata, "processed", "spot_price_annual.csv")
	, DataFrame)

# configuration parameters
# discount factor
δ = 0.95
# purchase agreement length
T = 20

# auction level variables
var_symbol = [
	  :auction_id
	, :lead_time
	]
auction_a = unique(auction_data[:, var_symbol])
l_a = auction_a[:, :lead_time]
A = length(auction_a[:, :auction_id])

# estimate spot market variance from spot price historical data (2001--2022)
spot_price = spot_price_data[:, :spot_price]
obs_spot_price = length(spot_price)
# spot price AR(1) model
# # random walk
# X = spot_price[1:end-1]
# mean reversion model
cons = repeat([1],obs_spot_price-1)
X = hcat(cons,spot_price[1:end-1])
Y = spot_price[2:end]
coefhat = (X'*X)\X'*Y
ρhat = coefhat[end]
e = Y-X*coefhat
σ2ξhat = e'e/(obs_spot_price-1)

# standard error calculation for spot price AR(1) model parameter estimates
# log-likelihood function for spot price AR(1) model
function ll_spotAR1(β)
	coef = β[1:2]
	σ2ξ = β[3]
	ll = [-((1/2)*log(σ2ξ) + (1/2)*(1/σ2ξ)*(spot_price[2:end][t] - coef[1] - coef[2]*spot_price[1:end-1][t])^2) for t=1:(obs_spot_price-1)]
	return(ll)
end
# standard errors for spot price AR(1) model parameter estimates
# variance covariance matrix for MLE with outer product approximation method
function mle_vcm_approx(ll, obs, estimate; perturb_scale = 1e-6)
    ∂ll_∂param_approx = zeros(obs,length(estimate))
    for element in 1:length(estimate)
        perturb_vec = zeros(length(estimate))
        perturb_vec[element] = estimate[element]*perturb_scale
        estimate_perturb = estimate + perturb_vec
        # gradient for each element
        ∂ll_∂param_approx[:,element] = (ll(estimate) - ll(estimate_perturb)) ./ perturb_vec[element]
    end    

    # information matrix for estimate
    info_estimate_approx = zeros(length(estimate),length(estimate))
    for i in 1:obs
        info_estimate_approx = info_estimate_approx + (∂ll_∂param_approx[i,:] * transpose(∂ll_∂param_approx[i,:]))
    end

    # variance-covariance matrix for estimate
    vcm_estimate_approx = (info_estimate_approx)^(-1)
    # return variance-covariance matrix
    return vcm_estimate_approx
end
vcm_spotAR1 = mle_vcm_approx(ll_spotAR1, obs_spot_price-1, vcat(coefhat,σ2ξhat))
spotAR1_se = sqrt.([vcm_spotAR1[i,i] for i in 1:size(vcm_spotAR1)[1]])
# # does not change much by perturbation scale
# vcm_spotAR1 = mle_vcm_approx(ll_spotAR1, obs_spot_price-1, vcat(coefhat,σ2ξhat), perturb_scale = 1e-3)
# sqrt.([vcm_spotAR1[i,i] for i in 1:size(vcm_spotAR1)[1]])

# function that calculates discounted spot market variance given
# ρ, AR(1) coefficient
# σ2ξ, AR(1) residual variance
# δ, discount factor
# T, purchase agreement length
# l, lead time
function func_spot_variance(ρ,σ2ξ,δ,T,l)
	l_fc = [floor(l),ceil(l)]
	l_plus_T = [Int(l_fc[k]+T) for k=1:2]
	var_prep = zeros(2)
	for k in 1:2
		for t in 1:(l_plus_T[k]-1)
			if t <= l_fc[k]
				var_prep[k] += (((δ^(l_fc[k]))*(ρ^(l_fc[k]-t))*((1-(δ^T)*(ρ^T))/(1-δ*ρ)))^2)
			else
				var_prep[k] += (((δ^t)*(1-(δ^(l_plus_T[k]-t))*(ρ^(l_plus_T[k]-t)))/(1-δ*ρ))^2)
			end
		end
	end
	spot_variance = l == l_fc[1] ? (σ2ξ/(T^2)) * var_prep[1] :
		(σ2ξ/(T^2)) * ((l_fc[2]-l)*var_prep[1] + (l-l_fc[1])*var_prep[2])
	return spot_variance
end

# variance of the average discounted spot price
σ2rhat = [func_spot_variance(ρhat,σ2ξhat,δ,T,l_a[a]) for a=1:A]
# sensitivity checks
# # random walk
# σ2rhat = [func_spot_variance(ρhat,σ2ξhat,δ,T,l_a[a]) for a=1:A]

# calculate standard error of spot market variance
spot_variance_se = Array{Float64}(undef, A)
for a in 1:A
	spot_variance_Tlfix(param) = func_spot_variance(param[1],param[2],δ,T,l_a[a])
	spot_variance_grad = grad(central_fdm(5,1), spot_variance_Tlfix, vcat(ρhat,σ2ξhat))[1]
	spot_variance_se[a] = sqrt(spot_variance_grad'vcm_spotAR1[2:3,2:3]spot_variance_grad)
end

# save results
# output directory
mkpath(joinpath("output", "result"))
# variable names
variable = vcat(
	# AR1 parameter values
	  ["A","rho","var_xi"]
	# spot market variance for each auction
	, [string("var_r_a", a) for a in 1:A]
)
# estimates
est = vcat(
	  coefhat
	, σ2ξhat
	, σ2rhat
)
# standard errors
se = vcat(
	  spotAR1_se
	, spot_variance_se
)
spot_result = (
	  variable = variable  
	, est = est
	, se = se
)
CSV.write(joinpath("output", "result", "spot_result.csv")
	, spot_result)
