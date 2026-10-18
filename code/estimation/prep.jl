#=
Brazil new energy auction structural model estimation
preparation code

@author: Konan Hara

@date: September 13, 2026

@software: Julia version 1.12

=#

# package loading
using DataFrames
using CSV
using Serialization
using Interpolations
using Statistics
using Distributions
using Random
using Roots
using LinearAlgebra
using FastGaussQuadrature
using Tables

# set working directory to project folder
cd(dirname(dirname((@__DIR__))))

# path to data
pathtodata = readline(open("data.txt", "r"))

# load auction data
auction_data = CSV.read(joinpath(pathtodata, "processed", "generation_auction_wind20.csv")
	, DataFrame)
# remove observations with the possibility of partial bid quantity (partial_bid_quantity)
auction_data = auction_data[auction_data.partial_bid_quantity.==0, :]

# load spot market variance
spot_result = CSV.read(joinpath("output", "result", "spot_result.csv")
	, DataFrame)
# spot market variance for each auction
spot_variance = filter(row -> occursin(r"var_r", row.variable), spot_result)[:,:est]

# configuration parameters
# auction id
a_sb_a = unique(auction_data[auction_data.descending_clock .== 0,:auction_id])
a_dc_a = unique(auction_data[auction_data.descending_clock .== 1,:auction_id])
# discount factor
δ = 0.95
# purchase agreement length
T = 20
# initial value for descending-clock FOC root search
x0 = 0.5
# parameter bounds for GMM
# homogeneous risk aversion
# (αγ, αr0, αr1, αη)
lower_homo = [log(0.05), -20.0, -1.5, log(0.002)]
upper_homo = [log(3.00), 40.0, 1.5, log(0.30)]
# heterogeneous risk aversion
# (αγ0, αγ1, αr0, αr1, αη)
lower_hetero = [-10.0, -5.0, -20.0, -1.5, log(0.002)]
upper_hetero = [10.0, 5.0, 40.0, 1.5, log(0.30)]
# number of points for GMM numerical integration
Sgmm = 10
# buffer to make sure that log(0) and 1/0 do not happen in optimization
small_num = 1e-6
# number of bootstrap
bs_num = 200
# set seed starting point for bootstrap draws
seed_bs_start = 1
# set seeds for bootstrap draws
Random.seed!(seed_bs_start)
seed_bs_list = rand(UInt32, bs_num)
# smoothing parameter for nonparametric winning probability estimation
# smoothing parameter for bid price is about 1/30 of the level of a typical bid (60)
# Reference: Ryan (2022)
hb = 2.0
# smoothing parameter for procurement capacity is about 1/30 of the level of a typical procurement quantity (450)
hD = 15.0
# number of auction draws to calculate winning probability function
SW = 2000
# set seed for auction draws to calculate winning probability function
seed_winprob = 1
# winning probability function generator parameters
# winning probability bounds
W_lb = 0.01
W_ub = 1.0 - W_lb
# tolerance for convergence to the positive lower asymptote
W_asymptote_tol = 0.001
# bid price initial search grid
b_range = 20:2:80
# grid interval
grid_interval = 0.2
# number of draws to simulate (γ,c) distribution
Sγc = 10000
# set seed for (γ,c) distribution simulation
seed_γc = 1

# unpack parameters
σ2r_sb_a = spot_variance[a_sb_a]
σ2r_dc_a = spot_variance[a_dc_a]

# descending-clock auction indicator
I_dc = auction_data[:,:descending_clock]

# auction level variables

# sealed-bid auction
auction_data_sb = auction_data[I_dc .== 0,:]
var_symbol_sb_a = [
	  :auction_id
	, :year_decimal
	, :lead_time
	, :number_participant
	, :total_capacity
	, :procurement_quantity
    , :bid_share_lb_auction
    , :quantity_minimum_unit
    , :spot_price
	, :mc_operation
	, :contract_type
	]
auction_sb_a = unique(auction_data_sb[:, var_symbol_sb_a])
t_sb_a = auction_sb_a[:,:year_decimal]
l_sb_a = auction_sb_a[:,:lead_time]
N_sb_a = auction_sb_a[:,:number_participant]
D_sb_a = auction_sb_a[:,:procurement_quantity]
qlb_sb_a = auction_sb_a[:, :bid_share_lb_auction]
Qmu_sb_a = auction_sb_a[:, :quantity_minimum_unit]
m_sb_a = auction_sb_a[:,:mc_operation]
r_sb_a = auction_sb_a[:,:spot_price]
qtd_sb_a = auction_sb_a[:, :contract_type] .== "QTD"
A_sb = size(auction_sb_a)[1]
δtilde_sb_a = [sum([δ^t for t=l_sb_a[a]:1.0:(l_sb_a[a]+T-1)])/T for a=1:A_sb]

# descending-clock auction
auction_data_dc = auction_data[I_dc .== 1,:]
var_symbol_dc_a = [
	  :auction_id
	, :year_decimal
	, :lead_time
	, :number_participant
	, :number_winner
	, :procurement_quantity
	, :spot_price
	, :mc_operation
	, :contract_type
	, :clearing_price
	]
auction_dc_a = unique(auction_data_dc[:, var_symbol_dc_a])
t_dc_a = auction_dc_a[:,:year_decimal]
l_dc_a = auction_dc_a[:,:lead_time]
N_dc_a = auction_dc_a[:,:number_participant]
NW_dc_a = auction_dc_a[:,:number_winner]
D_dc_a = auction_dc_a[:,:procurement_quantity]
r_dc_a = auction_dc_a[:, :spot_price]
m_dc_a = auction_dc_a[:,:mc_operation]
qtd_dc_a = auction_dc_a[:, :contract_type] .== "QTD"
p_dc_a = auction_dc_a[:,:clearing_price]
A_dc = size(auction_dc_a)[1]
δtilde_dc_a = [sum([δ^t for t=l_dc_a[a]:1.0:(l_dc_a[a]+T-1)])/T for a=1:A_dc]

# bidder level variables

# sealed-bid auction
# sample size
obs_sb = size(auction_data_sb)[1]
# auction variables
a_sb_id = auction_data_sb[:, :auction_id]
a_sb = Int.(indexin(a_sb_id, a_sb_a))
t_sb = auction_data_sb[:,:year_decimal]
l_sb = auction_data_sb[:,:lead_time]
N_sb = auction_data_sb[:,:number_participant]
Qmu_sb = auction_data_sb[:, :quantity_minimum_unit]
r_sb = auction_data_sb[:,:spot_price]
m_sb = auction_data_sb[:,:mc_operation]
qtd_sb = auction_data_sb[:, :contract_type] .== "QTD"
# bidder variables
# Auction-recorded physical assurance defines the portfolio share, quantity,
# and feasible bounds.
function build_auction_bidder_matrix(auction_data_local)
	return hcat(
		  Float64.(auction_data_local.bid_share)
		, Float64.(auction_data_local.bid_price)
		, Float64.(auction_data_local.capacity)
		, Float64.(auction_data_local.bid_share_lb)
		, Float64.(auction_data_local.bid_share_ub)
	)
end
sb_mat = build_auction_bidder_matrix(auction_data_sb)
prior_assurance_100_sb = Float64.(auction_data_sb.prior_started_wind_capacity) ./ 100.0
first_positive_prior_time = minimum(t_sb[prior_assurance_100_sb .> 0])
structural_zero_prior_sb_a = t_sb_a .< first_positive_prior_time

# descending-clock auction
# sample size
obs_dc = size(auction_data_dc)[1]
# auction variables
a_dc_id = auction_data_dc[:, :auction_id]
a_dc = A_sb .+ Int.(indexin(a_dc_id, a_dc_a))
t_dc = auction_data_dc[:,:year_decimal]
l_dc = auction_data_dc[:,:lead_time]
Qmu_dc = auction_data_dc[:, :quantity_minimum_unit]
r_dc = auction_data_dc[:,:spot_price]
m_dc = auction_data_dc[:,:mc_operation]
qtd_dc = auction_data_dc[:, :contract_type] .== "QTD"
# bidder variables
dc_mat = build_auction_bidder_matrix(auction_data_dc)
prior_assurance_100_dc = Float64.(auction_data_dc.prior_started_wind_capacity) ./ 100.0

# other variables used in bootstrap
a_sb_list = unique(a_sb)
id_sb = 1:obs_sb
a_dc_list = unique(a_dc)
id_dc = 1:obs_dc
obs_a = zeros(Int,A_sb+A_dc)
for a in a_sb_list
	obs_a[a] = sum(a_sb .== a)
end
for a in a_dc_list
	obs_a[a] = sum(a_dc .== a)
end

# create sealed-bid/descending-clock combined vectors for convenience
a_vec = [a_sb; a_dc]
m_vec = [m_sb; m_dc]
qtd_vec = [qtd_sb; qtd_dc]
obs = obs_sb+obs_dc

# Draw firms globally so one bootstrap multiplicity applies to every retained
# observation of the same firm across auctions and auction formats.
function draw_cluster_multiplicities(ids)
	levels = sort(unique(ids))
	counts = Dict(level => 0 for level in levels)
	for level in rand(levels, length(levels))
		counts[level] += 1
	end
	return [counts[id] for id in ids], count(>(0), values(counts))
end

function draw_global_firm_weights()
	firm_ids = [auction_data_sb.parent; auction_data_dc.parent]
	multiplicities, active_firms = draw_cluster_multiplicities(firm_ids)
	weights = Float64.(multiplicities)
	weights_sb = weights[1:obs_sb]
	weights_dc = weights[(obs_sb + 1):end]
	weight_sum = sum(weights)
	diagnostics = (
		  bootstrap_weight_sum=weight_sum
		, bootstrap_nonzero_observations=count(>(0), weights)
		, bootstrap_effective_sample_size=weight_sum^2 / sum(abs2, weights)
		, bootstrap_active_auctions_sb=length(unique(a_sb_id[weights_sb .> 0]))
		, bootstrap_active_auctions_dc=length(unique(a_dc_id[weights_dc .> 0]))
		, bootstrap_active_firms=active_firms
		, bootstrap_total_firms=length(unique(firm_ids))
	)
	return weights, diagnostics
end
