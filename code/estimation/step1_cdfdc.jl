#=
Brazil new energy auction structural model estimation
inverse cdf of descending-clock auction bid price truncated from above at purchase agreement price

@author: Konan Hara

@date: April 1, 2026

@software: Julia version 1.12

=#

# preparation code
include("prep.jl")

# load specified parameters
params_specified = CSV.read(joinpath("output", "result", "params_specified.csv"), DataFrame)

# unpack parameters
μpDhat = filter(row -> occursin(r"mu_pD_", row.variable), params_specified)[:, :value]
σ2pDhat = filter(row -> row.variable == "sigma2_pD", params_specified)[:, :value][1]

# configuration parameters
# descending-clock auction bid price cdf lower bound
b_dc_lb = 0
# bid price cdf grid interval
grid_interval = 0.01

# descending-clock auction purchase agreement price
p_dc = dc_mat[:, 2]

# inverse cdf of descending-clock auction bid price truncated from above at purchase agreement price
inv_cdf_b_dc = Array{Any}(undef, A_dc)
# float version of binomial(n,k) function
binomial_float(n, k) = prod((n + 1 - i) / i for i in 1:min(k, n - k); init=1.0)
# set b = 0 to cdf = 0
# this approximation is fine since F_b(0) is very close to zero
println("confirm whether the following are very close to zero")
for a = 1:A_dc
    # bid price grid truncated from above at clearing price
    b_dc_grid = b_dc_lb:grid_interval:p_dc_a[a]
    # conversion from order statistics cdf(x) to original cdf(y)
    function cdf_conversion(y, x)
        cdf(Normal(μpDhat[a], sqrt(σ2pDhat)), x) -
        sum(binomial_float(N_dc_a[a], j) * (y^j) * ((1 - y)^(N_dc_a[a] - j)) for j = (NW_dc_a[a]+1):N_dc_a[a])
    end
    # cdf of bid price up to clearing price
    cdf_b_dc_grid = zeros(length(b_dc_grid))
    for (index, b) in enumerate(b_dc_grid)
        cdf_b_dc_grid[index] = find_zero(cdf_conversion, (0, 1), p=b)
    end
    # set b = 0 to cdf = 0
    println(string("a = ", a, ": ", cdf_b_dc_grid[1]))
    cdf_b_dc_grid[1] = 0.0
    # cdf of bid price truncated from above at purchase agreement price
    inv_cdf_b_dc[a] = Array{Any}(undef, sum(a_dc .== a + A_sb))
    for (i, p) in enumerate(p_dc[a_dc.==a+A_sb])
        b_dc_grid_p = b_dc_grid[b_dc_grid.<=p]
        cdf_b_dc_grid_p = cdf_b_dc_grid[b_dc_grid.<=p]
        cdf_b_dc_trunc_grid = cdf_b_dc_grid_p / cdf_b_dc_grid_p[end]
        # linear interpolation of the inverse cdf of bid price
        Interpolations.deduplicate_knots!(cdf_b_dc_trunc_grid, move_knots=true)
        inv_cdf_b_dc[a][i] = linear_interpolation(cdf_b_dc_trunc_grid, b_dc_grid_p, extrapolation_bc=Throw())
    end
end
# a = 1:5.475072824997575e-88
# a = 2:1.6055653171525148e-6
# a = 3:1.5503150785772023e-26
# a = 4:3.8013339296156355e-6
# a = 5:5.472327646465157e-31
# a = 6:6.032762107532248e-7
# a = 7:0.0
# a = 8:0.0
# -> good
inv_cdf_b_dc = reduce(vcat, inv_cdf_b_dc)
# save
mkpath(joinpath(pathtodata, "intermediate"))
serialize(joinpath(pathtodata, "intermediate", "inv_cdf_b_dc.dat"), inv_cdf_b_dc)
