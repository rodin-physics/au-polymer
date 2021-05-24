using Distributed
using Plots

@everywhere include("../calculations/calc_settings.jl")


## Calculation
location_list = range(-7,7,step = 1)
full_map = @showprogress pmap(x -> pmap(y -> spectral_bulk(x, Location(y, 0), s), location_list), ωs)
full_grid = reshape([(full_map...)...], length(location_list), :)

distx = b1[1] .* location_list
disty = b1[2] .* location_list

# Plotting

plotly()
heatmap(ωs.-μ,
        distx .* 0.0529177249,
        full_grid,
        c = :jet, clims = (0.03, 0.08),
        xlabel = "Energy(eV)", ylabel = "Position(nm)")
