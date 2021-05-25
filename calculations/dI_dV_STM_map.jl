using Distributed
using Plots
using PlotlyBase
using LaTeXStrings

@everywhere include("../calculations/calc_settings.jl")

## Calculation
YS = repeat(ys, 1, length(ys))
XS = permutedims(YS)

res = @showprogress pmap(
    (x, y) -> spectral_bulk(slice_energy, Location(x, y), s),
    XS,
    YS,
)

xs = reduce(vcat, XS)
ys = reduce(vcat, YS)

X = b1[1] .* xs + b2[1] .* ys
Y = b1[2] .* xs + b2[2] .* ys
signal = reduce(vcat, res)

# Plotting

plotly()
scatter(
    X,
    Y,
    zcolor = signal,
    markerstrokewidth = 0.0,
    markersize = 5.85,
    markershape = :hexagon,
    color = :oslo,
    colorbar = true,
    aspect_ratio = 1,
    leg = false,
    # title = "ω = $exp_energy eV",
    framestyle = :box,
    xlabel = "Distance (a₀)",
    ylabel = "Distance (a₀)",
    clim = (0.02,0.088),
    xrange = (-100, 100),
    yrange = (-100, 100),
)

# Plots.savefig("test.pdf")
