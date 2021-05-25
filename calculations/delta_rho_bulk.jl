using Distributed
using Plots
using Colors, ColorSchemes

@everywhere include("../calculations/calc_settings.jl")

## Calculation
YS = repeat(ys, 1, length(ys))
XS = permutedims(YS)

res = @showprogress pmap((x, y) -> δρR(Location(x, y), s), XS, YS)

xs = reduce(vcat, XS)
ys = reduce(vcat, YS)

X = b1[1] .* xs + b2[1] .* ys
Y = b1[2] .* xs + b2[2] .* ys

signal = reduce(vcat, res)
newsignal = signal .* 1000

## Plotting
myred = Colors.RGB(255 / 255, 55 / 255, 0 / 255)
myblue = Colors.RGB(105 / 255, 141 / 255, 247 / 255)
mywhite = Colors.RGB(1, 1, 1)

drho_scheme = ColorScheme([
    get(ColorScheme([myred, mywhite, myblue]), i) for i = 0.0:0.01:1.0
])

loadcolorscheme(
    :drho_scheme,
    [get(drho_scheme, i) for i = 0.0:0.01:1.0],
    "custom scheme",
)

# Plotting δρ in charge density
plotly()
scatter(
    X,
    Y,
    marker_z = newsignal,
    markerstrokewidth = 0.0,
    markersize = 5.3,
    markershape = :hexagon,
    color = cgrad([myred, mywhite, myblue], [0.41, 0.60]),
    aspect_ratio = 1,
    colorbar = true,
    colorbar_title = "δρ (10⁻³ electron/orbital)",
    framestyle = :box,
    xlabel = "Distance (a₀)",
    ylabel = "Distance (a₀)",
    clim = (-10.0, 10.0),
    label = :none,
    title = "μ = $μ eV",
    xrange = (-100, 100),
    yrange = (-100, 100),
)

# Plot unit cells with nonzero local potential
scatter!(
    X2,
    Y2,
    markersize = 3.5,
    markerstrokewidth = 0.001,
    markershape = :xcross,
    markercolor = :black,
    label = :none,
)
