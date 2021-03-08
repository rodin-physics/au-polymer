using Distributed
using Plots
using Colors, ColorSchemes

@everywhere include("calculations/calc_settings.jl")

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
# pyplot()
# plot(
#     # xlims = (-pos_lim, pos_lim),
#     # ylims = (-pos_lim, pos_lim),
#     aspectratio = 1,
#     title = "μ=$μ",
# )
# heatmap!(xs, ys, res)
plot(fontfamily = "sans-serif")
scatter(
    X,
    Y,
    marker_z = newsignal,
    markerstrokecolor = :white,
    markerstrokewidth = 0.0,
    markersize = 2.72,
    markershape = :hexagon,
    # color = :blackbody,
    color = cgrad([myred, mywhite, myblue], [0.41, 0.60]),
    aspect_ratio = 1,
    colorbar = true,
    # title = "QD Energy = $ϵ eV",
    framestyle = :box,
    annotate = (
        350,
        15,
        Plots.text(
            L"\delta \rho \, (10^{-3} \, \textrm{electron/orbital})",
            10,
            :dark,
            rotation = 90,
        ),
    ),
    xlabel = L"\textrm{Distance}\, (a_0)",
    ylabel = L"\textrm{Distance}\, (a_0)",
    guidefontsize = 10,
    left_margin = 3Plots.mm,
    right_margin = 9Plots.mm,
    # clim = (-0.0084,0.01504),
    clim = (-10.0, 10.0),
    label = :none,
    title = "μ = $μ eV",
    top_margin = 0.1Plots.mm,
    xtickfont = font(7, "serif"),
    ytickfont = font(7, "serif"),
)
scatter!(
    X_QD,
    Y_QD,
    aspectratio = 1,
    markersize = 2.0,
    markerstrokewidth = 1.5,
    markershape = :x,
    markercolor = :black,
    label = L"\textrm{Localised State}",
    legend = :bottomright,
    legendfontsize = 7,
)
