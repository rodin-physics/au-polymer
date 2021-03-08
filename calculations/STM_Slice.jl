using Distributed
using Plots
using PlotlyBase

@everywhere include("calculations/calc_settings.jl")

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

## Plotting
plotly()
plot(
    xlabel = L"\textrm{Distance} (a_0)",
    ylabel = L"\textrm{Distance} (a_0)",
    grid = :false,
    aspect_ratio = 1,
    leg = false,
    extra_plot_kwargs = KW(
        :include_mathjax => "cdn",
        :yaxis => KW(:automargin => true),
        :xaxis => KW(:domain => "auto"),
    ),
)
scatter!(
    X,
    Y,
    zcolor = signal,
    markerstrokecolor = :white,
    markerstrokewidth = 0.0,
    markersize = 5.85,
    markershape = :hexagon,
    color = :thermal,
    colorbar = true,
    # # title = "Barrier height = $U_val eV, ω = $exp_energy eV",
    # framestyle = :box,
    # annotate = (130, 130, L"\mathcal{A} (\omega)"),
    # # xlabel = L"\textrm{Distance} (a_0)",
    # # ylabel = L"\textrm{Distance} (a_0)",
    # xlabel = "\$\\mathrm{Distance} (a_0)\$",
    # ylabel = "\$\\mathrm{Distance} (a_0)\$",
    # top_margin = 12Plots.mm,
    # right_margin = 2Plots.mm,
    # # clim = (0.03,0.055),
    # xrange = (-120, 120),
    # yrange = (-120, 120),
)
Plots.savefig("test.pdf")
## Change the text and text colour accordingly
annotate!(90, 100, Plots.text("1.32 eV", 10, :white))


# savefig("first_try.pdf")
# signal
# pgfplotsx()
# y = rand(100)
# plot(
#     0:10:100,
#     rand(11, 4),
#     lab = "lines",
#     w = 3,
#     palette = cgrad(:grays),
#     fill = 0,
#     α = 0.6,
# )
# scatter!(
#     y,
#     zcolor = abs.(y .- 0.5),
#     m = (:heat, 0.8, Plots.stroke(1, :green)),
#     ms = 10 * abs.(y .- 0.5) .+ 4,
#     lab = "grad",
# )
#
