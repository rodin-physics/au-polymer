using Distributed
using GLMakie
using CairoMakie

CairoMakie.activate!()
@everywhere include("calculations/calc_settings.jl")

YS = repeat(ys, 1, length(ys))
XS = permutedims(YS)

res = @showprogress pmap(
        (x, y) -> spectral_bulk(slice_energy, Location(x, y), sRing),
        XS,
        YS,
)

xs = reduce(vcat, XS)
ys = reduce(vcat, YS)

X = b1[1] .* xs + b2[1] .* ys
Y = b1[2] .* xs + b2[2] .* ys
signal = reduce(vcat, res)

## Plotting
fig = Figure(resolution = (1800, 1800))
ax =
        fig[1, 1] = Axis(
                fig,
                xlabel = "x/a₀",
                ylabel = "y/a₀",
                xlabelpadding = 0,
                ylabelpadding = 0,
                xlabelsize = 16,
                ylabelsize = 16,
                xticklabelsize = 14,
                yticklabelsize = 14,
                aspect = AxisAspect(1),
                xticklabelfont = "serif-roman",
                yticklabelfont = "serif-roman",
                xlabelfont = "serif-italic",
                ylabelfont = "serif-italic",
        )

sc = CairoMakie.scatter!(
        ax,
        X,
        Y,
        color = signal,
        strokewidth = 0,
        marker = :hexagon,
        markersize = 9,
        colormap = :oslo,
)

CairoMakie.xlims!(ax, [-120, 120])
CairoMakie.ylims!(ax, [-120, 120])

# ax.attributes
# cbar = Colorbar(fig[1,2],
#         # width = 10,
#         # limits = sc.colorrange,
#         # limits = [minimum(signal), maximum(signal)],
#         colormap = :oslo,
#         ticklabelfont = "serif-roman",
#         ticklabelsize = 14
#         # label = "A",height = Relative(1/2)
#
# )
# tightlimits!(ax)
fig
CairoMakie.save("Sp_Fun13.pdf", fig)
# fig = Figure(resolution = (1200, 900))
#
# fig[1, 1] = Axis(fig, title = "My column has size Relative(2/3)")
# fig[1, 2] = Axis(fig, title = "My column has size Auto()")
# fig[1, 3] = Colorbar(fig, width = 30)


# colsize!(fig.layout, 1, Relative(2/3))
# fig
#
