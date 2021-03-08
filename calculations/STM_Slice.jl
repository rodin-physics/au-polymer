using Distributed
# using PyPlot
using Plots

@everywhere include("/Users/harshitramahalingam/Documents/Senior Semester 2/Capstone and Research Seminar/BP_Quantum_Dots/calculation/calc_settings.jl")

YS = repeat(ys, 1, length(ys))
XS = permutedims(YS)

res = @showprogress pmap(
    (x, y) -> spectral_bulk(slice_energy, Location(x, y), s),
    XS,
    YS,
)

# res = @showprogress pmap(
#     (x, y) -> spectral_bulk_add_potential(slice_energy, Location(x, y), s, 0.25),
#     XS,
#     YS,
# )
xs = reduce(vcat, XS)
ys = reduce(vcat, YS)

X = b1[1] .* xs + b2[1] .* ys
Y = b1[2] .* xs + b2[2] .* ys
signal = reduce(vcat, res)

## Plotting
pyplot()
Plots.scatter(
    X,
    Y,
    marker_z = signal,
    markerstrokecolor = :white,
    markerstrokewidth = 0.0,
    markersize = 5.85,
    markershape = :hexagon,
    color = :thermal,
    aspect_ratio = 1,
    leg = false,
    colorbar = true,
    # title = "Barrier height = $U_val eV, ω = $exp_energy eV",
    framestyle = :box,
    annotate = (130, 130, L"\mathcal{A} (\omega)"),
    # xlabel = L"\textrm{Distance} (a_0)",
    # ylabel = L"\textrm{Distance} (a_0)",
    xlabel = "\$\\mathrm{Distance} (a_0)\$",
    ylabel = "\$\\mathrm{Distance} (a_0)\$",
    top_margin = 12Plots.mm,
    right_margin = 2Plots.mm,
    # clim = (0.03,0.055),
    xrange = (-120, 120),
    yrange = (-120, 120),
    grid = :false,
)
## Change the text and text colour accordingly
annotate!(90, 100, Plots.text("1.32 eV", 10, :white))


# savefig("first_try.pdf")
