using Distributed
using Plots

## Change parameters in calc_settings.jl
@everywhere include("calculation/calc_settings.jl")

YS = repeat(ys, 1, length(ys))
XS = permutedims(YS)

tip_pot = 0.3                       # Value of STM tip potential

## Animation loop
anim = @animate for ii in range(0.00, 2.0, step = 0.005)
    res = @showprogress pmap((x, y) -> spectral_bulk(ii, Location(x, y), s), XS, YS)

    xs = reduce(vcat, XS)
    ys = reduce(vcat, YS)

    X = b1[1] .* xs + b2[1] .* ys
    Y = b1[2] .* xs + b2[2] .* ys
    signal = reduce(vcat, res)

    shift_val = round(ii - μ, digits = 2)

    scatter(
        X,
        Y,
        marker_z = signal,
        markerstrokecolor = :white,
        markerstrokewidth = 0.001,
        markersize = 3.5,
        color = :thermal,
        aspect_ratio = 1,
        leg = false,
        colorbar = true,
        title = "Barrier height = $U_val eV, ω = $shift_val eV",
        framestyle = :box,
        annotate = (285, 182, L"\mathcal{A} (\omega)"),
        xlabel = L"\textrm{Distance}(a_0)",
        ylabel = L"\textrm{Distance}(a_0)",
        right_margin = 4Plots.mm,
        # clim = (0.01,0.2)
    )
    savefig("animation/anim_frame-$(ii/.005).png")
end

# saving gif
gif(anim, "spectral_fn_anim.gif", fps = 30)
