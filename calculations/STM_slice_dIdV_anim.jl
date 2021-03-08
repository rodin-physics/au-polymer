using Distributed
using Plots

## Change parameters in calc_settings.jl
@everywhere include("calculation/calc_settings.jl")

YS = repeat(ys, 1, length(ys))
XS = permutedims(YS)

dI_dV = pmap(x -> spectral_bulk(x, locator(0, R * sin(π / 5)), s), ωs)
ω_step = 0.005

## Animation loop
anim = @animate for ii = 0:400
        res = @showprogress pmap(
                (x, y) -> spectral_bulk(ii * ω_step, Location(x, y), s),
                XS,
                YS,
        )

        xs = reduce(vcat, XS)
        ys = reduce(vcat, YS)

        X = b1[1] .* xs + b2[1] .* ys
        Y = b1[2] .* xs + b2[2] .* ys
        signal = reduce(vcat, res)

        shift_val = round(ii * ω_step - μ, digits = 2)

        s1 = plot(
                ωs .- μ,
                dI_dV,
                xlabel = "Energy ω (eV)",
                ylabel = "Spectral Fn",
                # label = "U = $U_val eV",
                framestyle = :box,
                linecolor = :blue,
                xrange = (-0.5, 2.0),
                leg = false,
                title = "Barrier height = $U_val eV",
        )
        scatter!(
                [ii * ω_step - μ],
                [spectral_bulk(ii * ω_step, locator(0, R * sin(π / 5)), s)],
        )


        s2 = scatter(
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
                xlims = (-250, 250),
                ylims = (-150, 150),
                framestyle = :box,
                # annotate = (285, 182, L"\mathcal{A} (\omega)"),
                # xlabel = L"\textrm{Distance}(a_0)",
                # ylabel = L"\textrm{Distance}(a_0)",
                right_margin = 4Plots.mm,
                # clim = (0.01,0.2)
        )
        l = @layout([a{0.5h}; b])
        plot(s1, s2, layout = l)

        savefig("animation/fig8_Middle_of_Ring/anim_frame_$ii.png")
end

# saving gif
gif(anim, "spectral_fn_anim.gif", fps = 30)
