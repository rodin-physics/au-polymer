using Distributed
using LaTeXStrings

@everywhere include("calculations/calc_settings.jl")

dI_dV_ring_center = pmap(x -> spectral_bulk(x, Location(0, 0), sRing), ωs)
dI_dV_ring_off_center =
        pmap(x -> spectral_bulk(x, locator(R / 2, 0), sRing), ωs)
#
# dI_dV_fig8_center = pmap(x -> spectral_bulk(x, Location(0, 0), sFig8), ωs)
# dI_dV_fig8_off_center =
#     pmap(x -> spectral_bulk(x, locator(0, R * sin(π / 5)), sFig8), ωs)

## Plotting dI/dV curves
fig = Figure(resolution = (1280, 800))

ax =
        fig[1, 1] = Axis(
                fig,
                ylabel = "A",
                xlabel = "ω (eV)",
                xlabelpadding = 0,
                ylabelpadding = 0,
                xlabelsize = 16,
                ylabelsize = 16,
                xticklabelsize = 14,
                yticklabelsize = 14,
                aspect = AxisAspect(1.5),
                xticklabelfont = "serif-roman",
                yticklabelfont = "serif-roman",
                xlabelfont = "serif-roman",
                ylabelfont = "serif-italic",
        )
line_center = CairoMakie.lines!(ax, ωs .- μ, dI_dV_ring_center, color = :red)
line_off_center =
        CairoMakie.lines!(ax, ωs .- μ, dI_dV_ring_off_center, color = :blue)
CairoMakie.lines!([-.36, -.36], [-1, 2])
CairoMakie.lines!([0.2, 0.2], [-1, 2])
CairoMakie.lines!([0.9, 0.9], [-1, 2])
CairoMakie.lines!([1.3, 1.3], [-1, 2])
CairoMakie.xlims!(ax, [-0.5, 2])
CairoMakie.ylims!(ax, [0, 0.12])

fig

CairoMakie.save("dI_dV_Ring.pdf", fig)
