using Distributed
using LaTeXStrings

@everywhere include("calculations/calc_settings.jl")

dI_dV_ring_center = pmap(x -> spectral_bulk(x, Location(0, 0), sRing), ωs)
dI_dV_ring_off_center =
    pmap(x -> spectral_bulk(x, locator(R / 2, 0), sRing), ωs)

    
gr()
plot(
    ωs .- μ,
    dI_dV_ring_center,
    xlabel = L"Energy $\omega$ (eV)",
    ylabel = "Spectral Fn",
    label = "Center",
    framestyle = :box,
    linecolor = :red,
    xrange = (-0.5, 2.0),
)
plot!(
    ωs .- μ,
    dI_dV_ring_off_center,
    label = "Off-center",
    framestyle = :box,
    linecolor = :blue,
    xrange = (-0.5, 2.0),
)
