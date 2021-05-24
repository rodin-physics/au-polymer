using Distributed
using Plots

@everywhere include("../calculations/calc_settings.jl")

## Calculation - center and off-center of potential shape
dI_dV = @showprogress pmap(x -> spectral_bulk(x, Location(0, 0), s), ωs)
dI_dV2 = @showprogress pmap(x -> spectral_bulk(x, Location(0, 3), s), ωs)

# Plotting spectral function
plotly()
plot(ωs .- μ, dI_dV,
    xlabel = "Energy ω (eV)",
    ylabel = "Spectral Fn",
    label = "Center",
    framestyle = :box,
    linecolor = :blue,
    linewidth = 2,
    xrange = (-0.5, 2.0))

plot!(ωs .- μ, dI_dV2,
    label = "Off-Center",
    linecolor = :red,
    linewidth = 2)

##  Energies of infinite circular potential well
besselzeros = [
    [2.40483, 5.52008, 8.65373, 11.7915, 14.9309, 18.0711, 21.2116, 24.3525, 27.4935, 30.6346, 33.7758, 36.9171, 40.0584, 43.1998, 46.3412],
    [3.83171, 7.01559, 10.1735, 13.3237, 16.4706, 19.6159, 22.7601, 25.9037, 29.0468, 32.1897, 35.3323, 38.4748, 41.6171, 44.7593, 47.9015],
    [5.13562, 8.41724, 11.6198, 14.796, 17.9598, 21.117, 24.2701, 27.4206, 30.5692, 33.7165, 36.8629, 40.0084, 43.1535, 46.298, 49.4422]]

function bound_state(n, l, rad)
    return (ryd)*(besselzeros[l+1][n])^2/(sqrt(mx * my) * rad^2)
end

ns = range(1,3, step = 1)

energies_l0 = map(x -> bound_state(x,0,R),ns)
energies_l1 = map(x -> bound_state(x,1,R),ns)
energies_l2 = map(x -> bound_state(x,2,R),ns)

# Plotting energies of infinite circular potential well

# plot!(energies_l0, seriestype = "vline", linestyle = :dash, linecolor = :black, label = "l = 0", linewidth = 1.5)
#
# plot!(energies_l1, seriestype = "vline", linestyle = :dash, linecolor = :orange, label = "l = 1", linewidth = 1.5)
#
# plot!(energies_l2, seriestype = "vline", linestyle = :dash, linecolor = :green, label = "l = 2", linewidth = 1.5)
