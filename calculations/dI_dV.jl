using Distributed
using Plots

@everywhere include("/Users/harshitramahalingam/Documents/Senior Semester 2/Capstone and Research Seminar/BP_Quantum_Dots/calculation/calc_settings.jl")
tip_pot = 0.15

dI_dV = pmap(x -> spectral_bulk(x, Location(0, 0), s), ωs)
# dI_dV_with_tip = pmap(x -> spectral_bulk_add_potential(x, Location(0, 0), s, tip_pot), ωs)
# dI_dV2 = pmap(x -> spectral_bulk(x, Location(0, 5), s), ωs)


## Plotting dI/dV curves
#
plot(
    ωs,
    dI_dV,
    xlabel = "Energy ω (eV)",
    ylabel = "Spectral Fn",
    label = "U = $U_val eV",
    framestyle = :box,
    linecolor = :blue,
    xrange = (-0.5, 2.0),
)

# plot!(
#     ωs.-μ,
#     dI_dV_with_tip,
#     label = "U = $U_val eV; TP = $tip_pot eV",
# )

# plot(xlabel = L"\textrm{Energy} \; \omega \; (\textrm{eV})",
#      ylabel = L"\textrm{Spectral Fn} \; \; \mathcal{A}(\omega)",
#      framestyle = :box,
#      xrange = (0.0, 4.1),
#      yrange = (0.0, 1.5),
#      guidefont = ("Times new roman", 10),
#      legendfontsize = 7,
#      title = L"\textrm{36-point Ring with 2 layers: U = 2.0 eV}")


## Bound States - Particle in a Ring

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

## Plotting bound state energies

# plot!(energies_l0, seriestype = "vline", linestyle = :dash, linecolor = :black, label = "l = 0", linewidth = 1.5)
#
# plot!(energies_l1, seriestype = "vline", linestyle = :dash, linecolor = :blue, label = "l = 1", linewidth = 1.5)
#
# plot!(energies_l2, seriestype = "vline", linestyle = :dash, linecolor = :red, label = "l = 2", linewidth = 1.5)
#
# plot!(ωs, dI_dV, linecolor = :orange,linewidth = 2, label = "U = $U_val eV",)
