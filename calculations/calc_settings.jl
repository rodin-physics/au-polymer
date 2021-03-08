include("../src/density.jl")
using ProgressMeter
using Plots
using LaTeXStrings

## System Parameters
μ = 0.46
T = 0.0
slice_energy = 0.15

## --- Change parameters here ---
R = 21.8 / a0       # Radius of localized state ring and potential ring

# Local Potential
U_val = 0.5
numP = 80
Potential = mkRing(numP, U_val, R)
## System and Calculation Settings
s = BulkSystem(μ, Potential)
STM_Bias = -0.30

nPts = 1500;
pos_lim = 30;
xs = ys = -pos_lim:pos_lim

# Energy range to calculate dI/dV
ωs = range(0.0, 2.0 + μ, length = nPts)

## Visualize the positions of the potential prior to calculations
# x_positions2 = map(y -> y.loc.v1, Potential)
# y_positions2 = map(y -> y.loc.v2, Potential)
#
# X2 = b1[1] .* x_positions2 + b2[1] .* y_positions2
# Y2 = b1[2] .* x_positions2 + b2[2] .* y_positions2
# scatter(
#     X2,
#     Y2,
#     aspectratio = 1,
#     title = L"\textrm{Position of UCs with Local Potential}",
#     xlabel = L"\textrm{Distance}\, (a_0)",
#     ylabel = L"\textrm{Distance}\, (a_0)",
#     label = L"\textrm{Local Potential}",
#     framestyle = :box,
# )
