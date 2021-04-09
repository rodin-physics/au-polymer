include("../src/general.jl")
using ProgressMeter
using LaTeXStrings

## System Parameters
μ = 0.46
slice_energy = μ + 1.3

## Ring parameters
R = 21.8 / a0

# Local Potential
U_val = 0.5
numP = 48
Potential_Ring = mkRing(numP, U_val, R)

# Potential =
#     mkFig8(80, 0.1, 19 / a0, 19 / a0, (5.5) / a0, 5.5/ a0)

## Figure 8 parameters
# R = 21.8 / a0
#
# # Local Potential
# U_val = 0.5
# numP = 80
# Potential_Fig8 = mkFig8(numP, U_val, R)
## System and Calculation Settings
sRing = BulkSystem(μ, Potential_Ring)
# sFig8 = BulkSystem(μ, Potential_Fig8)

nPts = 1500;
pos_lim = 45;
xs = ys = -pos_lim:pos_lim

# Energy range to calculate dI/dV
ωs = range(0.0, 2.5 + μ, length = nPts)

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

#
# #
# Potential =
#     mkFig8(50, 0.1, 19 / a0, 19 / a0, (5.5) / a0, 5.5/ a0)
#
# # Potential = mkRing(30, 0.1, 20.5 / a0)
# #
# # gr()
# x_positions2 = map(y -> y.loc.v1, Potential)
# y_positions2 = map(y -> y.loc.v2, Potential)
#
# X2 = b1[1] .* x_positions2 + b2[1] .* y_positions2
# Y2 = b1[2] .* x_positions2 + b2[2] .* y_positions2
# scatter(
#     X2,
#     Y2,
#     aspectratio = 1,
#     # title = L"\textrm{Position of UCs with Local Potential}",
#     # xlabel = L"\textrm{Distance}\, (a_0)",
#     # ylabel = L"\textrm{Distance}\, (a_0)",
#     # label = L"\textrm{Local Potential}",
#     # framestyle = :box,
# )
# #
#
# Potential =
#     mkFig8(50, 0.1, 19 / a0, 19 / a0, (5.5) / a0, 5.5/ a0)

# Potential = mkRing(30, 0.1, 20.5 / a0)
#
# # gr()
# x_positions2 = map(y -> y.loc.v1, Potential)
# y_positions2 = map(y -> y.loc.v2, Potential)
#
# X2 = b1[1] .* x_positions2 + b2[1] .* y_positions2
# Y2 = b1[2] .* x_positions2 + b2[2] .* y_positions2
# scatter(
#     X2,
#     Y2,
#     aspectratio = 1,
#     # title = L"\textrm{Position of UCs with Local Potential}",
#     # xlabel = L"\textrm{Distance}\, (a_0)",
#     # ylabel = L"\textrm{Distance}\, (a_0)",
#     # label = L"\textrm{Local Potential}",
#     # framestyle = :box,
# )
# #
