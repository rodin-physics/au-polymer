include("../src/density.jl")
using ProgressMeter
using Plots
using LaTeXStrings

## System Parameters
μ = 0.46
T = 0.0

## --- Change parameters here ---
exp_energy = -0.33
slice_energy = exp_energy + μ     # Energy measured from band edge

U_val = 0.6                     # Value of potential
numP = 54                         # Number of points
potential_shape = "symCassOval"  #"Ring" OR "symCassOval" OR "asymCassOval"
vis_Potential = true              # Visualise potential shape

# Parameters of potential rings
R = 20.5 / a0            # Radius of circular potential ring

R_sCO = 19.0 / a0        # Radius of circle in Cassini Oval
vert_sCO = 15.8 / a0     # Length of vertical separation in Cassini Oval

R_aCO1 = 20.3 / a0       # Radius of upper ring
R_aCO2 = 18.3 / a0       # Radius of lower ring
vert_aCO1 = 13.5 / a0    # Vertical separation of upper ring
vert_aCO2 = 12.0 / a0    # Vertical separation of lower ring

if potential_shape == "Ring"
    Potential = mk_Ring(numP, U_val, R)
elseif potential_shape == "symCassOval"
    Potential = mk_symCassOval(numP, U_val, R_sCO, vert_sCO)
elseif potential_shape == "asymCassOval"
    Potential = mk_asymCassOval(numP, U_val, R_aCO1, R_aCO2, vert_aCO1, vert_aCO2)
end

## System and Calculation Settings
s = BulkSystem(μ, Potential)

nPts = 1500;
pos_lim = 30;
xs = ys = -pos_lim:pos_lim

# Energy range to calculate dI/dV
ωs = range(0.0, 2.0 + μ, length = nPts)

# Visualize the positions of the potential prior to calculations
if vis_Potential == true
    x_positions2 = map(y -> y.loc.v1, Potential)
    y_positions2 = map(y -> y.loc.v2, Potential)

    X2 = b1[1] .* x_positions2 + b2[1] .* y_positions2
    Y2 = b1[2] .* x_positions2 + b2[2] .* y_positions2
    scatter(
        X2,
        Y2,
        aspectratio = 1,
        title = "Position of UCs with Local Potential",
        xlabel = "Distance (a₀)",
        ylabel = "Distance (a₀)",
        label = "Local Potential",
        framestyle = :box,
        aspect_ratio = 1,
        leg = false,
        xrange = (-60, 60),
        yrange = (-60, 60)
    )
end
