using Distributed
using Plots
using DelimitedFiles

@everywhere include("/Users/harshitramahalingam/Documents/Senior Semester 2/Capstone and Research Seminar/BP_Quantum_Dots/src/density.jl")

## Read in data - change file path accordingly
data122 = readdlm("/Users/harshitramahalingam/Documents/Senior Semester 2/Capstone and Research Seminar/BP_Quantum_Dots/dIdV ExpData/Bias-Spectroscopy00122.dat")

data357 = readdlm("/Users/harshitramahalingam/Documents/Senior Semester 2/Capstone and Research Seminar/BP_Quantum_Dots/dIdV ExpData/Bias-Spectroscopy00357.dat")

data362 = readdlm("/Users/harshitramahalingam/Documents/Senior Semester 2/Capstone and Research Seminar/BP_Quantum_Dots/dIdV ExpData/Bias-Spectroscopy00362.dat")

data406 = readdlm("/Users/harshitramahalingam/Documents/Senior Semester 2/Capstone and Research Seminar/BP_Quantum_Dots/dIdV ExpData/Bias-Spectroscopy00406.dat")

## Data processing
V_data122 = data122[2:end,1]
dIdV_data122 = data122[2:end,3]

V_data357 = data357[2:end,1]
dIdV_data357 = data357[2:end,3]

V_data362 = data362[2:end,1]
dIdV_data362 = data362[2:end,3]

V_data406 = data406[2:end,1]
dIdV_data406 = data406[2:end,3]

## Parameters - change here
ind= 380
μ = 0.46
T = 0.0
R = 21.8 / a0
ring_R = [R]
numP = 24                               # Number of points in ring

Hoppings = Hop[]
LocalizedStates = LocalizedState[]

nPts = 1500;
ωs = range(-μ, 2.0+μ, length = nPts)
tip_pot = 0.3                           # Value for STM tip potential
shiftby = 0.1                           # Value to shift x-axis by


## Animation loop
anim = @animate for ii = range(0.1, 1.5, step = 0.1)

    Potential = [(
        map(
            y -> map(
                x -> LocalPotential(
                    ii,
                    locator(
                        y * cos(π / (numP / 2) * x + π / numP),
                        y * sin(π / (numP / 2) * x + π / numP),
                    ),
                ),
                1:numP,
            ),
            ring_R,
        )...
    )...]

    s = BulkSystem(μ, T, LocalizedStates, Hoppings, Potential)


    dI_dV = pmap(x -> spectral_bulk(x, Location(0, 0), s), ωs)
    dI_dV_with_tip = pmap(x -> spectral_bulk_add_potential(x, Location(0, 0), s, tip_pot), ωs)

    plot(
        ωs.-μ.+shiftby,
        dI_dV,
        xlabel = "Energy ω (eV)",
        ylabel = "Spectral Fn",
        label = "U = $ii eV",
        framestyle = :box,
        linecolor = :blue,
        xrange = (-0.5, 2.0),
        title = "Shift of +0.1eV"
    )
    plot!(
        ωs.-μ.+shiftby,
        dI_dV_with_tip,
        label = "U = $ii eV; TP = $tip_pot eV",
    )
    factor = dI_dV[end]/((dIdV_data406[1]/dIdV_data406[ind])-1.0)
    plot!(V_data406[1:ind], (dIdV_data406[1:ind]./dIdV_data406[ind].-1.0).*factor, linecolor = :black, label = "Exp. - Ring Center")
end

# saving gif
gif(anim, "dIdV_anim.gif", fps = 2)
