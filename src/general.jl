using QuadGK
using SpecialFunctions
using LinearAlgebra

## Parameters
const ħ = 1;
const ryd = 13.606;     # Rydberg constant in eV
const a0 = 0.529;       # Bohr radius in angstroms

const mx = 0.27;        # Effective mass in the x and y direction in m_e -
const my = 0.27;        # anisotropy tuning

# Lattice Basis vectors in Bohr radii
b1 = 2.88 / a0 * [1, 0]
b2 = 2.88 / a0 * [cos(π / 3), sin(π / 3)]
# Area of the unit cell in Bohr radii squared
UC_area = abs.(b1[1] * b2[2] - b1[2] * b2[1]);

# Upper limit of propagator integral
C = 4 * π / UC_area / sqrt(mx * my) * ryd;

## Integration
const ν = 1e-4;         # Relative tolerance for integration
const η = 1e-6;         # Small number for imaginary parts
const NumEvals = 1e5;   # Max number of integrals in quadgk

## Location struct to label bulk unit cells.
struct Location
    v1::Int64   # Coefficient of b1
    v2::Int64   # Coefficient of b2
end

# Function to turn continuous coordinates into the unit cell ones
function locator(x, y)
    v = round.(inv(hcat(b1, b2)) * [x, y])
    return Location(v[1], v[2])
end

# Introduce addition and subtraction for the Location struct
function Base.:+(l1::Location, l2::Location)
    return Location(l1.v1 + l2.v1, l1.v2 + l2.v2)
end

function Base.:-(l1::Location, l2::Location)
    return Location(l1.v1 - l2.v1, l1.v2 - l2.v2)
end

# Potential modification
struct LocalPotential
    U::Float64                  # Additive potential
    loc::Location               # Location of the unit cell
end

# The structure describing the system
struct BulkSystem
    μ::Float64                          # Chemical potential
    potential::Vector{LocalPotential}   # Local potential
end


## Propagator
function Ξ(r::Location, z::ComplexF64)
    position = r.v1 * b1 + r.v2 * b2
    md = √(mx * position[1]^2 + my * position[2]^2)
    if md == 0.0
        return (-1 / C) * (log(Complex(1 - (C / z))))
    else
        return (-2 / C) * besselk(0, (md * √(-z / ryd)))
    end
end

function propagator_matrix(z, Coords::Vector{Location})
    len_coords = length(Coords)
    out = zeros(ComplexF64, len_coords, len_coords)
    for ii = 1:len_coords
        @inbounds for jj = ii:len_coords
            out[ii, jj] = Ξ(Coords[ii] - Coords[jj], z)
            out[jj, ii] = out[ii, jj]
        end
    end
    return out
end

function spectral_bulk(ω, R::Location, s::BulkSystem)
    z = ω + 1im * η
    potential = s.potential
    pristine_spectral = -Ξ(Location(0, 0), ω + 1im * η) / π |> imag
    n_UCs = length(potential)
    prop_mat = propagator_matrix(z, map(x -> x.loc, potential))
    PropVectorR = map(x -> Ξ(x.loc - R, z), potential)
    U = map(x -> x.U, potential) |> Diagonal
    D = U * inv(Diagonal(ones(n_UCs)) .- prop_mat * U)
    correction_spectral =
        -((transpose(PropVectorR)*D*PropVectorR)[1]) / π |> imag
    return (pristine_spectral + correction_spectral)
end

## Potential shapes

function mkRing(numP::Int, U::Float64, R::Float64)
    res = map(
        x -> LocalPotential(
            U,
            locator(
                R * cos(2 * π / numP * x + π / numP),
                R * sin(2 * π / numP * x + π / numP),
            ),
        ),
        1:numP,
    )
    return res
end

# function mkFig8(numP::Int, U::Float64, R::Float64)
#     Potential_Top = map(
#         x -> LocalPotential(
#             U,
#             locator(R * cos(x - π / 5), R * sin(x - π / 5) + R * sin(π / 5)),
#         ),
#         range(0, 7 * π / 5, length = Integer(numP / 2 + 1)),
#     )
#
#     Potential_Bottom = map(
#         x -> LocalPotential(
#             U,
#             locator(R * cos(x - π / 5), -R * sin(x - π / 5) - R * sin(π / 5)),
#         ),
#         range(0, 7 * π / 5, length = Integer(numP / 2 + 1)),
#     )
#     Potential_Bottom = Potential_Bottom[2:end-1]
#     Potential = vcat(Potential_Top, Potential_Bottom)
# end

function mkFig8(
    numP::Int,
    U::Float64,
    R1::Float64,
    R2::Float64,
    d1::Float64,
    d2::Float64,
)
    θ1 = acos(d1 / R1)
    θ2 = acos(d2 / R2)
    arc = 2 * (π - θ1) * R1 + 2 * (π - θ2) * R2
    numP1 = (numP - 2) * 2 * (π - θ1) * R1 / arc + 2 |> round |> Integer
    numP2 = (numP - 2) * 2 * (π - θ2) * R2 / arc + 2 |> round |> Integer

    Potential_Right = map(
        x -> LocalPotential(U, locator(d1 + R1 * cos(x), R1 * sin(x))),
        range(-π + θ1, π - θ1, length = numP1),
    )
    Potential_Left = map(
        x -> LocalPotential(U, locator(-d2 + R2 * cos(x), R2 * sin(x))),
        range(θ2, 2 * π - θ2, length = numP2),
    )
    Potential_Left = Potential_Left[2:end-1]

    Potential = vcat(Potential_Right, Potential_Left)
end
