using QuadGK
using SpecialFunctions
using LinearAlgebra

## Parameters
const ryd = 13.606;     # Rydberg constant in eV
const a0 = 0.529;       # Bohr radius in angstroms

const mx = 0.27;        # Effective mass in x direction in m_e
const my = 0.27;        # Effective mass in y direction in m_e

# Lattice Basis vectors in Bohr radii
b1 = 2.88 / a0 * [1, 0]
b2 = 2.88 / a0 * [cos(π / 3), sin(π / 3)]
# Area of the unit cell in Bohr radii squared
UC_area = abs.(b1[1] * b2[2] - b1[2] * b2[1]);

# Upper limit of propagator integral
C = 4 * π / UC_area / sqrt(mx * my) * ryd;

## Integration parameters
const ν = 1e-4;         # Relative tolerance for integration
const η = 1e-6;         # Small number for imaginary parts
const NumEvals = 1e5;   # Max number of integrals in quadgk

# Location struct to label bulk unit cells.
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


## Propagator functions
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

## Shapes of Potential Profiles
function mk_Ring(numP::Int, U::Float64, R::Float64)
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

function mk_symCassOval(numP::Int, U::Float64, R1::Float64, vert::Float64)
    Potential_Top = map(
        x -> LocalPotential(
            U,
            locator(R1 * cos(x - π / 5), R1 * sin(x - π / 5) + vert * sin(π / 5)),
        ),
        range(0, 7 * π / 5, length = Integer(numP / 2 + 1)),
    )

    Potential_Bottom = map(
        x -> LocalPotential(
            U,
            locator(R1 * cos(x - π / 5), -R1 * sin(x - π / 5) - vert * sin(π / 5)),
        ),
        range(0, 7 * π / 5, length = Integer(numP / 2 + 1)),
    )
    Potential_Bottom = Potential_Bottom[2:end-1]
    Potential = vcat(Potential_Top, Potential_Bottom)
    deleteat!(Potential, findall(x->x== LocalPotential(0.6, Location(-5,-1)),Potential))
    deleteat!(Potential, findall(x->x== LocalPotential(0.6, Location(6,-1)),Potential))
end


function mk_asymCassOval(numP::Int, U::Float64, R1::Float64, R2::Float64, vert1::Float64, vert2::Float64)
    Potential_Top = map(
        x -> LocalPotential(
            U,
            locator(R1 * cos(x - π / 5), R1 * sin(x - π / 5) + vert1 * sin(π / 5)),
        ),
        range(0, 7 * π / 5, length = Integer(numP / 2 + 3)),
    )

    Potential_Bottom = map(
        x -> LocalPotential(
            U,
            locator(R2 * cos(x - π / 5), -R2 * sin(x - π / 5) - vert2 * sin(π / 5)),
        ),
        range(0, 7 * π / 5, length = Integer(numP / 2 - 1)),
    )
    Potential_Bottom = Potential_Bottom[2:end-1]
    Potential = vcat(Potential_Top, Potential_Bottom)
end
