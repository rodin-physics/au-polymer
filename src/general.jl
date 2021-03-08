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
