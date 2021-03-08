include("general.jl")

# Integrand used to calculate the local density in the bulk
function δρ_calc(loc::Location, z, potential::Vector{LocalPotential})
    n_UCs = length(potential)
    prop_mat = propagator_matrix(z, map(x -> x.loc, potential))
    PropVectorR = map(x -> Ξ(x.loc - loc, z), potential)
    U = map(x -> x.U, potential) |> Diagonal
    D = U * inv(Diagonal(ones(n_UCs, n_UCs)) .- prop_mat * U)
    return (transpose(PropVectorR)*D*PropVectorR)[1]
end

function δρR(R::Location, s::BulkSystem)
    μ = s.μ
    potential = s.potential
    f_int(x) = (1 / 2 * π) * δρ_calc(R, μ + 1im * x, potential)
    res = quadgk(f_int, 0, Inf, maxevals = NumEvals, rtol = ν)
    return 2 * real(res[1][1])
end

function spectral_bulk(ω, R::Location, s::BulkSystem)
    μ = s.μ
    potential = s.potential
    pristine_spectral = -Ξ(Location(0, 0), ω + 1im * η) / π |> imag
    correction_spectral = -δρ_calc(R, ω + 1im * η, potential) / π |> imag
    return (pristine_spectral + correction_spectral)
end
