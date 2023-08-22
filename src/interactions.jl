
abstract type InteractionPotential end
struct Moliere <: InteractionPotential end
struct KR_C <: InteractionPotential end
struct Lenz_Jensen <: InteractionPotential end
struct ZBL <: InteractionPotential end
struct TRIDYN <: InteractionPotential end
struct WW <: InteractionPotential end
struct Coulomb <: InteractionPotential end
struct Lennard_Jones_12_6 <: InteractionPotential
  σ::Float64
  ε::Float64
end
struct Lennard_Jones_65_6 <: InteractionPotential
  σ::Float64
  ε::Float64
end
struct Morse <: InteractionPotential
  D::Float64
  α::Float64
  r₀::Float64
end
struct aber <: InteractionPotential
  s::Int
end


#Analytic solutions to outermost root of the interaction potential.
crossing_point_doca(pot::Union{Lennard_Jones_12_6,Lennard_Jones_65_6}) = pot.σ
crossing_point_doca(pot::Morse) = pot.r₀ - log(2) / pot.α
crossing_point_doca(::Type{WW}) = 50.0 * u"A"
crossing_point_doca(::Type{InteractionPotential}) = 10.0 * u"A"

# TODO more types with union
function screened_coulomb(r::LengthConcrete, a::LengthConcrete, Za::Int, Zb::Int, interaction_potential::Moliere)
  # @debug "screened_coulomb"
  Za * Zb * (e^2) / (4π * ε0 * r * ϕ_screening(upreferred(r / a), interaction_potential))
end

# Coulombic interaction potential.
function coulomb(r::LengthConcrete, Za::Int, Zb::Int)
  # @debug "coulomb"
  Za * Zb * (e^2) / (4π * ε0 * r)
end

# todo rest
interaction_potential_fn(r::LengthConcrete, a::LengthConcrete, Za::Int, Zb::Int, interaction_potential::Moliere) = screened_coulomb(r, a, Za, Zb, interaction_potential)

interaction_potential_fn(r::LengthConcrete, ::Real, Za::Int, Zb::Int, ::Type{Coulomb}) = coulomb(r, Za, Zb)

moliere(xi::Real) = 0.35 * exp(-0.3 * xi) + 0.55 * exp(-1.2 * xi) + 0.10 * exp(-6.0 * xi)

#TODO
# Screening functions for screened-coulomb interaction potentials.
ϕ_screening(xi::Real, ::Type{Moliere}) = moliere(xi)
# ϕ_screening(xi::Real, ::Type{KR_C}) = kr_c(xi)
# ϕ_screening(xi::Real, ::Type{ZBL}) = zbl(xi)
# ϕ_screening(xi::Real, ::Type{LENZ_JENSEN}) = lenz_jensen(xi)
# ϕ_screening(xi::Real, ::Type{TRIDYN}) = kr_c(xi)

diff_moliere(xi::Real) = -0.35 * 0.3 * exp(-0.3 * xi) - 0.55 * 1.2 * exp(-1.2 * xi) - 0.10 * 6.0 * exp(-6.0 * xi)

# First derivative w.r.t. `r` of the screening functions for screened-coulomb interaction potentials.
dϕ_screening(xi::Real, ::Type{Moliere}) = diff_moliere(xi)
# dϕ_screening(xi::Real, ::Type{KR_C}) = diff_kr_c(xi)
# dϕ_screening(xi::Real, ::Type{ZBL}) = diff_zbl(xi)
# dϕ_screening(xi::Real, ::Type{LENZ_JENSEN}) = diff_lenz_jensen(xi)
# dϕ_screening(xi::Real, ::Type{TRIDYN}) = diff_kr_c(xi)

#TODO
function screening_length(Za::Int, Zb::Int, ::Type{Moliere})
  # @debug "screening_length"
  # returns length
  0.8853 * a_0 * (√Za + √Zb)^(-2 / 3)
end

# TODO
# Coefficients of inverse-polynomial interaction potentials.
# function polynomial_coefficients(relative_energy::Real, impact_parameter::Real, ::Type{Lennard_Jones_12_6)
# end

# #Interaction potential between two particles a and b at a distance `r`.
# function interaction_potential(r::LengthConcrete, a::Float64, Za::Float64, Zb::Float64, pot::Lennard_Jones_12_6)
#   lennard_Jones_12_6_potential(r, pot.σ, pot.ε)
# end
#
# function interaction_potential(r::LengthConcrete, a::Float64, Za::Float64, Zb::Float64, pot::Lennard_Jones_65_6)
#   lennard_Jones_65_6_potential(r, pot.σ, pot.ε)
# end
#
# function interaction_potential(r::LengthConcrete, a::Float64, Za::Float64, Zb::Float64, pot::Morse)
#   morse_potential(r, pot.D, pot.α, pot.r₀)
# end
#
# function interaction_potential(r::LengthConcrete, a::Float64, Za::Float64, Zb::Float64, pot::WW)
#   tungsten_tungsten_cubic_spline(r)
#
# end
# function interaction_potential(r::LengthConcrete, a::Float64, Za::Float64, Zb::Float64, pot::Coulomb)
#   coulomb_potential(r, Za, Zb)
# end

#Analytic energy threshold above which the distance of closest approach function is guaranteed to have only one root.
function energy_threshold_single_root(pot::Type{Union{Lennard_Jones_12_6,Lennard_Jones_65_6,Morse,WW,Coulomb}})
  Inf
end

function energy_threshold_single_root(pot::Union{Moliere,KR_C,Lenz_Jensen,ZBL,TRIDYN})
  0.0
end

first_screening_radius(::Type{Moliere}) = 0.3

# Distance of closest approach function for screened coulomb potentials.
# Nonlinear equation to determine distance of closest approach
# TODO expando to other interaction potentials
doca_function(x0::Real, β::Real, reduced_energy, interaction_potential::Type{Moliere}) = x0 - ϕ_screening(x0, interaction_potential) / reduced_energy - β^2 / x0

# First derivative w.r.t. `r` of the distance of closest approach function for screened coulomb potentials.
diff_doca_function(x0::Real, β::Real, reduced_energy::EnergyConcrete, interaction_potential::InteractionPotential) = β^2 / x0^2 - dϕ_screening(x0, interaction_potential) / reduced_energy + 1

# TODO rest of interaction potentials 
function distance_of_closest_approach_function(a::LengthConcrete, Za::Int, Zb::Int, relative_energy::EnergyConcrete, impact_parameter::LengthConcrete, interaction_potential::Type{Moliere})
  # @debug "distance_of_closest_approach_function"
  reduced_energy = LINDHARD_REDUCED_ENERGY_PREFACTOR * a * relative_energy / (Za * Zb)
  β = impact_parameter / a |> NoUnits
  f(r::LengthConcrete) = doca_function(upreferred(r / a), β, reduced_energy, interaction_potential)
  return f
end
