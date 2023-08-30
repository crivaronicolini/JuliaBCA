using Unitful, StaticArrays
import Unitful: Length, Energy, Mass
using DimensionfulAngles
using DimensionfulAngles: Angle

Unitful.preferunits(u"Å")
const LengthConcrete = typeof(1.0u"μm")
const EnergyConcrete = typeof(1.0u"eV")
const MassConcrete = typeof(1.0u"u")
const AngleConcrete = typeof(1.0u"radᵃ")

mutable struct Vec3{T} <: FieldVector{3,T}
  x::T
  y::T
  z::T
end

# StaticArrays.similar_type(::Type{<:Vec3}, ::Type{T}, s::Size{(3,)}) where {T} = Vec3{T}

# function Base.show(io::Base.IO, v::Vec3)
#   println("show Vec3($(v.x), $(v.y), $(v.z))")
# end

struct TrajectoryElement
  E::EnergyConcrete
  pos::Vec3{LengthConcrete}
end

struct EnergyLoss
  En::EnergyConcrete
  Ee::EnergyConcrete
  pos::Vec3{LengthConcrete}
end

abstract type AbstractParticle end

# Particle object. Particles in rustbca include incident ions and material atoms.
# defaults to Particle::new from rustbca
@kwdef mutable struct Particle <: AbstractParticle
  const m::MassConcrete
  const Z::Int
  E::EnergyConcrete
  const Ec::EnergyConcrete
  const Es::EnergyConcrete
  pos::Vec3{LengthConcrete}
  dir::Vec3{Float64}
  pos_old::Vec3{LengthConcrete} = pos
  dir_old::Vec3{Float64} = dir
  pos_origin::Vec3{LengthConcrete} = pos
  energy_origin::EnergyConcrete = E
  asymptotic_deflection::LengthConcrete = 0.0 * u"μm"
  stopped::Bool = false
  left::Bool = false
  incident::Bool = false
  first_step::Bool = incident
  trajectory::Vector{TrajectoryElement} = [TrajectoryElement(E, pos)]
  energies::Vector{EnergyLoss} = [EnergyLoss(0.0u"eV", 0.0u"eV", pos)]
  track_trajectories::Bool = false
  number_collision_events::Int = 0
  backreflected::Bool = false
  interactionindex::Int8 = 1
  weight::Float64 = 1.0
  tag::Int8 = 0
  tracked_vector::Vec3{Float64} = Vec3(0.0, 0.0, 0.0)
  is_from_target::Bool = true
end

function default_incident(m::Mass, Z::Int, E::Energy, Ec::Energy, Es::Energy, x::Length, dir::Vector; track_trajectories=false)
  # @debug "default_incident"
  dir = Vec3(dir...)

  # @assert isapprox(abs(dir.x / dir_mag), 1, atol=eps(Float64), rtol=0) "Input error: incident direction cannot round to exactly (1, 0, 0) due to gimbal lock. Use a non-zero y-component."
  @assert !(normalized(dir) ≈ [1, 0, 0]) "Input error: incident direction cannot round to exactly (1, 0, 0) due to gimbal lock. Use a non-zero y-component."

  # dir_mag = norm(dir)
  @assert E > zero(E) "Input error: incident energy $E; must be greater than zero."

  Particle(m=Float64(m), Z=Z, E=Float64(E), Ec=Float64(Ec), Es=Float64(Es), pos=Vec3(ustrip(Float64(x)), 0.0, 0.0) .* unit(x), dir=normalized(dir), incident=true, track_trajectories=track_trajectories, is_from_target=false)
end

function add_trajectory!(particle::Particle)
  # @debug "add_trajectory!" particle.trajectory particle.E particle.pos
  if particle.track_trajectories
    push!(particle.trajectory, TrajectoryElement(particle.E, particle.pos))
  end
end

# MOVED to main.jl
# # If `track_energy_losses`, add the most recent electronic and nuclear energy loss terms and (x, y, z) to the energy loss tracker.
# function energy_loss!(particle::Particle, options::Options, En::EnergyConcrete, Ee::EnergyConcrete)
#   if particle.incident && options.track_energy_losses
#     push!(particle.energies, EnergyConcreteLoss(Ee, En, particle.pos))
#   end
# end


# CONCRETE TYPES DAN CERO CODE WARNTYPES
# @kwdef mutable struct Particle4
#   const m::typeof(1.0u"u")
#   const Z::Int
#   E::typeof(1.0u"eV")
#   const Ec::typeof(1.0u"eV")
#   const Es::typeof(1.0u"eV")
#   pos::Vec3{typeof(1.0u"μm")}
#   dir::Vec3{Float64}
# end
# p2 = Particle4(1u"u", 2, 3.3u"eV", 4u"eV", 5u"eV", Vec3(1.2, 2.0, 1.0)u"μm", Vec3(1.0, 0.0, 0.0))

# function get_momentum(particle::Particle4)
#   speed = sqrt(2.0 * particle.E / particle.m)
#   return particle.m * speed .* particle.dir
# end

function get_momentum(particle::Particle)
  speed = sqrt(2.0 * particle.E / particle.m)
  return particle.m * speed .* particle.dir
end

function dot(u::FieldVector, v::FieldVector)
  sum(u .* v)
end

function norm(v::FieldVector)
  √(sum(v .^ 2))
end

function norm(v)
  √(sum(v .^ 2))
end

function normalized(v::FieldVector)
  v ./ norm(v)
end

#Rotate a particle by deflection ψ at an azimuthal angle ϕ
function rotate!(particle::Particle, ψ::AngleConcrete, ϕ::AngleConcrete)
  # @debug "rotate!"
  #Particle direction update formula (2) from the original TRIDYN paper, see Moeller and Eckstein 1988
  cosα, cosβ, cosγ = particle.dir
  sinα = √(1 - cosα^2)
  sinϕ, cosϕ = sincos(ϕ)
  sinψ, cosψ = sincos(ψ)

  particle.dir = normalized(Vec3(cosψ * cosα + sinψ * cosϕ * sinα,
    cosψ * cosβ - (cosϕ * cosα * cosβ - sinϕ * cosγ) * sinψ / sinα,
    cosψ * cosγ - (cosϕ * cosα * cosγ + sinϕ * cosβ) * sinψ / sinα))
end

function advance!(particle::Particle, mfp::LengthConcrete, asymptotic_deflection::LengthConcrete)
  # @debug "advance!"
  if particle.E > particle.Ec
    add_trajectory!(particle)
  end

  #Update previous position
  particle.pos_old = particle.pos

  #In order to keep average denisty constant, must add back previous asymptotic deflection
  # TODO ver que onda el codigo original, esta sumando dist con angulos
  # distance_traveled = mfp + particle.asymptotic_deflection - asymptotic_deflection
  # voy a emparcharlo con lo que supongo que es
  distance_traveled = mfp + particle.asymptotic_deflection - asymptotic_deflection

  #dir has been updated, so use previous direction to advance in space
  particle.pos += particle.dir_old * distance_traveled
  particle.asymptotic_deflection = asymptotic_deflection

  #Update previous direction
  particle.dir_old = particle.dir

  return distance_traveled
end

# change particle direction by refraction
function surface_refraction!(particle::Particle, normal::Vec3, Es::EnergyConcrete)
  # @debug "surface_refraction!"
  E = particle.E
  cosθ = dot(particle.dir, normal)

  u = √(E / (E + Es)) .* particle.dir .+ ((-√E * cosθ + √(E * cosθ^2 + Es)) / √(E + Es)) * normal
  particle.dir = u
  particle.E += Es
end


# Calculate the refraction angle based on the surface binding energy of the material.
function refraction_angle(cosθ, energy_old, energy_new)
  # @debug "refraction_angle"
  cosθ = cosθ > 1 ? sign(cosθ) : cosθ
  sinθ_0 = sqrt(1 - cosθ^2)
  sinθ_1 = sinθ_0 * sqrt(energy_old / energy_new)
  Δsinθ = asin(sinθ_1) - asin(sinθ_0)
  @assert !isnan(Δsinθ) "Numerical error: refraction returned NaN."
  sign = -sign(cosθ)
  return sign * Δsinθ
end
