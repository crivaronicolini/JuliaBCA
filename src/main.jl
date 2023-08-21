using Roots

@enum MeanFreePathModel LIQUID GASEOUS
@derived_dimension StoppingPower (Unitful.𝐋^4 * Unitful.𝐌 * Unitful.𝐓^(-2))
function dist(u::Vec3, v::Vec3)
  √sum((u - v) .^ 2)
end

function scattering_integrand_mw(x::Real, β::Real, reduced_energy, interaction_potential::Type{T} where {T<:InteractionPotential})
  (1 - ϕ_screening(x, interaction_potential) / (x * reduced_energy) - β^2 / x^2)^(-1 / 2)
end

abstract type ScatteringIntegral end
struct Mendenhall_weller <: ScatteringIntegral end

# function Mendenhall_weller(Za::Int, Zb::Int, ma::Mass, mb::Mass, E0::Energy, impactparameter::Length, x0::Length, interaction_potential::Union{Type{Moliere}})
function Mendenhall_weller(Za::Int, Zb::Int, ma::Mass, mb::Mass, E0::Energy, impactparameter::Length, x0::Real, interaction_potential::Type{T} where {T<:InteractionPotential})
  # Lindhard screening length and reduced energy
  a = screening_length(Za, Zb, interaction_potential)
  relative_energy = E0 * mb / (ma + mb)
  reduced_energy = LINDHARD_REDUCED_ENERGY_PREFACTOR * a * relative_energy / (Za * Zb)
  β = impactparameter / a |> NoUnits

  f = scattering_integrand_mw
  # Scattering integral quadrature from Mendenhall and Weller 2005
  λ0 = ((1 / 2) + (β^2 / (2x0^2)) - (dϕ_screening(x0, interaction_potential) / 2reduced_energy))^(-1 / 2)
  # rust bca uses different formula and constants TODO check
  α = (1 / 12) *
      (1 + λ0
       + 5 * (0.4206 * f(x0 / 0.9072, β, reduced_energy, interaction_potential)
              +
              0.9072 * f(x0 / 0.4206, β, reduced_energy, interaction_potential)))


  return π * (1 - β * α / x0)
end

@kwdef struct Options
  name::String
  track_trajectories::Bool = false
  track_recoils::Bool = false
  track_recoil_trajectories::Bool = false
  write_buffer_size::Int = 3
  weak_collision_order::Int = 3
  suppress_deep_recoils::Bool = false
  high_energy_free_flight_paths::Bool = false
  electronic_stopping_mode::ElectronicStoppingMode = LOWENERGYNONLOCAL
  mean_free_path_model::MeanFreePathModel = LIQUID
  interaction_potential::Vector{Vector{Type{T} where {T<:InteractionPotential}}} = [[Moliere]]
  scattering_integral::Type{S} where {S<:ScatteringIntegral} = Mendenhall_weller
  # root_finder::Vector{Vector{Function}}
  num_threads::Int = 8
  num_chunks::Int = 8
  use_hdf5::Bool = false
  track_displacements::Bool = false
  track_energy_losses::Bool = false
  accelerated_ions::Bool = false
end

struct CollisionResult
  θ::Angle
  ψ::Angle
  ψ_recoil::Angle
  E_recoil::Energy
  asymptotic_deflection::Length
  normalized_distance_of_closest_approach::Real
end

struct CollisionGeometry
  ϕ_azimuthal::Angle
  impactparameter::Length
  mfp::Length
end

# FROM particle.jl
# If `track_energy_losses`, add the most recent electronic and nuclear energy loss terms and (x, y, z) to the energy loss tracker.
function energy_loss!(particle::Particle, En::Energy, Ee::Energy, options::Options)
  @debug "energy_loss!"
  if particle.incident && options.track_energy_losses
    push!(particle.energies, EnergyLoss(Ee, En, particle.pos))
  end
end

# CollisionGeometry(ϕs::Vector{Real}, impactparameters::Vector{Real}, mfp::Real) = [CollisionGeometry(ϕs[k], impactparameters[k], mfp) for k in 1:length(ϕs)]


# For a particle in a material, determine the mean free path and choose the azimuthal angle and impact parameter.
# The mean free path can be exponentially distributed (gaseous) or constant (amorphous solid/liquid). Azimuthal angles are chosen uniformly. Impact parameters are chosen for collision partners distributed uniformly on a disk of density-dependent radius.
function determine_mfp_phi_impact_parameter(particle1::Particle, target::Target, options::Options)
  @info "determine_mfp_phi_impact_parameter"

  mfp = meanfreepath(particle1.pos, target)

  # Each weak collision gets its own aziumuthal angle in annuli around collision point
  collision_order = options.weak_collision_order
  ϕs_azimuthal = 2π * rand(collision_order) * u"radᵃ"

  if options.high_energy_free_flight_paths
    #TODO
    @warn "high_energy_free_flight_paths is not implemented"

  else
    #If not using free flight paths, use weak collision model
    pmax = mfp / √π

    #cylindrical geometry
    impactparameters = pmax .* sqrt.(rand(collision_order))

    #Atomically rough surface - scatter initial collisions
    if particle1.first_step
      mfp *= rand()
      particle1.first_step = false
    end

    if options.mean_free_path_model == GASEOUS
      mfp *= -log(rand())
    end

    return [CollisionGeometry(ϕs_azimuthal[k], impactparameters[k], mfp) for k in eachindex(ϕs_azimuthal)]
    # return CollisionGeometry(ϕs_azimuthal, impactparameters, mfp)
  end
end

# For a particle in a material, determine the mean free path and choose the azimuthal angle and impact parameter.
# The mean free path can be exponentially distributed (gaseous) or constant (amorphous solid/liquid). Azimuthal angles are chosen uniformly. Impact parameters are chosen for collision partners distributed uniformly on a disk of density-dependent radius.
function determine_mfp_phi_impact_parameter_deterministic(particle1::Particle, target::Target, options::Options)
  @info "determine_mfp_phi_impact_parameter"

  mfp = meanfreepath(particle1.pos, target)

  # Each weak collision gets its own aziumuthal angle in annuli around collision point
  collision_order = options.weak_collision_order
  ϕs_azimuthal = 2π * fill(0.2, length(collision_order)) * u"radᵃ"

  if options.high_energy_free_flight_paths
    #TODO
    @warn "high_energy_free_flight_paths is not implemented"

  else
    #If not using free flight paths, use weak collision model
    pmax = mfp / √π

    #cylindrical geometry
    impactparameters = pmax .* sqrt.(fill(0.2, length(collision_order)))

    #Atomically rough surface - scatter initial collisions
    if particle1.first_step
      mfp *= 0.2
      particle1.first_step = false
    end

    if options.mean_free_path_model == GASEOUS
      mfp *= -log(0.2)
    end

    return [CollisionGeometry(ϕs_azimuthal[k], impactparameters[k], mfp) for k in eachindex(ϕs_azimuthal)]
    # return CollisionGeometry(ϕs_azimuthal, impactparameters, mfp)
  end
end

function distance_of_closest_approach(particle1::Particle, particle2::Particle, collisiongeometry::CollisionGeometry, options::Options)
  @info "distance_of_closest_approach"
  Za, Ma = particle1.Z, particle1.m
  Zb, Mb = particle2.Z, particle2.m
  E0 = particle1.E

  relative_energy = E0 * Mb / (Ma + Mb)
  p = collisiongeometry.impactparameter

  interaction_potential = options.interaction_potential[particle1.interactionindex][particle2.interactionindex]

  # TODO coulomb tiene una parte especial

  # inline newton and friends, use Roots.jl
  a = screening_length(Za, Zb, interaction_potential)
  f = distance_of_closest_approach_function(a, Za, Zb, relative_energy, p, interaction_potential)

  try
    x0 = find_zero(f, p)
    return upreferred(x0 / uconvert(u"Å", a)) #adimensional
  catch
    error("failed to calculate closest approach on find_zero")
  end
end

#For a particle in a material, and for a particular binary collision geometry, choose a species for the collision partner.
function choose_collision_partner(particle::Particle, target::Target, collisiongeometry::CollisionGeometry, options::Options)
  x, y, z = particle.pos
  ϕ, impactparameter, mfp = collisiongeometry.ϕ_azimuthal, collisiongeometry.impactparameter, collisiongeometry.mfp

  sinϕ, cosϕ = sincos(ϕ)
  cosdir = cosx, cosy, cosz = particle.dir
  sinx = √(1 - cosx^2)

  recoil = Vec3(x + mfp * cosx - impactparameter * cosϕ * sinx,
    y + mfp * cosy - impactparameter * (sinϕ * cosz - cosϕ * cosy * cosx) / sinx,
    z + mfp * cosz + impactparameter * (sinϕ * cosy - cosϕ * cosx * cosz) / sinx)

  species_index, Z, m, Ec, Es, interactionindex = choose(recoil, target)
  newparticle = Particle(m=m, Z=Z, E=0.0 * Ec, Ec=Ec, Es=Es, pos=recoil, dir=cosdir, track_trajectories=options.track_recoil_trajectories, interactionindex=interactionindex, weight=particle.weight, tag=particle.tag, tracked_vector=particle.tracked_vector)

  return species_index, newparticle
end

#For a particle in a material, and for a particular binary collision geometry, choose a species for the collision partner.
function choose_collision_partner_deterministic(particle::Particle, target::Target, collisiongeometry::CollisionGeometry, options::Options)
  x, y, z = particle.pos
  ϕ, impactparameter, mfp = collisiongeometry.ϕ_azimuthal, collisiongeometry.impactparameter, collisiongeometry.mfp

  sinϕ, cosϕ = sincos(ϕ)
  cosdir = cosx, cosy, cosz = particle.dir
  sinx = √(1 - cosx^2)

  recoil = Vec3(x + mfp * cosx - impactparameter * cosϕ * sinx,
    y + mfp * cosy - impactparameter * (sinϕ * cosz - cosϕ * cosy * cosx) / sinx,
    z + mfp * cosz + impactparameter * (sinϕ * cosy - cosϕ * cosx * cosz) / sinx)

  species_index, Z, m, Ec, Es, interactionindex = choose_deterministic(recoil, target)
  newparticle = Particle(m=m, Z=Z, E=0.0 * Ec, Ec=Ec, Es=Es, pos=recoil, dir=cosdir, track_trajectories=options.track_recoil_trajectories, interactionindex=interactionindex, weight=particle.weight, tag=particle.tag, tracked_vector=particle.tracked_vector)

  return species_index, newparticle
end

# Calculate the binary collision result of two particles for a given binary collision geometry. If the calculation fails, return an Error.
function calculate_binary_collision(particle1::Particle, particle2::Particle, collisiongeometry::CollisionGeometry, options::Options)
  @info "calculate_binary_collision"
  Za, Zb = particle1.Z, particle2.Z
  ma, mb = particle1.m, particle2.m
  E0 = particle1.E

  interaction_potential = options.interaction_potential[particle1.interactionindex][particle2.interactionindex]
  scattering_integral = options.scattering_integral

  a = screening_length(Za, Zb, interaction_potential)
  x0 = distance_of_closest_approach(particle1, particle2, collisiongeometry, options)

  θ = scattering_integral(Za, Zb, ma, mb, E0, collisiongeometry.impactparameter, x0, interaction_potential) * u"radᵃ"
  if θ.val == NaN
    error("θ is NaN, integration error")
  end

  # See Eckstein 1991 for details on center of mass and lab frame angles
  asymptotic_deflection = typeof(interaction_potential) != Coulomb ? x0 * a * sin(θ / 2) : 0.0u"μm"
  ψ = abs(atan(sin(θ), (ma / mb + cos(θ)))) * u"radᵃ"
  ψ_recoil = abs(atan(sin(θ), 1 - cos(θ))) * u"radᵃ"
  E_recoil = E0 * (sin(θ / 2))^2 * (4ma * mb) / (ma + mb)^2

  return CollisionResult(θ, ψ, ψ_recoil, E_recoil, asymptotic_deflection, x0)
end


# Oen-Robinson local electronic energy loss for a collision between particles a and b.
function oen_robinson_loss(Za::Int, Zb::Int, Se::StoppingPower, x0::Real, interaction_potential::Type{T} where {T<:InteractionPotential})
  @info "oen_robinson_loss"
  #Oen-Robinson local electronic stopping power
  a = screening_length(Za, Zb, interaction_potential)

  #d1 is the first (smallest) interior constant of the screening function
  d1 = first_screening_radius(interaction_potential)

  # TODO ISSUE? no tiene las ctes numericas que tiene el paper
  d1^2 * Se * exp(-d1 * x0) / (2π * a^2)
end


function update_particle_energy!(particle::Particle, target::Target, distance_traveled::Length, E_recoil::Energy, x0::Real, strong_collision_Z::Int, strong_collision_index::Int, options::Options)
  @info "update_particle_energy!"
  #If particle energy  drops below zero before electronic stopping calcualtion, it produces NaNs
  @assert !isnan(E_recoil) "Numerical error: recoil Energy is NaN. E0 = $(particle.E) Za = $(particle.Z) Ma = $(particle.m) x0 = $x0 Zb = $strong_collision_Z delta-x = $distance_traveled"

  particle.E -= E_recoil
  @assert !isnan(particle.E) "Numerical error: particle energy is NaN following collision."
  if particle.E < 0u"eV"
    particle.E = 0.0u"eV"
  end

  pos = particle.pos
  e_stop_mode = options.electronic_stopping_mode

  if inside(particle.pos, target)
    interaction_potential = options.interaction_potential[particle.interactionindex][property_atpos(:interactionindex, pos, target)[strong_collision_index]]

    e_stop_powers = electronic_stopping_cross_sections(particle, target, e_stop_mode)
    # ISSUE? El original usa n = number density, no electronic density
    # lo cambié siguiendo wikipedia
    # ADDENDUM no dan las dimensiones con eso, vuelvo a lo que tenia
    # PS: sigues sin dar, por las unidades de n. Multiplico por Å para que de
    n = property_atpos(:density, pos, target) * u"Å"

    #TODO electronic stopping mode
    ΔE = e_stop_mode == INTERPOLATED ? sum(n .* e_stop_powers) :
         e_stop_mode == LOWENERGYLOCAL ? sum(n .* e_stop_powers) * distance_traveled :
         e_stop_mode == LOWENERGYNONLOCAL ? oen_robinson_loss(particle.Z, strong_collision_Z, e_stop_powers[strong_collision_index], x0, interaction_potential) :
         e_stop_mode == LOWENERGYEQUIPARTITION ? 0.5sum(n .* e_stop_powers) * distance_traveled + 0.5oen_robinson_loss(particle.Z, strong_collision_Z, e_stop_powers[strong_collision_index], x0, interaction_potential) : error("No electronic stopping mode supplied")

    @info ΔE
    particle.E += -ΔE

    energy_loss!(particle, E_recoil, ΔE, options)

  elseif E_recoil > 0u"eV"
    energy_loss!(particle, E_recoil, 0u"eV", options)
  end
end

function boundary_condition_planar!(particle::Particle, target::Target)
  @info "boundary_condition_planar"
  pos = particle.pos
  E = particle.E
  Ec = particle.Ec

  if !inside_simulation_boundary(pos, target)
    particle.left = true
    add_trajectory!(particle)
  end

  if (E < Ec) && !particle.left && inside_simulation_boundary(pos, target)
    particle.stopped = true
    add_trajectory!(particle)
  end
end


function singleionbca(particle::Particle, target::Target, options::Options)
  particles = [particle]
  particle_output = Particle[]
  particle_index = length(particles)

  while particle_index > 0
    @info "start new loop"
    @info "particle index $particle_index"
    @info "particle len $(length(particles))"
    particle1 = pop!(particles)
    while !particle1.stopped && !particle1.left
      @info "while !particle1.stopped && !particle1.left"

      collisiongeometries = determine_mfp_phi_impact_parameter(particle1, target, options)

      if options.accelerated_ions
        if !particle1 in target
          point = closest_point(particle1.pos, target)
          distance_to_target = dist(point, particle1.pos)
        else
          distance_to_target = 0.0
        end
      end

      total_energy_loss = 0.0u"eV"
      total_asymptotic_deflection = 0.0u"μm"
      normalized_distance_of_closest_approach = 0.0
      strong_collision_Z = 0
      strong_collision_index = 1

      # collision loop
      # aca pone un take(weak_collison_order +1)
      for (k, collisiongeometry) in enumerate(collisiongeometries)
        @info "for (k, collisiongeometry) in enumerate(collisiongeometries)"

        species_index, particle2 = choose_collision_partner(particle1, target, collisiongeometry)
        #If recoil location is inside, proceed with binary collision loop
        if inside(particle2.pos, target) && inside_energybarrier(particle1.pos, target)
          @info "if inside(particle2.pos, target) && inside_energybarrier(particle1.pos, target)"

          #Determine scattering angle from binary collision
          bca_result = calculate_binary_collision(particle1, particle2, collisiongeometry, options)

          # Only use 0th order collision for local electronic stopping
          if k == 0
            normalized_distance_of_closest_approach = bca_result.normalized_distance_of_closest_approach
            strong_collision_Z = particle2.Z
            strong_collision_index = species_index
          end

          # Energy transfer to recoil
          particle2.E = bca_result.E_recoil - average_property_atpos(:Eb, particle2.pos, target)
          #TODO por que el origen es ese?
          particle2.energy_origin = particle2.E

          # Accumulate energy losses and asymptotic deflections for primary particle
          total_energy_loss += bca_result.E_recoil
          total_asymptotic_deflection += bca_result.asymptotic_deflection

          rotate!(particle1, bca_result.ψ, collisiongeometry.ϕ_azimuthal)
          rotate!(particle2, -bca_result.ψ_recoil, collisiongeometry.ϕ_azimuthal)

          particle2.dir_old = particle2.dir

          #Only track number of strong collisions, i.e., k = 0
          if bca_result.ψ > 0.0 * unit(bca_result.ψ) && k == 0
            particle1.number_collision_events += 1
          end

          #Deep recoil suppression
          #See Eckstein 1991 7.5.3 for recoil suppression function
          if options.track_recoils && options.suppress_deep_recoils
            #TODO
          end
          # If transferred energy > cutoff energy, add recoil to particle vector
        elseif options.track_recoils && particle2.E > particle1.Ec
          push!(particles, particle2)
        end

      end
      #Advance particle in space and track total distance traveled
      if options.accelerated_ions
        #TODO whats up with the energy barrier thickness method?
        #TODO check math
        distance_traveled = advance!(particle1, collisiongeometries[1].mfp + distance_to_target - target.top_energy_barrier_thickness, total_asymptotic_deflection)
      else
        distance_traveled = advance!(particle1, collisiongeometries[1].mfp, total_asymptotic_deflection)
      end
      #Subtract total energy from all simultaneous collisions and electronic stopping
      update_particle_energy!(particle1, target, distance_traveled, total_energy_loss, normalized_distance_of_closest_approach, strong_collision_Z, strong_collision_index, options)

      #Check boundary conditions on leaving and stopping
      boundary_condition_planar!(particle1, target)

      #Set particle index to topmost particle
      particle_index = length(particles)
      @info "particle len after new index $(length(particles))"
    end
    push!(particle_output, particle1)
    @info "particle output len after push $(length(particle_output))"
  end
  return particle_output
end

function bca(particles::Vector{Particle}, target::Target, options::Options)
  l = length(particles)
  result = Particle[]
  for (i, p) in enumerate(particles)
    println("ion $i/$l")
    r = singleionbca(p, target, options)
    push!(result, r...)
  end
  result
end
