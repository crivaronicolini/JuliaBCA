using Unitful
import Unitful: Length, Energy, Mass, Wavenumber, c0, me, Na
# import PhysicalConstants.CODATA2018: c_0, e, ε_0, a_0, m_e
const ang = u"Å"

#from PhysicalConstants, adding units
const e = 1.602176634e-19 * u"C"
const a_0 = 5.29177210903e-11 * u"m"
# from wikipedia, CODATA2018, Molar mass constant
const Mu = (0.99999999965 * 10^(-3))u"kg/mol"

const ε0 = uconvert(u"F/m", Unitful.ε0)

# global BETHE_BLOCH_PREFACTOR = (e^2 / (4π * ε0))^2 * 4π / (me * c0^2)
const BETHE_BLOCH_PREFACTOR = (4π / (me * c0^2)) * (e^2 / (4π * ε0))^2
# TODO CHANGED e BY eV
const LINDHARD_SCHARFF_PREFACTOR = 1.212 * u"eV^(1/2)*Å^2"
const LINDHARD_REDUCED_ENERGY_PREFACTOR = 4334488014869623000000000000 * unit(ε0 / e^2)


@enum ElectronicStoppingMode INTERPOLATED LOWENERGYLOCAL LOWENERGYNONLOCAL LOWENERGYEQUIPARTITION
@derived_dimension AtomicConcentration Unitful.𝐋^3
@derived_dimension AtomicDensity inv(Unitful.𝐋^3)

# Data layout
# Target([ Layer1([ Material1, Material2]), Layer2([Material3, Material4]) )

# material properties, Materials are input to Layers and Targets
@kwdef struct Material
  m::Vector{Mass}
  Z::Vector{Int} # atomic numbers
  Eb::Vector{Energy} # target species bulk binding energies
  Ec::Vector{Energy} # target species cutoff energies
  Es::Vector{Energy} # target species surface binding energies (planar)
  density::Vector{AtomicDensity}
  # concentrations::Vector{AtomicConcentration}
  concentrations::Vector{Real} = inv.(uconvert.(NoUnits, density * u"Å^3"))# concentrations siempre se usa normalizado por amu
  interactionindex::Vector{Int} = ones(Int64, length(m)) # index of interaction matrix to use for each species. Defaults to [1], which means all interactions use the same potential that occupies i=1, j=1 of the interaction index matrix
  electronic_stopping_correction_factor::Vector{Real} = ones(Float64, length(m))
  # TODO chequear en new que todos los vectores tengan el mismo largo
end


Material(m::Mass, Z::Int, Eb::Energy, Ec::Energy, Es::Energy, d::AtomicDensity) = Material(m=[m], Z=[Z], Eb=[Eb], Ec=[Ec], Es=[Es], density=[d])

Material(m::Mass, Z::Int, Eb::Energy, Ec::Energy, Es::Energy, d::AtomicDensity, interactionindex::Int, electronic_stopping_correction_factor::Real) = Material(m=[m], Z=[Z], Eb=[Eb], Ec=[Ec], Es=[Es], density=[d], interactionindex=interactionindex, electronic_stopping_correction_factor=electronic_stopping_correction_factor)

# Materials sum to create a composite material, like a metal alloy
Base.:+(x::Material, y::Material) = Material(m=[x.m; y.m], Z=[x.Z; y.Z], Eb=[x.Eb; y.Eb], Ec=[x.Ec; y.Ec], Es=[x.Es; y.Es], density=[x.density; y.density], concentrations=[x.concentrations; y.concentrations], interactionindex=[1, 1], electronic_stopping_correction_factor=[1.0, 1.0])

# function Base.show(io::Base.IO, m::Material)
#   println("$(length(m.m))-component Material")
#   for f in fieldnames(Material)
#     println(io, f, "= [", strfield(m, f), ']')
#   end
# end


# Material(c::ChemElem) = Material(c.atomic_mass, c.number, 1u"eV", 1u"eV", 1u"eV", 1u"amu/cm^3")
# Material(c::ChemElem, Eb::Energy, Ec::Energy, Es::Energy) = Material(c.atomic_mass, c.number, Eb, Ec, Es, 1u"amu/cm^3")

struct Layer
  material::Material
  thickness::Length
end

abstract type Target end
abstract type Target3D <: Target end
abstract type Target2D <: Target end
abstract type Target1D <: Target end
abstract type Target0D <: Target end

# Layer(m::Material, t::Length) = Layer([m], [t])

# Base.:+(x::Layer, y::Layer) = Layer([x.materials...; y.materials...], [x.thicknesess...; y.thicknesess...])

function strfield(x::Any, f::Any)
  join(getproperty(x, Symbol(f)), " ")
end

# function Base.show(io::Base.IO, l::Layer)
#   println("$(length(l.thicknesess))-material Layer")
#   println("materials =")
#   for m in l.materials
#     for f in fieldnames(Material)
#       println(io, "\t", f, "= [", strfield(m, f), ']')
#     end
#     println()
#   end
#   print("thicknesess= [", strfield(l, "thicknesess"), ']')
# end

struct Disk <: Target3D
  layers::Vector{Layer}
  diameter::Length
  interfaces::Vector{Length}
  top_energy_barrier_thicknesess::Length
  bottom_energy_barrier_thicknesess::Length
  function Disk(l::Layer, d::Length)
    Disk([l], d)
  end
  function Disk(ls::Vector{Layer}, d::Length)
    top_thick, bot_thick = energy_barrier_thicknesess(ls)
    interfaces = vcat([0.0u"μm"], accumulate(+, [l.thickness for l in ls]))
    new(ls, d, interfaces, top_thick, bot_thick)
  end
end

# Disk(l::Layer, d::Length) = Disk(layer=l, diameter=d)
Disk(m::Material, t::Length, d::Length) = Disk([Layer(m, t)], d)

Base.:+(x::Disk, y::Disk) = x.diameter == y.diameter ? Disk([x.layers..., y.layers...], x.diameter) : error("Diameters must be equal")

# function Base.show(io::Base.IO, d::Disk)
#   nlayers = (length(d.layer.thicknesess))
#   println("$nlayers-layers Disk")
#   # println("layers =")
#   for l in d.layer
#     @show l
#   end
#   #   for f in fieldnames(Layer)
#   #     println(io, "\t", f, "= [", strfield(l, f), ']')
#   #   end
#   #   println()
#   # end
#   print("Diameter = $d.diameter")
#   print("Top Energy Barrier Thickness = $d.top_energy_barrier_thicknesess")
#   print("Bottom Energy Barrier Thickness = $d.bottom_energy_barrier_thicknesess")
# end

function radius(x::Length, y::Length)
  √(x^2 + y^2)
end

#this function depends on the specific target geometry, so target must be specific
function inside(pos::Vec3, disk::Disk)
  x, y, z = pos
  if 0 * x < x < disk.interfaces[end] && radius(y, z) < disk.diameter / 2
    return true
  else
    return false
  end
end

#
function inside_simulation_boundary(pos::Vec3, disk::Disk)
  x, y, z = pos
  bot_thickness = -10.0 * disk.bottom_energy_barrier_thicknesess
  top_thickness = 10.0 * disk.top_energy_barrier_thicknesess

  if bot_thickness <= x <= disk.interfaces[end] + top_thickness && radius(y, z) <= disk.diameter / 2 + top_thickness
    return true
  else
    return false
  end
end

function energy_barrier_thicknesess(layers::Vector{Layer})
  #buscar cuenta en algún libro, rustBCA usa la dimension de la longitud del layer (ver geometry.rs:157)
  dtop, dbottom = sum(layers[1].material.density), sum(layers[end].material.density)
  (dtop, dbottom) .^ (-1 / 3) .* (2 / √π)
end

function inside_energybarrier(pos::Vec3, disk::Disk)
  x, y, z = pos
  top_e, bot_e = disk.top_energy_barrier_thicknesess, disk.bottom_energy_barrier_thicknesess
  if -top_e <= x <= disk.interfaces[end] + bot_e && radius(y, z) <= disk.diameter + max(top_e, bot_e)
    return true
  else
    return false
  end
end

function layer_atpos(pos::Vec3, target::Target)
  layerindex = searchsortedlast(target.interfaces, pos.x)
  if layerindex == 0
    layerindex += 1
  end
  if layerindex >= length(target.interfaces)
    layerindex -= 1
  end
  return target.layers[layerindex]
end

function material_atpos(pos::Vec3, target::Target)
  layer_atpos(pos, target).material
end

function average_property_atpos(property::Symbol, pos::Vec3, target::Target)
  material = material_atpos(pos, target)
  sum(material.concentrations .* getfield(material, property))
end

function property_atpos(property::Symbol, pos::Vec3, target::Target)
  material = material_atpos(pos, target)
  getfield(material, property)
end

function properties_atpos(properties::Vector{Symbol}, pos::Vec3, target::Target)
  material = material_atpos(pos, target)
  answer = []
  for property in properties
    push!(answer, getfield(material, property))
  end
  return answer
end

# Determines the local mean free path from the formula sum(n(x, y))^(-1/3)
function meanfreepath(pos::Vec3, target::Target)
  total_number_density(pos, target)^(-1 / 3)
end

function total_number_density(pos::Vec3, target::Target)
  sum(property_atpos(:density, pos, target))
end

function cumulative_concentration_atpos(pos::Vec3, target::Target)
  cumsum(property_atpos(:concentrations, pos, target))
end

function electronic_density_atpos(pos::Vec3, target::Target)
  ρ, Z, m = properties_atpos([:density, :Z, :m], pos, target)
  Ar = m * u"1/u"#TODO relative atomic mass
  return @. Na * Z * ρ / Ar * Mu
end

# Choose the parameters of a target atom as a concentration-weighted random draw from the species in the triangle that contains or is nearest to (x, y).
function choose(recoil::Vec3, target::Target)
  m = layer_atpos(recoil, target).material
  for (componentidx, cumulative_concentration) in enumerate(cumulative_concentration_atpos(recoil, target))
    u = unit(cumulative_concentration)
    r = rand()
    if r * u < cumulative_concentration
      return componentidx, m.Z[componentidx],
      m.m[componentidx],
      m.Ec[componentidx],
      m.Es[componentidx],
      m.interactionindex[componentidx]
    end
  end
  @error "Input error: method choose() operation failed to choose a valid species. Check densities."
end

# Choose the parameters of a target atom as a concentration-weighted random draw from the species in the triangle that contains or is nearest to (x, y).
function choose_deterministic(recoil::Vec3, target::Target)
  m = layer_atpos(recoil, target).material
  for (componentidx, cumulative_concentration) in enumerate(cumulative_concentration_atpos(recoil, target))
    return componentidx, m.Z[componentidx],
    m.m[componentidx],
    m.Ec[componentidx],
    m.Es[componentidx],
    m.interactionindex[componentidx]
  end
  @error "Input error: method choose() operation failed to choose a valid species. Check densities."
end

function electronic_stopping_cross_sections(particle::Particle, target::Target, electronic_stopping_mode::ElectronicStoppingMode)
  E, m, Za = particle.E, particle.m, particle.Z
  pos = particle.pos
  stopping_powers = []
  # ck, ns, Zbs = properties_atpos([:electronic_stopping_correction_factor, :n, :Z], pos, target)
  ck_vec, Zbs = properties_atpos([:electronic_stopping_correction_factor, :Z], pos, target)
  # TODO ISSUE? el original devuelve un float en vez de un vector, revisar
  ck = ck_vec[1]

  for Zb in Zbs
    β = √(1 - (1 + E / (m * c0^2))^(-2))
    v = β * c0

    # This term is an empirical fit to the mean ionization potential
    Iₒ = Zb < 13 ? 12.0 + 7.0 / Zb : 9.76 + 58.5Zb^(-1.19)
    I = Zb * Iₒ * 10u"eV"

    #See Biersack and Haggmark - this looks like an empirical shell correction
    C = Zb < 3 ? 100Za / Zb : 5.0
    #Bethe stopping modified by Biersack and Varelas
    prefactor = uconvert(u"eV*Å^2", BETHE_BLOCH_PREFACTOR * Zb * Za^2 / β^2)
    # ISSUE no dan las unidades en el original, Zb*I0 deberia dar unidades de energia, no lo veo.
    # Wikipedia tiene la formula con I=10eV*Z, uso esos 10eV (https://en.wikipedia.org/wiki/Bethe_formula)
    # el codigo original usa I = Zb * I0 * e
    eb = 2me * v^2 / I
    # TODO ISSUE? la formula es distinta en el paper, no usan el prefactor de bethe bloch, usan:
    # (el Zb está para cancelar el de I)
    # S_high = log(eb + 1 + C / eb) * 8π * Za^2 * e^4 * Zb / (I * eb)
    # ORIGINAL:
    S_high = prefactor * log(eb + 1 + C / eb)

    #Lindhard-Scharff electronic stopping
    # TODO averiguar que onda la masa dividiendo. Saco el / e que dividia a la energia, las unidades quedan bien sin eso, y poniendo el sqrt(amu) a mano, aparece en la formula dela paper
    # S_low = (ck * 1.212 * √(E / m) * Za^(7 / 6) * Zb) / (Za^(2 / 3) + Zb^(2 / 3))^(3 / 2)
    S_low = (ck * LINDHARD_SCHARFF_PREFACTOR * √(E / m) * Za^(7 / 6) * Zb) / (Za^(2 / 3) + Zb^(2 / 3))^(3 / 2) * u"sqrt(u)"
    # Las unidades de S deberian ser eV*Å^2

    # TODO
    stopping_power = electronic_stopping_mode == INTERPOLATED ? 1 / (1 / S_high + 1 / S_low) : S_low
    # electronic_stopping_mode == LOWENERGYLOCAL ? S_low :
    # electronic_stopping_mode == LOWENERGYNONLOCAL ? S_low :
    # electronic_stopping_mode == LOWENERGYEQUIPARTITION ? S_low :
    push!(stopping_powers, stopping_power)
  end
  return stopping_powers
end

function closest_point(pos::Vec3, disk::Disk)
  x, y, z = pos
  radialdistance = abs(radius(y, z) - disk.diameter / 2)
  topdistance = abs(x - disk.interfaces[end])
  botdistance = abs(x)

  if topdistance < radialdistance
    if topdistance < botdistance
      return Vec3(topdistance, y, z)
    else
      return Vec3(botdistance, y, z)
    end
  else
    if radialdistance < botdistance
      y′, z′ = radius(y, z) .* (cos(y), sin(z))
      return Vec3(x, y′, z′)
    else
      return Vec3(botdistance, y, z)
    end

  end
end

