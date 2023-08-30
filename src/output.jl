using BCA
using DataFrames
using Unitful


function trajectory(p::Particle, unit=u"Å")
  [uconvert.(unit, (t.pos.x, t.pos.y)) for t in p.trajectory]
end

function trajectory(ps::Vector{Particle}, unit=u"Å")
  [[uconvert.(unit, (t.pos.x, t.pos.y)) for t in p.trajectory] for p in ps]
end

# function destinations(particles::Vector{Particle})
#   ma = vcat([[p.Z p.m p.E p.pos.x p.pos.y p.pos.z p.dir.x p.dir.y p.dir.z p.is_from_target] for p in particles]...)
#   E = ma[:, 3]
#   x = ma[:, 4]
#   is_from_target = Bool.(Int16.(ma[:, 10]))
#   sputtered = ma[is_from_target.&&E.>0u"eV", :]
#   reflected = ma[.!is_from_target.&&x.<0u"μm", :]
#   implanted = ma[.!is_from_target.&&x.>0u"μm", :]
#   return sputtered, reflected, implanted
# end

# function DataFrame(particles::Vector{Particle})
#   # matrix row for every particle, then vcat
#   df = DataFrame(vcat([[p.Z p.m p.E p.pos.x p.pos.y p.pos.z p.dir.x p.dir.y p.dir.z p.is_from_target] for p in particles]...), [:Z, :m, :E, :x, :y, :z, :ux, :uy, :uz, :is_from_target])
#   mapcols!(c -> uconvert.(unit(c[1]), c), df)
#   df.Z = Int16.(df.Z)
#   return df
# end

function destinations(df::DataFrame)
  sputtered = df[df.is_from_target.&&df.E.>0u"eV", :]
  reflected = df[.!df.is_from_target.&&df.x.<0u"μm", :]
  implanted = df[.!df.is_from_target.&&df.x.>0u"μm", :]
  return sputtered, reflected, implanted
end
