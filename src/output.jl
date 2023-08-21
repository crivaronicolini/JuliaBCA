using BCA
using Unitful


function trajectory(p::Particle, unit=u"Å")
  [uconvert.(unit, (t.pos.x, t.pos.y)) for t in p.trajectory]
end

function trajectory(ps::Vector{Particle}, unit=u"Å")
  [[uconvert.(unit, (t.pos.x, t.pos.y)) for t in p.trajectory] for p in ps]
end

# for t in ts
#   plot!(t)
# end
