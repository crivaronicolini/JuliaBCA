using RecipesBase
using Unitful
using BCA

@recipe function f(p::Particle, unit=u"Å")
  BCA.trajectory(p, unit)
end

@userplot TrajectoryPlot
@recipe function f(t::TrajectoryPlot)
  ps = t.args[1]
  target = t.args[2]
  unit = length(t.args) == 3 ? t.args[3] : u"Å"
  framestyle := :grid
  aspect_ratio --> true
  title --> "Ion trajectories"
  xlabel --> "x "
  ylabel --> "y "
  intf =target.interfaces
  for i in intf
    @series begin
      seriestype := :vline
      label := nothing
      if i==0.0u"Å" || i==intf[end]
        color := :gray
      else
        color := :lightgray
      end
      [i]
    end
  end
  Zs = Int64[]
  for p in ps
    @series begin
      seriestype := :path
      label := if p.Z ∉ Zs
        push!(Zs, p.Z)
        string(p.Z)
      else
        nothing
      end
      color := p.Z
      # z_order := p.Z
      BCA.trajectory(p, unit)
    end
  end
end


