using RecipesBase
using Unitful
using BCA

@recipe function f(p::Particle, unit=u"Å")
  BCA.trajectory(p, unit)
end

@userplot TrajectoryPlot
@recipe function f(t::TrajectoryPlot)
  ps = t.args[1]
  unit = length(t.args) == 2 ? t.args[2] : u"Å"
  # p = plot()
  # map(p -> plot!(p, unit), ps)
  # display(p)
  framestyle := :grid
  aspect_ratio --> true
  title --> "Ion trajectories"
  xlabel --> "x "
  ylabel --> "y "

  @series begin
    seriestype := :vline
    label := nothing
    color := :gray
    [0]
  end

  # @series begin
  #   seriestype := :shape
  #   color := :gray
  #   [(0, 100), (0, -100), (100, -100), (100, 100)]
  # end

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
      BCA.trajectory(p, unit)
    end
  end
end


