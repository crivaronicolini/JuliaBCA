module BCA
export Particle, Material, Layer, Disk, InteractionPotential, Vec3, Options
include("particle.jl")
include("target.jl")
include("interactions.jl")

include("main.jl")

include("output.jl")
include("plotrecipes.jl")
end # module BCA
