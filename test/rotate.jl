using Test
using BCA
using Unitful
using DimensionfulAngles

amu = u"u"
eV = u"eV"
um = u"μm"
cm = u"cm"
rad = u"radᵃ"
ang = u"Å"

@testset "rotate" begin
  pos = Vec3(0.0, 0.0, 0.0) .* um
  dir = Vec3(cos(π / 4), sin(π / 4), 0.0)
  ψ = (-π / 4)rad
  ϕ = 0rad

  particle = Particle(m=1amu, Z=1, E=1eV, Ec=1eV, Es=1eV, pos=pos, dir=dir)

  # Check that rotation in 2D works
  BCA.rotate!(particle, ψ, ϕ)
  @test particle.dir.x ≈ 0 atol = 1e-12
  @test particle.dir.y ≈ 1 atol = 1e-12

  # Check that rotating back by negative psi returns to the previous values
  BCA.rotate!(particle, -ψ, ϕ)
  @test particle.dir.x ≈ dir.x atol = 1e-12
  @test particle.dir.y ≈ dir.y atol = 1e-12

  # Check that azimuthal rotation by 180 degrees works correctly
  BCA.rotate!(particle, ψ, π * rad)
  @test particle.dir.x ≈ 1 atol = 1e-12
  @test particle.dir.y ≈ 0 atol = 1e-12

  # Check that particle direction vector remains normalized following rotations
  @test BCA.norm(particle.dir) ≈ 1
end
