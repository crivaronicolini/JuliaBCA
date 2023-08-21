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

@testset "particle advance" begin
  mfp = 1.0um
  asymptotic_deflection = 0.5um
  pos = Vec3(0.0, 0.0, 0.0) .* um
  dir = Vec3(cos(π / 4), sin(π / 4), 0.0)

  particle = Particle(m=1amu, Z=1, E=1eV, Ec=1eV, Es=1eV, pos=pos, dir=dir)

  distance_traveled = BCA.advance!(particle, mfp, asymptotic_deflection)

  @test particle.pos.x == um * dir.x / 2
  @test particle.pos.y == um * dir.y / 2
  @test particle.pos.z == um * 0.0
  @test distance_traveled == mfp - asymptotic_deflection
end
