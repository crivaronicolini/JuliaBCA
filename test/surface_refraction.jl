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

@testset "surface refraction" begin
  m = 1amu
  Z = 1
  E = 10.0eV
  Ec = 1.0eV
  Es = 5.76eV
  pos = Vec3(0.0, 0.0, 0.0) .* um
  dir = Vec3(√(1 / 2), √(1 / 2), 0.0)
  cosx, cosy, cosz = dir
  particle1 = Particle(m=m, Z=Z, E=E, Ec=Ec, Es=Es, pos=pos, dir=dir)

  # Test particle entering material and gaining energy

  # Eckstein's formulation for particle entering surface
  cosx_new = sqrt((E * cosx^2 + Es) / (E + Es))
  sinx = √(1.0 - cosx^2)
  sinx_new = √(1.0 - cosx_new^2)
  cosy_new = cosy * sinx_new / sinx
  cosz_new = cosz * sinx_new / sinx
  dir_mag = sqrt(cosx_new * cosx_new + cosy_new * cosy_new + cosz_new * cosz_new)
  cosx_new = cosx_new / dir_mag
  cosy_new = cosy_new / dir_mag
  cosz_new = cosz_new / dir_mag

  dir_new = Vec3(cosx_new, cosy_new, cosz_new)

  normal = Vec3(1.0, 0.0, 0.0)
  BCA.surface_refraction!(particle1, normal, Es)
  @debug "surface refraction"

  @debug dir_mag dir dir_new particle1.dir

  @test particle1.dir ≈ dir_new atol = 1e-12

  # Test particle leaving material and losing energy

  cosx_new = sqrt((particle1.E * particle1.dir.x^2 - Es) / (particle1.E - Es))
  sinx = sqrt(1.0 - particle1.dir.x^2)
  sinx_new = sqrt(1.0 - cosx_new^2)
  cosy_new = particle1.dir.y * sinx_new / sinx
  cosz_new = particle1.dir.z * sinx_new / sinx
  dir_new = Vec3(cosx_new, cosy_new, cosz_new)

  @debug particle1.dir dir_new

  normal = Vec3(1.0, 0.0, 0.0)
  BCA.surface_refraction!(particle1, normal, -Es)

  @debug "surface refraction"
  @debug particle1.dir dir_new

  @test particle1.dir ≈ dir_new atol = 1e-12
end
