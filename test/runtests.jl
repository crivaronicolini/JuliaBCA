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

@testset "all tests" begin
  include("deterministic_momentum.jl")
  include("momentum_conservation.jl")
  include("particle_advance.jl")
  include("rotate.jl")
  include("singleion.jl")
  include("surface_refraction.jl")
end

# @testset "surface refraction" begin
#   m = 1amu
#   Z = 1
#   E = 10.0eV
#   Ec = 1.0eV
#   Es = 5.76eV
#   pos = Vec3(0.0, 0.0, 0.0) .* um
#   dir = Vec3(√(1 / 2), √(1 / 2), 0.0)
#   cosx, cosy, cosz = dir
#   particle1 = Particle(m=m, Z=Z, E=E, Ec=Ec, Es=Es, pos=pos, dir=dir)
#
#   # Test particle entering material and gaining energy
#
#   # Eckstein's formulation for particle entering surface
#   cosx_new = sqrt((E * cosx^2 + Es) / (E + Es))
#   sinx = √(1.0 - cosx^2)
#   sinx_new = √(1.0 - cosx_new^2)
#   cosy_new = cosy * sinx_new / sinx
#   cosz_new = cosz * sinx_new / sinx
#   dir_mag = sqrt(cosx_new * cosx_new + cosy_new * cosy_new + cosz_new * cosz_new)
#   cosx_new = cosx_new / dir_mag
#   cosy_new = cosy_new / dir_mag
#   cosz_new = cosz_new / dir_mag
#
#   dir_new = Vec3(cosx_new, cosy_new, cosz_new)
#
#   normal = Vec3(1.0, 0.0, 0.0)
#   BCA.surface_refraction!(particle1, normal, Es)
#   @debug "surface refraction"
#
#   @debug dir_mag dir dir_new particle1.dir
#
#   @test particle1.dir ≈ dir_new atol = 1e-12
#
#   # Test particle leaving material and losing energy
#
#   cosx_new = sqrt((particle1.E * particle1.dir.x^2 - Es) / (particle1.E - Es))
#   sinx = sqrt(1.0 - particle1.dir.x^2)
#   sinx_new = sqrt(1.0 - cosx_new^2)
#   cosy_new = particle1.dir.y * sinx_new / sinx
#   cosz_new = particle1.dir.z * sinx_new / sinx
#   dir_new = Vec3(cosx_new, cosy_new, cosz_new)
#
#   @debug particle1.dir dir_new
#
#   normal = Vec3(1.0, 0.0, 0.0)
#   BCA.surface_refraction!(particle1, normal, -Es)
#
#   @debug "surface refraction"
#   @debug particle1.dir dir_new
#
#   @test particle1.dir ≈ dir_new atol = 1e-12
# end
#
# @testset "momentum conservation" begin
#   for energy_eV in [1]eV
#     # for energy_eV in [1, 10, 100, 1000, 10000]eV
#     m1 = 183.8amu
#     Z1 = 74
#     E1 = energy_eV
#     Ec1 = 1eV
#     Es1 = 1eV
#     pos = Vec3(0.0, 0.0, 0.0) .* um
#
#     material = Material(6.941amu, 3, 0eV, 1eV, 1eV, 0.06306u"1/Å^3")
#     target = Disk(material, 0.1um, 0.1um)
#
#     θ = 0.974194583091052rad
#     dir = Vec3(cos(θ), sin(θ), 0.0)
#
#     for high_energy_free_flight_paths in [false]
#       # for high_energy_free_flight_paths in [true, false]
#       for potential in [BCA.Moliere]
#         for scattering_int in [BCA.Mendenhall_weller]
#
#           @debug "Case:" energy_eV high_energy_free_flight_paths potential scattering_int
#           particle1 = Particle(m=m1, Z=Z1, E=energy_eV, Ec=Ec1, Es=Es1, pos=pos, dir=dir)
#
#           options = Options(
#             name="test",
#             track_recoils=true,
#             weak_collision_order=1,
#             high_energy_free_flight_paths=high_energy_free_flight_paths,
#             electronic_stopping_mode=BCA.INTERPOLATED,
#             interaction_potential=[[potential]],
#             scattering_integral=scattering_int
#           )
#
#           binary_collision_geometries = BCA.determine_mfp_phi_impact_parameter(particle1, target, options)
#           @debug "ϕ_azimuthal " binary_collision_geometries[1].ϕ_azimuthal binary_collision_geometries[1].impactparameter binary_collision_geometries[1].mfp
#
#           species_index, particle2 = BCA.choose_collision_partner(particle1, target, binary_collision_geometries[1])
#
#           mom1_0 = BCA.get_momentum(particle1)
#           mom2_0 = BCA.get_momentum(particle2)
#
#           initial_momentum = uconvert.(u"u*Å/s", mom1_0 + mom2_0)
#
#           bca_result = BCA.calculate_binary_collision(particle1, particle2, binary_collision_geometries[1], options)
#
#           @debug "bca result" bca_result.E_recoil bca_result.ψ bca_result.ψ_recoil
#
#           @debug "Initial energies" particle1.E particle2.E
#
#           # Energy transfer to recoil
#           particle2.E = bca_result.E_recoil - BCA.average_property_atpos(:Eb, particle2.pos, target)
#           @debug "after recoil energies" particle2.E
#
#           # Rotate particle 1, 2 by lab frame scattering angles
#           BCA.rotate!(particle1, bca_result.ψ, binary_collision_geometries[1].ϕ_azimuthal)
#           BCA.rotate!(particle2, -bca_result.ψ_recoil, binary_collision_geometries[1].ϕ_azimuthal)
#
#           # Subtract total energy from all simultaneous collisions and electronic stopping
#           BCA.update_particle_energy!(particle1, target, 0um, bca_result.E_recoil, 0.0, particle2.Z, species_index, options)
#
#           @debug "Final energies" particle1.E particle2.E
#
#           mom1_1 = BCA.get_momentum(particle1)
#           mom2_1 = BCA.get_momentum(particle2)
#
#           final_momentum = uconvert.(u"u*Å/s", mom1_1 + mom2_1)
#           error = 100(final_momentum - initial_momentum) / BCA.norm(initial_momentum)
#
#           @debug "momentum error" initial_momentum final_momentum error
#
#           @test initial_momentum.x.val ≈ final_momentum.x.val atol = 1e-12
#           @test initial_momentum.y.val ≈ final_momentum.y.val atol = 1e-12
#           @test initial_momentum.z.val ≈ final_momentum.z.val atol = 1e-12
#
#           @test !isnan(particle1.E)
#           @test !isnan(particle2.E)
#           @test !isnan(initial_momentum.x)
#           @test !isnan(initial_momentum.y)
#           @test !isnan(initial_momentum.z)
#         end
#       end
#     end
#   end
#
# end
#
# @testset "rotate" begin
#   pos = Vec3(0.0, 0.0, 0.0) .* um
#   dir = Vec3(cos(π / 4), sin(π / 4), 0.0)
#   ψ = (-π / 4)rad
#   ϕ = 0rad
#
#   particle = Particle(m=1amu, Z=1, E=1eV, Ec=1eV, Es=1eV, pos=pos, dir=dir)
#
#   # Check that rotation in 2D works
#   BCA.rotate!(particle, ψ, ϕ)
#   @test particle.dir.x ≈ 0 atol = 1e-12
#   @test particle.dir.y ≈ 1 atol = 1e-12
#
#   # Check that rotating back by negative psi returns to the previous values
#   BCA.rotate!(particle, -ψ, ϕ)
#   @test particle.dir.x ≈ dir.x atol = 1e-12
#   @test particle.dir.y ≈ dir.y atol = 1e-12
#
#   # Check that azimuthal rotation by 180 degrees works correctly
#   BCA.rotate!(particle, ψ, π * rad)
#   @test particle.dir.x ≈ 1 atol = 1e-12
#   @test particle.dir.y ≈ 0 atol = 1e-12
#
#   # Check that particle direction vector remains normalized following rotations
#   @test BCA.norm(particle.dir) ≈ 1
# end
#
# @testset "particle advance" begin
#   mfp = 1.0um
#   asymptotic_deflection = 0.5um
#   pos = Vec3(0.0, 0.0, 0.0) .* um
#   dir = Vec3(cos(π / 4), sin(π / 4), 0.0)
#
#   particle = Particle(m=1amu, Z=1, E=1eV, Ec=1eV, Es=1eV, pos=pos, dir=dir)
#
#   distance_traveled = BCA.advance!(particle, mfp, asymptotic_deflection)
#
#   @test particle.pos.x == um * dir.x / 2
#   @test particle.pos.y == um * dir.y / 2
#   @test particle.pos.z == um * 0.0
#   @test distance_traveled == mfp - asymptotic_deflection
# end
#
# @testset "test deterministic momentum" begin
#   for energy_eV in [1]eV
#     # for energy_eV in [1, 10, 100, 1000, 10000]eV
#     m1 = 183.8amu
#     Z1 = 74
#     E1 = energy_eV
#     Ec1 = 1eV
#     Es1 = 1eV
#     pos = Vec3(0.0, 0.0, 0.0) .* um
#
#     material = Material(6.941amu, 3, 0eV, 1eV, 1eV, 0.06306u"1/Å^3")
#     target = Disk(material, 0.1um, 0.1um)
#
#     θ = 0.974194583091052rad
#     dir = Vec3(cos(θ), sin(θ), 0.0)
#
#     for high_energy_free_flight_paths in [false]
#       # for high_energy_free_flight_paths in [true, false]
#       for potential in [BCA.Moliere]
#         for scattering_int in [BCA.Mendenhall_weller]
#
#           @debug "Case:" energy_eV high_energy_free_flight_paths potential scattering_int
#           particle1 = Particle(m=m1, Z=Z1, E=energy_eV, Ec=Ec1, Es=Es1, pos=pos, dir=dir)
#
#           options = Options(
#             name="test",
#             track_recoils=true,
#             weak_collision_order=1,
#             high_energy_free_flight_paths=high_energy_free_flight_paths,
#             electronic_stopping_mode=BCA.INTERPOLATED,
#             interaction_potential=[[potential]],
#             scattering_integral=scattering_int
#           )
#
#           binary_collision_geometries = BCA.determine_mfp_phi_impact_parameter_deterministic(particle1, target, options)
#
#           @debug "Phi: $(binary_collision_geometries[1].ϕ_azimuthal), p: $(binary_collision_geometries[1].impactparameter), mfp: $(binary_collision_geometries[1].mfp)"
#
#           @test binary_collision_geometries[1].ϕ_azimuthal ≈ 1.2566370614359172rad
#           @test binary_collision_geometries[1].impactparameter ≈ 0.6339019280131679ang
#           @test binary_collision_geometries[1].mfp ≈ 2.5123608152996693ang
#
#           species_index, particle2 = BCA.choose_collision_partner_deterministic(particle1, target, binary_collision_geometries[1])
#
#           @debug particle2
#
#           @test particle2.m ≈ 6.941amu
#           @test particle2.E ≈ 0eV
#           @test particle2.Ec ≈ 1eV
#           @test particle2.Es ≈ 1eV
#           @test uconvert.(ang, particle2.pos) ≈ Vec3(1.249484012440868, 2.1884053732017623, 0.6028765593289847)ang #recoil
#           @test particle2.dir ≈ Vec3(0.561834516559727, 0.8272496455134314, 0.0) #cosdir
#           @test particle2.asymptotic_deflection ≈ 0um
#
#           mom1_0 = BCA.get_momentum(particle1)
#           mom2_0 = BCA.get_momentum(particle2)
#
#           initial_momentum = uconvert.(u"u*Å/s", mom1_0 + mom2_0)
#
#           bca_result = BCA.calculate_binary_collision(particle1, particle2, binary_collision_geometries[1], options)
#
#           @debug "bca result" bca_result.E_recoil bca_result.ψ bca_result.ψ_recoil
#
#           @debug "Initial energies" particle1.E particle2.E
#
#           # Case: 1 false Moliere Potential Mendenhall-Weller 4-Point Lobatto Quadrature
#           #
#           # Phi: 1.2566370614359172 rad p: 0.6339019280131679 Angstrom mfp: 2.5123608152996693 Angstrom
#           #
#           # Particle 2: m: 6.941, E: 0, Ec: 1, Es: 1, pos: 1.249484012440868, 2.1884053732017623, 0.6028765593289847, dir: 0.561834516559727, 0.8272496455134314, 0, asy_defl: 0
#           #
#           # E_recoil: 0.13261896212622665 eV Psi: 0.017738493230885534 rad Psi_recoil: 0.23560429073267436 rad
#           #
#           # Initial Energies: 1 eV 0 eV
#
#         end
#       end
#     end
#   end
# end
#
#
# @testset "normal run" begin
#   boron = @test_nowarn Material(10.811amu, 5, 0eV, 1eV, 5.76eV, 0.065u"1/Å^3")
#   nitride = Material(14amu, 7, 0eV, 1eV, 0eV, 0.065u"1/Å^3")
#
#   target = @test_nowarn Disk(boron + nitride, 100um, 2u"inch")
#
#   opts = Options(name="hydrogen on boron-nitride", track_trajectories=true)
#
#   p = @test_nowarn BCA.default_incident(1.008amu, 1, 1000eV, 1eV, 10eV, 0um, [0.9999999999984769, 1.7453292519934434e-6, 0.0])
#   p.track_trajectories = true
#   result1 = BCA.singleionbca(p, target, opts)
# end
