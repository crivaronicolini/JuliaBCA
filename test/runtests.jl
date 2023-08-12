using Test
using BCA
using Unitful
using DimensionfulAngles

amu = u"u"
eV = u"eV"
um = u"μm"
cm = u"cm"
rad = u"radᵃ"

@testset "momentum conservation" begin
  for energy_eV in [1]eV
    # for energy_eV in [1, 10, 100, 1000, 10000]eV
    m1 = 183.8amu
    Z1 = 74
    E1 = energy_eV
    Ec1 = 1eV
    Es1 = 1eV
    pos = Vec3(0.0, 0.0, 0.0) .* um

    material = Material(6.941amu, 3, 0eV, 1eV, 1eV, 0.06306u"1/Å^3")
    target = Disk(material, 0.1um, 0.1um)

    θ = 0.974194583091052rad
    dir = Vec3(cos(θ), sin(θ), 0.0)

    for high_energy_free_flight_paths in [false]
      # for high_energy_free_flight_paths in [true, false]
      for potential in [BCA.Moliere]
        for scattering_int in [BCA.Mendenhall_weller]

          @info "Case:" energy_eV high_energy_free_flight_paths potential scattering_int
          particle1 = Particle(m=m1, Z=Z1, E=energy_eV, Ec=Ec1, Es=Es1, pos=pos, dir=dir)

          options = Options(
            name="test",
            track_recoils=true,
            weak_collision_order=1,
            high_energy_free_flight_paths=high_energy_free_flight_paths,
            electronic_stopping_mode=BCA.INTERPOLATED,
            interaction_potential=[[potential]],
            scattering_integral=scattering_int
          )
          binary_collision_geometries = BCA.choosecollisiongeometries(particle1, target, options)
          @info "ϕ_azimuthal " binary_collision_geometries[1].ϕ_azimuthal binary_collision_geometries[1].impactparameter binary_collision_geometries[1].mfp

          species_index, particle2 = BCA.choose_collision_partner(particle1, target, binary_collision_geometries[1])

          mom1_0 = BCA.get_momentum(particle1)
          mom2_0 = BCA.get_momentum(particle2)

          initial_momentum = uconvert.(u"u*Å/s", mom1_0 + mom2_0)

          bca_result = BCA.calculate_binary_collision(particle1, particle2, binary_collision_geometries[1], options)

          @info "bca result" bca_result.E_recoil bca_result.ψ bca_result.ψ_recoil

          @info "Initial energies" particle1.E particle2.E

          # Energy transfer to recoil
          particle2.E = bca_result.E_recoil - BCA.average_property_atpos(:Eb, particle2.pos, target)
          @info "after recoil energies" particle2.E

          # Rotate particle 1, 2 by lab frame scattering angles
          BCA.rotate!(particle1, bca_result.ψ, binary_collision_geometries[1].ϕ_azimuthal)
          BCA.rotate!(particle2, -bca_result.ψ_recoil, binary_collision_geometries[1].ϕ_azimuthal)

          # Subtract total energy from all simultaneous collisions and electronic stopping
          BCA.update_particle_energy!(particle1, target, 0um, bca_result.E_recoil, 0.0, particle2.Z, species_index, options)

          @info "Final energies" particle1.E particle2.E

          mom1_1 = BCA.get_momentum(particle1)
          mom2_1 = BCA.get_momentum(particle2)

          final_momentum = uconvert.(u"u*Å/s", mom1_1 + mom2_1)
          error = 100(final_momentum - initial_momentum) ./ final_momentum

          @info "momentum error" initial_momentum final_momentum error

          @test_broken initial_momentum.x.val ≈ final_momentum.x.val atol = 1e-12
          @test_broken initial_momentum.y.val ≈ final_momentum.y.val atol = 1e-12
          @test_broken initial_momentum.z.val ≈ final_momentum.z.val atol = 1e-12

          @test !isnan(particle1.E)
          @test !isnan(particle2.E)
          @test !isnan(initial_momentum.x)
          @test !isnan(initial_momentum.y)
          @test !isnan(initial_momentum.z)
        end
      end
    end
  end

end

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
