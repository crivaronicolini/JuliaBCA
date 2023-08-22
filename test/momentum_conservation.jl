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

@testset "momentum conservation" begin
  for energy_eV in [1.0]eV
    # for energy_eV in [1, 10, 100, 1000, 10000]eV
    m1 = 183.8amu
    Z1 = 74
    E1 = energy_eV
    Ec1 = 1.0eV
    Es1 = 1.0eV
    pos = Vec3(0.0, 0.0, 0.0) .* um

    material = Material(6.941amu, 3, 0.0eV, 1.0eV, 1.0eV, 0.06306u"1/Å^3")
    target = Disk(material, 0.1um, 0.1um)

    θ = 0.974194583091052rad
    dir = Vec3(cos(θ), sin(θ), 0.0)

    for high_energy_free_flight_paths in [false]
      # for high_energy_free_flight_paths in [true, false]
      for potential in [BCA.Moliere]
        for scattering_int in [BCA.Mendenhall_weller]

          # @debug "Case:" energy_eV high_energy_free_flight_paths potential scattering_int
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

          binary_collision_geometries = BCA.determine_mfp_phi_impact_parameter(particle1, target, options)
          # @debug "ϕ_azimuthal " binary_collision_geometries[1].ϕ_azimuthal binary_collision_geometries[1].impactparameter binary_collision_geometries[1].mfp

          species_index, particle2 = BCA.choose_collision_partner(particle1, target, binary_collision_geometries[1], options)

          mom1_0 = BCA.get_momentum(particle1)
          mom2_0 = BCA.get_momentum(particle2)

          initial_momentum = uconvert.(u"u*Å/s", mom1_0 + mom2_0)

          bca_result = BCA.calculate_binary_collision(particle1, particle2, binary_collision_geometries[1], options)

          # @debug "bca result" bca_result.E_recoil bca_result.ψ bca_result.ψ_recoil

          # @debug "Initial energies" particle1.E particle2.E

          # Energy transfer to recoil
          particle2.E = bca_result.E_recoil - BCA.average_property_atpos(:Eb, particle2.pos, target)
          # @debug "after recoil energies" particle2.E

          # Rotate particle 1, 2 by lab frame scattering angles
          BCA.rotate!(particle1, bca_result.ψ, binary_collision_geometries[1].ϕ_azimuthal)
          BCA.rotate!(particle2, -bca_result.ψ_recoil, binary_collision_geometries[1].ϕ_azimuthal)

          # Subtract total energy from all simultaneous collisions and electronic stopping
          BCA.update_particle_energy!(particle1, target, 0.0um, bca_result.E_recoil, 0.0, particle2.Z, species_index, options)

          # @debug "Final energies" particle1.E particle2.E

          mom1_1 = BCA.get_momentum(particle1)
          mom2_1 = BCA.get_momentum(particle2)

          final_momentum = uconvert.(u"u*Å/s", mom1_1 + mom2_1)
          error = 100(final_momentum - initial_momentum) / BCA.norm(initial_momentum)

          # @debug "momentum error" initial_momentum final_momentum error

          @test initial_momentum.x.val ≈ final_momentum.x.val atol = 1e-12
          @test initial_momentum.y.val ≈ final_momentum.y.val atol = 1e-12
          @test initial_momentum.z.val ≈ final_momentum.z.val atol = 1e-12

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
