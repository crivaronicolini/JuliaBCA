using Test
using BCA
using Unitful
using DimensionfulAngles

const amu = u"u"
const eV = u"eV"
const um = u"μm"
const cm = u"cm"
const rad = u"radᵃ"
const ang = u"Å"

@testset "test deterministic momentum" begin
  for energy_eV in [1]eV
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

          binary_collision_geometries = BCA.determine_mfp_phi_impact_parameter_deterministic(particle1, target, options)

          # @debug "Phi: $(binary_collision_geometries[1].ϕ_azimuthal), p: $(binary_collision_geometries[1].impactparameter), mfp: $(binary_collision_geometries[1].mfp)"

          @test binary_collision_geometries[1].ϕ_azimuthal ≈ 1.2566370614359172rad
          @test binary_collision_geometries[1].impactparameter ≈ 0.6339019280131679ang
          @test binary_collision_geometries[1].mfp ≈ 2.5123608152996693ang

          species_index, particle2 = BCA.choose_collision_partner_deterministic(particle1, target, binary_collision_geometries[1], options)

          # @debug particle2

          @test particle2.m ≈ 6.941amu
          @test particle2.E ≈ 0eV
          @test particle2.Ec ≈ 1eV
          @test particle2.Es ≈ 1eV
          @test uconvert.(ang, particle2.pos) ≈ Vec3(1.249484012440868, 2.1884053732017623, 0.6028765593289847)ang #recoil
          @test particle2.dir ≈ Vec3(0.561834516559727, 0.8272496455134314, 0.0) #cosdir
          @test particle2.asymptotic_deflection ≈ 0um

          mom1_0 = BCA.get_momentum(particle1)
          mom2_0 = BCA.get_momentum(particle2)

          initial_momentum = uconvert.(u"u*Å/s", mom1_0 + mom2_0)
          # @debug initial_momentum
          @test initial_momentum.x ≈ 1058100248191806ang * amu * u"1/s"
          @test initial_momentum.y ≈ 1557955286539063ang * amu * u"1/s"
          @test initial_momentum.z ≈ 0ang * amu * u"1/s"

          # @debug "Initial energies" particle1.E particle2.E
          @test particle1.E ≈ 1eV
          @test particle2.E ≈ 0eV

          bca_result = BCA.calculate_binary_collision(particle1, particle2, binary_collision_geometries[1], options)

          # cosa = uconvert(u"mm*eV", a * relative_energy / (Za * Zb))
          @debug "cosas" reduced_energy relative_energy β LINDHARD_REDUCED_ENERGY_PREFACTOR cosa
          # ┌ Info: cosas
          # │   reduced_energy = 7.015601472104279e12 F eV C^-2
          # │   relative_energy = 0.03638965927619127 eV
          # │   β = 6.419777603650638
          # │   LINDHARD_REDUCED_ENERGY_PREFACTOR = 4334488014869623000000000000 F C^-2 m^-1
          # └   cosa = 1.6185536672467421e-12 mm eV
          # a: 0.000000000009874203861091632
          # reduced_energy: 0.000001124023275206148
          # relative_energy: 0.036389659276191276
          # beta: 6.419777603650646
          # LINDHARD_PREFACTOR: 4334488014869623000000000000
          # reduced_energy_sin prefactor: 0.0000000000016185536672467421
          @debug a = 9.874203861091632e-12 m
          # a: 0.000000000009874203861091632

          # @debug "bca result" bca_result.E_recoil bca_result.ψ bca_result.ψ_recoil bca_result.θ bca_result.asymptotic_deflection bca_result.normalized_distance_of_closest_approach
          @test bca_result.E_recoil ≈ 0.13261896212622665eV atol = 1e-9eV
          @test bca_result.ψ ≈ 0.017738493230885534rad atol = 1e-9rad
          @test bca_result.ψ_recoil ≈ 0.23560429073267436rad atol = 1e-9rad
          @test bca_result.θ ≈ 2.6703840721244445rad atol = 1e-9rad
          @test bca_result.asymptotic_deflection ≈ 0.00000000029645864477063087u"m" atol = 1e-9u"m"
          @test bca_result.normalized_distance_of_closest_approach ≈ 30.876561734011396 atol = 1e-6

          # Energy transfer to recoil
          particle2.E = bca_result.E_recoil - BCA.average_property_atpos(:Eb, particle2.pos, target)
          # @debug "after recoil energies" particle2.E
          @test particle2.E ≈ 0.13261896212622665eV atol = 1e-10eV

          # Rotate particle 1, 2 by lab frame scattering angles
          BCA.rotate!(particle1, bca_result.ψ, binary_collision_geometries[1].ϕ_azimuthal)
          BCA.rotate!(particle2, -bca_result.ψ_recoil, binary_collision_geometries[1].ϕ_azimuthal)
          # @debug "rotated directions" particle1.dir particle2.dir
          @test particle1.dir.x ≈ 0.566280454808479
          @test particle1.dir.y ≈ 0.8240399680879582
          @test particle1.dir.z ≈ -0.016869424871612276
          @test particle2.dir.x ≈ 0.48664012151167524
          @test particle2.dir.y ≈ 0.844922987033493
          @test particle2.dir.z ≈ 0.2220057164072626

          # Subtract total energy from all simultaneous collisions and electronic stopping
          BCA.update_particle_energy!(particle1, target, 0.0um, bca_result.E_recoil, 0.0, particle2.Z, species_index, options)

          # @debug "Final energies" particle1.E particle2.E
          @test particle1.E ≈ 0.8673810378737734eV
          @test particle2.E ≈ 0.13261896212622665eV

          mom1_1 = BCA.get_momentum(particle1)
          mom2_1 = BCA.get_momentum(particle2)

          final_momentum = uconvert.(u"u*Å/s", mom1_1 + mom2_1)
          error = 100(final_momentum - initial_momentum) / BCA.norm(initial_momentum)

          # @debug "momentum error" initial_momentum final_momentum error

          @test initial_momentum.x.val ≈ final_momentum.x.val
          @test initial_momentum.y.val ≈ final_momentum.y.val
          @test initial_momentum.z.val ≈ final_momentum.z.val atol = 1e-2

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
