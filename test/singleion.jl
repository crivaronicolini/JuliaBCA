using Test
using Plots
using BCA
using Unitful
using DimensionfulAngles
using BenchmarkTools
amu = u"u"
eV = u"eV"
um = u"μm"
cm = u"cm"
rad = u"radᵃ"
ang = u"Å"

@testset "singleion normal run" begin
  boron = @test_nowarn Material(10.811amu, 5, 0eV, 1eV, 5.76eV, 0.065u"1/Å^3")
  nitride = Material(14amu, 7, 0eV, 1eV, 0eV, 0.065u"1/Å^3")
  target = @test_nowarn Disk(boron + nitride, 100um, 2u"inch")
  opts = Options(name="hydrogen on boron-nitride", track_trajectories=true)

  @benchmark begin
    p = BCA.default_incident(1.008amu, 1, 1000eV, 1eV, 10eV, 0um, [0.9999999999984769, 1.7453292519934434e-6, 0.0], track_trajectories=true)
    result1 = BCA.singleionbca(p, $target, $opts)
  end

end


@testset "multiple ion normal run" begin

  boron = @test_nowarn Material(10.811amu, 5, 0eV, 1eV, 5.76eV, 0.065u"1/Å^3")
  nitride = Material(14amu, 7, 0eV, 1eV, 0eV, 0.065u"1/Å^3")
  target = @test_nowarn Disk(boron + nitride, 100ang, 2u"inch")
  opts = Options(name="hydrogen on boron-nitride", track_trajectories=true, track_recoil_trajectories=true)
  H = @test_nowarn [BCA.default_incident(1.008amu, 1, 1000eV, 1eV, 10eV, 0um, [0.9999999999984769, 1.7453292519934434e-6, 0.0], track_trajectories=true) for _ in 1:20]
  B = [BCA.default_incident(1.008amu, 5, 1000eV, 1eV, 10eV, 0um, [0.9999999999984769, 1.7453292519934434e-6, 0.0], track_trajectories=true) for _ in 1:20]
  Zn = [BCA.default_incident(1.008amu, 30, 9000eV, 1eV, 10eV, 0um, [0.9999999999984769, 1.7453292519934434e-6, 0.0], track_trajectories=true) for _ in 1:20]
  ps = vcat(H, B, Zn)
  result = BCA.bca(ps, target, opts)
  trajectoryplot(result)

end

# using Logging
# debuglogger = ConsoleLogger(stderr, Logging.Debug)
# global_logger(debuglogger)

# infologger = ConsoleLogger(stderr, Logging.Info)
# global_logger(infologger)
