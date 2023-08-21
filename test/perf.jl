using Test
using Plots
using BCA
using Unitful
using DimensionfulAngles
using BenchmarkTools
const amu = u"u"
const eV = u"eV"
const um = u"μm"
const cm = u"cm"
const rad = u"radᵃ"
const ang = u"Å"

@testset "perf tracking" begin

  boron = Material(10.811amu, 5, 0eV, 1eV, 5.76eV, 0.065u"1/Å^3")
  nitride = Material(14amu, 7, 0eV, 1eV, 0eV, 0.065u"1/Å^3")
  target = Disk(boron + nitride, 100um, 2u"inch")
  opts = Options(name="hydrogen on boron-nitride", track_trajectories=true)
  p = BCA.default_incident(1.008amu, 1, 1000eV, 1eV, 10eV, 0um, [0.9999999999984769, 1.7453292519934434e-6, 0.0], track_trajectories=true)
  result1 = BCA.singleionbca(p, $target, $opts)

  boron = Material(10.811amu, 5, 0eV, 1eV, 5.76eV, 0.065u"1/Å^3")
  nitride = Material(14amu, 7, 0eV, 1eV, 0eV, 0.065u"1/Å^3")
  target = Disk(boron + nitride, 100um, 2u"inch")
  opts = Options(name="hydrogen on boron-nitride", track_trajectories=true)
  p = BCA.default_incident(1.008amu, 1, 1000eV, 1eV, 10eV, 0um, [0.9999999999984769, 1.7453292519934434e-6, 0.0], track_trajectories=true)
  @report_opt BCA.singleionbca(p, target, opts)

  boron = Material(10.811amu, 5, 0eV, 1eV, 5.76eV, 0.065u"1/Å^3")
  nitride = Material(14amu, 7, 0eV, 1eV, 0eV, 0.065u"1/Å^3")
  target = Disk(boron + nitride, 100um, 2u"inch")
  opts = Options(name="hydrogen on boron-nitride", track_trajectories=true)
  @benchmark begin
    p = BCA.default_incident(1.008amu, 1, 1000eV, 1eV, 10eV, 0um, [0.9999999999984769, 1.7453292519934434e-6, 0.0], track_trajectories=true)
    result1 = BCA.singleionbca(p, $target, $opts)
  end

  # WITHOUT CONST UNITS
  # BenchmarkTools.Trial: 83 samples with 1 evaluation.
  #  Range (min … max):  10.823 ms … 112.659 ms  ┊ GC (min … max): 0.00% … 0.00%
  #  Time  (median):     62.239 ms               ┊ GC (median):    0.00%
  #  Time  (mean ± σ):   60.614 ms ±  19.489 ms  ┊ GC (mean ± σ):  1.13% ± 4.48%
  #
  #                               ▂     ▂█         ▄                
  #   ▄▁▁▄▁▁▄▄▁▁▁▄▁▄▄█▁▆▁▄▁▆▆▄▄▁▆▆█▆▆▆▄███▆▆▄█▆█▄▆▆█▄▄▄▁▄▆▄▄▁▄▄▁▁▄ ▁
  #   10.8 ms         Histogram: frequency by time         98.7 ms <
  #
  #  Memory estimate: 796.48 KiB, allocs estimate: 31622.

  # WITH CONST UNITS
  # julia> BenchmarkTools.Trial: 84 samples with 1 evaluation.
  #  Range (min … max):  11.339 ms … 93.135 ms  ┊ GC (min … max): 0.00% … 0.00%
  #  Time  (median):     61.120 ms              ┊ GC (median):    0.00%
  #  Time  (mean ± σ):   59.693 ms ± 16.665 ms  ┊ GC (mean ± σ):  1.03% ± 2.97%
  #
  #                               ▂      ▂█ ▂█▂▅▂ ▅  ▅ ▂   ▅  ▂    
  #   ▅▅▁▁▁▁▅▅▁▁▁█▁▁▁▁█▁▁▅▁▁▁▅▅▁▅▁██▁██████▅████████▁█▅█▅▅███▁█▁█ ▁
  #   11.3 ms         Histogram: frequency by time        85.1 ms <
  #
  #  Memory estimate: 888.23 KiB, allocs estimate: 34972.

  # WITH CONST EVERYTHING
  # BenchmarkTools.Trial: 87 samples with 1 evaluation.
  #  Range (min … max):  10.824 ms … 126.559 ms  ┊ GC (min … max): 0.00% … 0.00%
  #  Time  (median):     58.351 ms               ┊ GC (median):    0.00%
  #  Time  (mean ± σ):   57.506 ms ±  20.547 ms  ┊ GC (mean ± σ):  1.01% ± 4.32%
  #
  #                            ▆▃ ▃▆ ▁▃▃▃█▁▃   ▃                    
  #   ▇▁▁▄▁▄▄▁▇▄▄▄▄▁▁▄▁▁▄▄▇▄▁▄▁██▇██▄███████▄▁▄█▁▄▇▄▄▁▁▁▄▇▁▄▁▁▁▇▁▇ ▁
  #   10.8 ms         Histogram: frequency by time         99.2 ms <
  #
  #  Memory estimate: 899.38 KiB, allocs estimate: 36396.

end

# using Logging
# debuglogger = ConsoleLogger(stderr, Logging.Debug)
# global_logger(debuglogger)

# infologger = ConsoleLogger(stderr, Logging.Info)
# global_logger(infologger)
