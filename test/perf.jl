using Revise
using Test
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


boron = Material(10.811amu, 5, 0eV, 1eV, 5.76eV, 0.065u"1/Å^3")
nitride = Material(14amu, 7, 0eV, 1eV, 0eV, 0.065u"1/Å^3")
target = Disk(boron + nitride, 100um, 2u"inch")
opts = Options(name="hydrogen on boron-nitride", track_trajectories=true)
p = BCA.default_incident(1.008amu, 1, 1000eV, 1eV, 10eV, 0um, [0.9999999999984769, 1.7453292519934434e-6, 0.0], track_trajectories=true)
result1 = BCA.singleionbca(p, target, opts)

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

# WITH typeof(1.0u"eV") and friends
# BenchmarkTools.Trial: 151 samples with 1 evaluation.
#  Range (min … max):   2.723 ms … 58.826 ms  ┊ GC (min … max): 0.00% … 0.00%
#  Time  (median):     34.141 ms              ┊ GC (median):    0.00%
#  Time  (mean ± σ):   33.169 ms ± 10.130 ms  ┊ GC (mean ± σ):  1.15% ± 3.97%
#
#                                     ▄▁▇  ▄▂▄ ▂  █ ▅▁  ▂        
#   ▃▁▁▅▁▃▁▁▃▅▁▅▁▃▃▃▅▁▁▁█▁▁▁▃▁█▅▅▆▃██▆████▆█████▅▆█▅██▃▆█▆▃▁▃▃▆ ▃
#   2.72 ms         Histogram: frequency by time        49.6 ms <
#
#  Memory estimate: 258.58 KiB, allocs estimate: 8723.

# MEAN FREE PATH FIX
# BenchmarkTools.Trial: 316 samples with 1 evaluation.
#  Range (min … max):  523.277 μs … 37.877 ms  ┊ GC (min … max): 0.00% … 24.93%
#  Time  (median):      16.204 ms              ┊ GC (median):    0.00%
#  Time  (mean ± σ):    15.817 ms ±  5.950 ms  ┊ GC (mean ± σ):  1.88% ±  5.89%
#
#                           ▁ ▃▁▃█▃▇▃▁▆▂▁▃ ▃                      
#   ▄▃▃▃▃▃▃▄▄▃▄▆▆▅▃▅▅▆▃▅▄▃▅▅█████████████████▃█▄▅▅▅▄▃▅▃▃▃▃▃▁▁▃▁▃ ▄
#   523 μs          Histogram: frequency by time         30.4 ms <
#
#  Memory estimate: 179.11 KiB, allocs estimate: 6082.

# TARGET GET PROPERTY FIXES
# Range (min … max):  219.533 μs … 20.290 ms  ┊ GC (min … max): 0.00% … 61.08%
# Time  (median):       4.834 ms              ┊ GC (median):    0.00%
# Time  (mean ± σ):     4.883 ms ±  2.105 ms  ┊ GC (mean ± σ):  3.86% ±  9.17%
#
#                     ▃▃▄▄▆█▅▂▁                                  
#  ▃▃▄▄▄▄▄▅▃▄▄▃▄▅▄▆▅▆████████████▅▅▄▅▄▃▃▃▃▃▃▃▂▂▂▂▃▂▂▂▁▃▂▃▁▂▃▃▃▃ ▄
#  220 μs          Histogram: frequency by time           12 ms <
#
# Memory estimate: 87.55 KiB, allocs estimate: 2390.

# CUT DYNAMIC DISPATCH
# BenchmarkTools.Trial: 2373 samples with 1 evaluation.
#  Range (min … max):  41.190 μs … 15.742 ms  ┊ GC (min … max): 0.00% … 80.36%
#  Time  (median):      2.083 ms              ┊ GC (median):    0.00%
#  Time  (mean ± σ):    2.100 ms ±  1.098 ms  ┊ GC (mean ± σ):  4.31% ±  8.51%
#
#              ▁▃██▆▄▁                                           
#   ▃▃▃▄▄▅▄▄▅▅▇███████▇▄▄▃▃▃▂▂▂▂▂▂▂▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▂▁▂▂▂ ▃
#   41.2 μs         Histogram: frequency by time        8.59 ms <
#
#  Memory estimate: 35.94 KiB, allocs estimate: 485.

using Profile
using PProf
Profile.clear()
using ProfileCanvas

boron = Material(10.811amu, 5, 0eV, 1eV, 5.76eV, 0.065u"1/Å^3")
nitride = Material(14amu, 7, 0eV, 1eV, 0eV, 0.065u"1/Å^3")
target = Disk(boron + nitride, 100ang, 2u"inch")
opts = Options(name="hydrogen on boron-nitride", track_trajectories=true, track_recoil_trajectories=true)

# H = (BCA.default_incident(1.008amu, 1, 1000eV, 1eV, 10eV, 0um, [0.9999999999984769, 1.7453292519934434e-6, 0.0], track_trajectories=true) for _ in 1:20)
# B = (BCA.default_incident(1.008amu, 5, 1000eV, 1eV, 10eV, 0um, [0.9999999999984769, 1.7453292519934434e-6, 0.0], track_trajectories=true) for _ in 1:20)
# Zn = (BCA.default_incident(1.008amu, 30, 9000eV, 1eV, 10eV, 0um, [0.9999999999984769, 1.7453292519934434e-6, 0.0], track_trajectories=true) for _ in 1:20)
# ps = Base.Iterators.flatten((H, B, Zn))

H = [BCA.default_incident(1.008amu, 1, 1000eV, 1eV, 10eV, 0um, [0.9999999999984769, 1.7453292519934434e-6, 0.0], track_trajectories=true) for _ in 1:20]
B = [BCA.default_incident(1.008amu, 5, 1000eV, 1eV, 10eV, 0um, [0.9999999999984769, 1.7453292519934434e-6, 0.0], track_trajectories=true) for _ in 1:20]
Zn = [BCA.default_incident(1.008amu, 30, 9000eV, 1eV, 10eV, 0um, [0.9999999999984769, 1.7453292519934434e-6, 0.0], track_trajectories=true) for _ in 1:200]
ps = vcat(H, B, Zn)

# profview
ps = [BCA.default_incident(1.008amu, 30, 9000eV, 1eV, 10eV, 0um, [0.9999999999984769, 1.7453292519934434e-6, 0.0], track_trajectories=true) for _ in 1:5000];
@profview BCA.run(ps, target, opts)

# profview single thread
ps = [BCA.default_incident(1.008amu, 30, 9000eV, 1eV, 10eV, 0um, [0.9999999999984769, 1.7453292519934434e-6, 0.0], track_trajectories=true) for _ in 1:5000];
@profview BCA.bca(ps, target, opts)

# profile
ps = [BCA.default_incident(1.008amu, 30, 9000eV, 1eV, 10eV, 0um, [0.9999999999984769, 1.7453292519934434e-6, 0.0], track_trajectories=true) for _ in 1:5000];
@profile BCA.bca(ps, target, opts)
pprof(webport=1113)

@code_warntype BCA.bca(ps, target, opts)

using JET
@report_opt BCA.bca(ps, target, opts)

# just run
ps = [BCA.default_incident(1.008amu, 30, 9000eV, 1eV, 10eV, 0um, [0.9999999999984769, 1.7453292519934434e-6, 0.0], track_trajectories=true) for _ in 1:5000];
result = BCA.run(ps, target, opts);


# using Logging
# debuglogger = ConsoleLogger(stderr, Logging.Debug)
# global_logger(debuglogger)

# infologger = ConsoleLogger(stderr, Logging.Info)
# global_logger(infologger)
