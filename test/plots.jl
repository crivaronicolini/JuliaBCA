using Test
using BCA
using Unitful
using DimensionfulAngles
using DataFrames
using Plots

const amu = u"u"
const eV = u"eV"
const um = u"μm"
const cm = u"cm"
const rad = u"radᵃ"
const ang = u"Å"

boron = Material(10.811amu, 5, 0eV, 1eV, 5.76eV, 0.065u"1/Å^3")
nitride = Material(14amu, 7, 0eV, 1eV, 0eV, 0.065u"1/Å^3")
target = Disk(boron + nitride, 200um, 2u"inch") + Disk(nitride, 100um, 2u"inch")
opts = Options(name="H,Zn on boron-nitride, nitride", track_trajectories=true, track_recoil_trajectories=true)

len = 100
Zn = [BCA.default_incident(1.008amu, 30, 9000eV, 1eV, 10eV, 0um, [0.9999999999984769, 1.7453292519934434e-6, 0.0], track_trajectories=true) for _ in 1:len];
H = [BCA.default_incident(1.008amu, 1, 1000eV, 1eV, 10eV, 0um, [0.9999999999984769, 1.7453292519934434e-6, 0.0], track_trajectories=true) for _ in 1:len];
Fe = [BCA.default_incident(1.008amu, 26, 1000eV, 1eV, 10eV, 0um, [0.9999999999984769, 1.7453292519934434e-6, 0.0], track_trajectories=true) for _ in 1:len];
ps = vcat(H, Zn, Fe);

result = BCA.run(ps, target, opts);

trajectoryplot(result, target)

df = DataFrame(result)
sputtered, reflected, implanted = BCA.destinations(df)
stephist(implanted.x .|> ang, title="implanted Ions", bins=100, density=true)
stephist(sputtered.E, title="Sputtered Ions Enegy", bins=100, density=true)


# one of every element
ps = vcat([[BCA.default_incident(1.008amu, Z, 1000eV, 1eV, 10eV, 0um, [0.9999999999984769, 1.7453292519934434e-6, 0.0], track_trajectories=true) for Z in 1:155] for _ in 1:50]...);
result = BCA.run(ps, target, opts);

trajectoryplot(result, target)

sputtered, reflected, implanted = BCA.destinations(result)
df = DataFrame(result)
sputtered, reflected, implanted = destinations(df)
stephist(implanted.x .|> ang, title="implanted Ions", bins=100, density=true)
stephist(sputtered.E, title="Sputtered Ions Enegy", bins=100, density=true)


