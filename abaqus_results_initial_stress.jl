#=
Script to visualize stresses at elements (we import values at integration points from abaqus)
=#

# Startup
cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__)

using HDF5
using DelimitedFiles
using CairoMakie
CairoMakie.activate!()

# Directories and options
abaqus_dir = "C:/Users/lsalo/abaqus/simulation_data/wilmington/output/plastic_ini/extracted/S/ENTIRE_DOMAIN/"
dir_in = "C:/Users/lsalo/julia/output/wilmington/"
figs_dir = "C:/Users/lsalo/julia/graphics/figs/wilmington/test/" 

abaqus_fn = "Geostatic.h5"

# Read stress data
file = h5open(join([abaqus_dir, abaqus_fn]), "r")
abq_out = [read(file[name]) for name in ["Sxx", "Syy", "Szz", "Sxy", "Sxz", "Syz"]]
abq_out = hcat(abq_out...)
stress_ini = readdlm(join([dir_in, "initial_stress_rotated_criticallyStressed_mu0.6.inp"]), ',')

# Plot comparison
fig_size = (900, 400)
f = Figure(size=fig_size)
ax = Axis(f[1, 1], xlabel = L"\sigma^\text{eff}_{11} \text{ (I.C.)} \; \text{ [MPa]}", ylabel = L"\sigma^\text{eff}_{11} \text{ (Geostatic)} \; \text{ [MPa]}", xlabelsize = 16, ylabelsize = 16)
sc1 = scatter!(ax, 1e-6*stress_ini[:, 2], 1e-6*abq_out[:, 1], markersize = 6, color=RGBAf(0,0,0,0), strokecolor=:black, strokewidth=0.5, label=L"\text{ Data}")
ln1 = lines!(ax, [-650, 0], [-650, 0], color = RGBAf(0.8, 0.8, 0.8, 1), linewidth = 1.25, linestyle = :solid, label=L" 1:1")
xlims!(ax, -600, 0)
ylims!(ax, -600, 0)
axislegend(position = :lt, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8)) 

ax = Axis(f[1, 2], xlabel = L"\sigma^\text{eff}_{22} \text{ (I.C.)} \; \text{ [MPa]}", ylabel = L"\sigma^\text{eff}_{22} \text{ (Geostatic)} \; \text{ [MPa]}", xlabelsize = 16, ylabelsize = 16)
sc1 = scatter!(ax, 1e-6*stress_ini[:, 3], 1e-6*abq_out[:, 2], markersize = 6, color=RGBAf(0,0,0,0), strokecolor=:blue, strokewidth=0.5, label=L"\text{ Data}")
ln1 = lines!(ax, [-250, 0], [-250, 0], color = RGBAf(0.8, 0.8, 0.8, 1), linewidth = 1.25, linestyle = :solid, label=L" 1:1")
xlims!(ax, -211, 0)
ylims!(ax, -211, 0)

ax = Axis(f[1, 3], xlabel = L"\sigma^\text{eff}_{33} \text{ (I.C.)} \; \text{ [MPa]}", ylabel = L"\sigma^\text{eff}_{33} \text{ (Geostatic)} \; \text{ [MPa]}", xlabelsize = 16, ylabelsize = 16)
sc1 = scatter!(ax, 1e-6*stress_ini[:, 4], 1e-6*abq_out[:, 3], markersize = 6, color=RGBAf(0,0,0,0), strokecolor=:red, strokewidth=0.5, label=L"\text{ Data}")
ln1 = lines!(ax, [-250, 0], [-250, 0], color = RGBAf(0.8, 0.8, 0.8, 1), linewidth = 1.25, linestyle = :solid, label=L" 1:1")
xlims!(ax, -191, 0)
ylims!(ax, -191, 0)
f

target_dpi=300 
target_size_in=(9,6)
ppu = target_size_in[1]*target_dpi/fig_size[1]
save(join([figs_dir, "ic_vs_geostatic.png"]), f, px_per_unit = ppu)