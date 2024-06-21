#=
Script to visualize stress paths of reservoir points, and plastic cap evolution
=#

cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__)

using LinearAlgebra
using HDF5
using DelimitedFiles
using CairoMakie
CairoMakie.activate!()

using Revise
lib_dir = "C:/Users/lsalo/julia/jl-tools/"
includet(joinpath(lib_dir, "modules", "arrayOps.jl"))
using .ArrayOps
includet(joinpath(lib_dir, "modules", "abaqus.jl"))
using .Abaqus
includet(joinpath(lib_dir, "modules", "fileio.jl"))
using .Fileio


# Directories and options
tstep = 2629746 / 31556952
t = 1936 .+ vcat(0, cumsum(repeat([tstep], inner=469)))
t = t[vcat(1, 13:12:length(t))]

nruns = 4
uz_top_max = zeros(length(t), nruns)
for run=1:nruns
    abaqus_dir = join(["C:/Users/lsalo/abaqus/simulation_data/wilmington/plastic-", run, "/extracted/U/TOPBC/"])
    print(abaqus_dir)
    abaqus_fn = readdir(abaqus_dir)
    abaqus_fn = sort(abaqus_fn, by=sort_steps)

    # Get rebound
    # Top node of max positive u
    nsteps = length(abaqus_fn)
    for tstep in 1:nsteps
        file = readdlm(join([abaqus_dir, abaqus_fn[tstep]]), ',', skipstart=1)
        uz_top_max[tstep, run] = maximum(file[:, end])
    end
end

dir_in = "C:/Users/lsalo/julia/output/wilmington/"
figs_dir = "C:/Users/lsalo/julia/graphics/figs/wilmington/test/"

# -------------------------------------- PLOTS ---------------------------------------------------------
# 
# Plot 2: (a) max(u_z) vs t, comparison with subsidence data
#         (b) u_z top + faults view
#
#data_u3max = [1947 1951 1954 1958 1965; -2.4 -4 -6.1 -7.9 -8.8]        # Kosloff et al., Tectonophysics (1980)
data_uzmax = [1975; 0.45]                                               # Kosloff et al., Tectonophysics (1980)
f = Figure(size=(600, 300))
ax = Axis(f[1, 1], xlabel = L"t \; \text{ [y]}", ylabel = L"\max \; u_z \; \text{ [m]}", xlabelsize = 18, ylabelsize = 18, xticks=1935:5:1980, yticks =0:.5:4)
sc1 = scatter!(ax, data_uzmax[1, :], data_uzmax[2, :], markersize=10, color=RGBAf(1,1,1,0), marker=:rect, strokewidth=1, strokecolor=:black, label=L"\text{ Field data}")
ln1 = lines!(ax, t[1:37], uz_top_max[1:37, 1], color=:black, linewidth=1, label=L"\text{Plastic-1}")
ln1 = lines!(ax, t[1:37], uz_top_max[1:37, 2], color=:blue, linewidth=1, label=L"\text{Plastic-2}")
ln1 = lines!(ax, t[1:37], uz_top_max[1:37, 3], color=:green, linewidth=1, label=L"\text{Plastic-3}")
ln1 = lines!(ax, t, uz_top_max[:, 4], color=:red, linewidth=1, label=L"\text{Plastic-4}")
xlims!(ax, 1935, 1980)
ylims!(ax, 0, 4)
axislegend(position = :lt, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8)) 
f
save(join([figs_dir, "u3_surface_rebound_plastic-1_plastic-4.pdf"]), f, pt_per_unit = 1)