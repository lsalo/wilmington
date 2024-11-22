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

nruns = 1
idx_run = [7]
uz_top_subsMax = zeros(length(t), nruns)
uz_top_rebMax = zeros(length(t), nruns)
uz_top_subs = zeros(length(t), nruns)
for run=1:nruns
    abaqus_dir = join(["C:/Users/lsalo/abaqus/simulation_data/wilmington/plastic-", idx_run[run], "/extracted/U/TOPBC/"])
    print(abaqus_dir)
    abaqus_fn = readdir(abaqus_dir)
    abaqus_fn = sort(abaqus_fn, by=sort_steps)

    # get values from t=300 first
    id_t300 = findfirst(abaqus_fn .== "step_300_nodes.dat")
    file = readdlm(join([abaqus_dir, abaqus_fn[id_t300]]), ',', skipstart=1)
    uz_top_300 = file[:, end]
    id_minUz_t300 = argmin(uz_top_300)

    # Get max subsidence and rebound
    nsteps = length(abaqus_fn)
    get_rebound = false
    for tstep in 1:nsteps
        file = readdlm(join([abaqus_dir, abaqus_fn[tstep]]), ',', skipstart=1)

        # Max subsidence and rebound
        uz_top_subsMax[tstep, run] = minimum(file[:, end])
        if get_rebound
            uz_top_rebMax[tstep, run] = maximum(file[:, end] - uz_top_300)
        end
        if abaqus_fn[tstep] == "step_300_nodes.dat"
            get_rebound = true
        end

        # Subsidence/rebound history at node of maximum subsidence at t=1961
        uz_top_subs[tstep, run] = file[id_minUz_t300, end]
    end

end

dir_in = "C:/Users/lsalo/julia/output/wilmington/"
figs_dir = "C:/Users/lsalo/julia/graphics/figs/wilmington/test/"

# -------------------------------------- PLOTS ---------------------------------------------------------
# 
# Subs and rebound
field_data = [1947 1951 1954 1958 1965 1975; -2.4 -4 -6.1 -7.9 -8.8 -8.3] 
f = Figure(size=(600, 300))
ax = Axis(f[1, 1], xlabel = L"t \; \text{ [y]}", ylabel = L" u_z \; \text{ [m]}", xlabelsize = 18, ylabelsize = 18, xticks=1935:5:1980, yticks =-10:2:0)
sc1 = scatter!(ax, field_data[1, :], field_data[2, :], markersize=10, color=RGBAf(1,1,1,0), marker=:rect, strokewidth=1, strokecolor=:black, label=L"\text{ Field data}")
ln1 = lines!(ax, t, uz_top_subs[:, 1], color=:black, linewidth=1, label=L"\text{Plastic-7}")
xlims!(ax, 1935, 1980)
ylims!(ax, -10, 0)
axislegend(position = :lt, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8)) 
f

# Subs and rebound, discrete points
field_data = [1958 1965 1975;-7.9 -8.8 -7.4] 
f = Figure(size=(600, 300))
ax = Axis(f[1, 1], xlabel = L"t \; \text{ [y]}", ylabel = L" u_z \; \text{ [m]}", xlabelsize = 18, ylabelsize = 18, xticks=1935:5:1980, yticks =-16:2:0)
sc1 = scatter!(ax, field_data[1, :], field_data[2, :], markersize=12, color=RGBAf(1,1,1,0), marker=:rect, strokewidth=1, strokecolor=:black, label=L"\text{ Field data}")
sc2 = scatter!(ax, [1961, 1975], [-6.2, -0.77], markersize=10, color=RGBAf(1,1,1,0), marker=:circle, strokewidth=1, strokecolor=:blue, label=L"\text{Plastic-7}")
sc3 = scatter!(ax, [1961, 1975], [-5.7, -1.34], markersize=10, color=RGBAf(1,1,1,0), marker=:circle, strokewidth=1, strokecolor=:red, label=L"\text{Plastic-14}")
sc4 = scatter!(ax, [1961, 1975], [-15.4, -3.4], markersize=10, color=RGBAf(1,1,1,0), marker=:circle, strokewidth=1, strokecolor=:green, label=L"\text{Plastic-31}")
sc5 = scatter!(ax, [1961, 1975], [-8.7, -1.5], markersize=10, color=RGBAf(1,1,1,0), marker=:circle, strokewidth=1, strokecolor=:cyan, label=L"\text{Plastic-34}")
sc6 = scatter!(ax, [1961, 1975], [-7.3, -3], markersize=10, color=RGBAf(1,1,1,0), marker=:circle, strokewidth=1, strokecolor=:magenta, label=L"\text{Plastic-41}")
sc7 = scatter!(ax, [1961, 1975], [-7.9, -2.8], markersize=12, color=:orange, marker=:circle, strokewidth=1, strokecolor=:orange, label=L"\text{Plastic-48}")
sc8 = scatter!(ax, [1961, 1975], [-8.3, -3.3], markersize=12, color=:purple, marker=:circle, strokewidth=1, strokecolor=:purple, label=L"\text{Plastic-54}")
xlims!(ax, 1935, 1980)
ylims!(ax, -16, 0)
axislegend(position = :lt, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8)) 
f
save(join([figs_dir, "u3_tests.pdf"]), f, pt_per_unit = 1)

# Plot 2: (a) max(u_z) vs t, comparison with subsidence data
#         (b) u_z top + faults view
#
#data_u3max = [1947 1951 1954 1958 1965; -2.4 -4 -6.1 -7.9 -8.8]        # Kosloff et al., Tectonophysics (1980)
data_uzmax = [1975; 0.45]                                               # Kosloff et al., Tectonophysics (1980)
f = Figure(size=(600, 300))
ax = Axis(f[1, 1], xlabel = L"t \; \text{ [y]}", ylabel = L"\max \; u_z \; \text{ [m]}", xlabelsize = 18, ylabelsize = 18, xticks=1935:5:1980, yticks =0:.5:4)
sc1 = scatter!(ax, data_uzmax[1, :], data_uzmax[2, :], markersize=10, color=RGBAf(1,1,1,0), marker=:rect, strokewidth=1, strokecolor=:black, label=L"\text{ Field data}")
ln1 = lines!(ax, t[1:37], uz_top_max[1:37, 1], color=:black, linewidth=1, label=L"\text{Plastic-7}")
#ln1 = lines!(ax, t[1:37], uz_top_max[1:37, 2], color=:blue, linewidth=1, label=L"\text{Plastic-2}")
#ln1 = lines!(ax, t[1:37], uz_top_max[1:37, 3], color=:green, linewidth=1, label=L"\text{Plastic-3}")
#ln1 = lines!(ax, t, uz_top_max[:, 4], color=:red, linewidth=1, label=L"\text{Plastic-4}")
xlims!(ax, 1935, 1980)
ylims!(ax, 0, 4)
axislegend(position = :lt, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8)) 
f
#save(join([figs_dir, "u3_surface_rebound_plastic-1_plastic-4.pdf"]), f, pt_per_unit = 1)