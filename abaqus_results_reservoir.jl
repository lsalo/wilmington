#=
Script to visualize point displacements, stress paths of reservoir points, and plastic cap evolution
=#

cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__)

using Distributed
addprocs(12)
@everywhere using SharedArrays

using LinearAlgebra
@everywhere using HDF5
using DelimitedFiles
using CairoMakie
CairoMakie.activate!()

@everywhere  lib_dir = "C:/Users/lsalo/julia/jl-tools/"
@everywhere include(joinpath(lib_dir, "modules", "arrayOps.jl"))
@everywhere using .ArrayOps
include(joinpath(lib_dir, "modules", "abaqus.jl"))
using .Abaqus
include(joinpath(lib_dir, "modules", "fileio.jl"))
using .Fileio


## Directories and options
sim_name = "plastic-92"
abaqus_dir_el = "C:/Users/lsalo/abaqus/simulation_data/wilmington/plastic-92/extracted/S/ENTIRE_DOMAIN/"
abaqus_fn_el = sort(readdir(abaqus_dir_el), by=sort_steps)

abaqus_dir_nod = "C:/Users/lsalo/abaqus/simulation_data/wilmington/plastic-92/extracted/U/TOPRANGERSURFACE/"
abaqus_fn_nod = sort(readdir(abaqus_dir_nod), by=sort_steps)
abaqus_dir_nod2 = "C:/Users/lsalo/abaqus/simulation_data/wilmington/plastic-92/extracted/S/TOPRANGERSURFACE/"
abaqus_fn_nod2 = sort(readdir(abaqus_dir_nod2), by=sort_steps)
abaqus_dir_nod3 = "C:/Users/lsalo/abaqus/simulation_data/wilmington/plastic-92/extracted/U/TOPBC/"
abaqus_fn_nod3 = sort(readdir(abaqus_dir_nod3), by=sort_steps)

dir_inp = "C:/Users/lsalo/Wilmington/AbaqusGrids/"
fn_inp  = "WilmingtonMesh_AllFaults_Rotated_Coarse.inp"
dir_in = "C:/Users/lsalo/julia/output/wilmington/"
figs_dir = "C:/Users/lsalo/julia/graphics/figs/wilmington/test/"

tn = 470
tstep = 2629746 / 31556952
t = 1936 .+ vcat(0, cumsum(repeat([tstep], inner=470)))
t = t[vcat(1, 13:12:length(t))]
node_select = "max_subsidence"


## Get nodal data
# Max diplacement node at top Ranger surface
u_topranger_fin = readdlm(join([abaqus_dir_nod, abaqus_fn_nod[end]]), ',', skipstart=1)
if node_select == "max_subsidence"
    id_u = argmin(u_topranger_fin[:, end])
elseif node_select == "max_rebound"
    id_u = argmax(u_topranger_fin[:, end])
end
nod_u = Int(u_topranger_fin[id_u, 1])

n = length(abaqus_fn_nod)
u_z = Vector{Float32}(undef, n)
for i in 1:n
    file = readdlm(join([abaqus_dir_nod, abaqus_fn_nod[i]]), ',', skipstart=1)
    u_z[i] = file[id_u, end]
end

# Max displacement at the surface (may not be the same node in all steps)
u_z_top = Vector{Float32}(undef, n)
for i in 1:n
    file = readdlm(join([abaqus_dir_nod3, abaqus_fn_nod3[i]]), ',', skipstart=1)
    if node_select == "max_subsidence"
        u_z_top[i] = minimum(file[:, end])
    elseif node_select == "max_rebound"
        u_z_top[i] = maximum(file[:, end])
    end
end

#s_33 = Vector{Int32}(undef, n)
s_33 = SharedArray{Int32}(n)
@distributed for i in 1:n
    file = h5open(join([abaqus_dir_nod2, abaqus_fn_nod2[i]]), "r")
    abq_out = [read(file[name]) for name in ["NODE_LABEL", "Sxx", "Syy", "Szz", "Sxy", "Sxz", "Syz"]]
    abq_out = Int.(hcat(abq_out...))
    s_nod = reduce_matrix_with_repeated_ids(abq_out, "mean") # currently very slow, needs to be revised / run in parallel
    s_33[i] = s_nod[id_u, 4]
    print("$i/$n completed\n")
end


## Get mesh info (id of Ranger elements containing selected node)
# Note that we assume that Ranger element ids start at 1, and data are sorted
nodes_inp, elements_inp, elset_names_inp, nset_inp, nset_names_inp = get_connectivity_from_inp(dir_inp, fn_inp, "tetrahedral")
idx = vec(any(elements_inp[1][:, 2:5] .== nod_u, dims=2))
id_elements = elements_inp[1][idx, 1]


## Get Ranger element with closest centroid to node above, and stress evolution
centroids = readdlm(join([dir_in, "element_centroids.dat"]), ',')
centroids_els = centroids[id_elements, 2:end]
coord_nod_u = nodes_inp[findfirst(nodes_inp[:, 1] .== nod_u), 2:end]
dist = mapslices(row -> norm(row - coord_nod_u), centroids_els, dims=2)
id_element = id_elements[argmin(dist)]

n = length(abaqus_fn_el)
s = Vector{}(undef, n)
p_eff = Vector{Float32}(undef, n)
q = Vector{Float32}(undef, n)
for i in 1:n
    file = h5open(join([abaqus_dir_el, abaqus_fn_el[i]]), "r")
    abq_out = [read(file[name]) for name in ["Sxx", "Syy", "Szz", "Sxy", "Sxz", "Syz"]]
    abq_out = hcat(abq_out...)
    s[i] = abq_out[id_element, :]
    p_eff[i] = -(1/3)*sum(s[i][1:3]) 
    q[i] = sqrt( 0.5*((s[i][1] - s[i][2])^2 + (s[i][2] - s[i][3])^2 + (s[i][3] - s[i][1])^2) + 6*(s[i][4]^2 + s[i][5]^2 + s[i][6]^2))
end


# -------------------------------------- PLOTS ---------------------------------------------------------
# 
# Plot 1: (a) stress path of ranger cell closest to node with specified u_z (e.g., max(u_z))
#         (b) u_z vs s33 (top ranger)

regime = "strike-slip"

# Plot parameters
# Cap model parameters
d = 1e-3        # MPa
b = 55        # 55, 58.6, 52.41 deg
R = 0.3         # 1.227
W = 0.27        # [-]
D = 0.7143      # 1/kbar
e_p = [0, 0.0025, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2]

# Hardening curve [Pa]
# p_eff_el below should be very close to p_eff[1] above
if regime == "reverse"
    s_ini_el = readdlm(join([dir_in, "initial_stress_rotated_criticallyStressed_mu0.6.inp"]), ',')
elseif regime == "normal"
    s_ini_el = readdlm(join([dir_in, "initial_stress_rotated_criticallyStressed_mu0.6_normal.inp"]), ',')
elseif regime == "strike-slip"
    s_ini_el = readdlm(join([dir_in, "initial_stress_rotated_criticallyStressed_mu0.6_strike-slip.inp"]), ',')
end
p_eff_ini_el = -(1/3)*sum(s_ini_el[id_element, 2:4])
#p_b = 1.18*p_eff_ini_el .+ 0.25*((1/3)*(1/D)*log.(e_p./W .+ 1)*1e8) # plastic-7
#p_b = 1.19*p_eff_ini_el .+ 0.25*((1/3)*(1/D)*log.(e_p./W .+ 1)*1e8) # plastic-41
#p_b = 1.19*p_eff_ini_el .+ 0.1*((1/3)*(1/D)*log.(e_p./W .+ 1)*1e8) # plastic-62
#p_b = 1.18*p_eff_ini_el .+ 0.1*((1/3)*(1/D)*log.(e_p./W .+ 1)*1e8) # plastic-70, plastic-77
#p_b = 1.17*p_eff_ini_el .+ 0.1*((1/3)*(1/D)*log.(e_p./W .+ 1)*1e8) # plastic-71
#p_b = 1.15*p_eff_ini_el .+ 0.1*((1/3)*(1/D)*log.(e_p./W .+ 1)*1e8) # plastic-73
#p_b = 1.08*p_eff_ini_el .+ 0.1*((1/3)*(1/D)*log.(e_p./W .+ 1)*1e8) # plastic-75
#p_b = 1.41*p_eff_ini_el .+ 0.1*((1/3)*(1/D)*log.(e_p./W .+ 1)*1e8) # plastic-75
#p_b = 1.02*p_eff_ini_el .+ 0.1*((1/3)*(1/D)*log.(e_p./W .+ 1)*1e8) # plastic-78
#p_b = 1.35*p_eff_ini_el .+ 0.1*((1/3)*(1/D)*log.(e_p./W .+ 1)*1e8)  # plastic-80
#p_b = 1.35*p_eff_ini_el .+ 0.6*((1/3)*(1/D)*log.(e_p./W .+ 1)*1e8)  # plastic-84
#p_b = 1.35*p_eff_ini_el .+ 0.7*((1/3)*(1/D)*log.(e_p./W .+ 1)*1e8)  # plastic-85, plastic-88, plastic-89
p_b = 1.2*p_eff_ini_el .+ 0.25*((1/3)*(1/D)*log.(e_p./W .+ 1)*1e8)  # plastic-92

# Shear envelope and cap
e_p0 = 0
p_Fs = Vector{}(0:10:50)
q_Fs = p_Fs.*tand(b) .+ d

p_a = (1e-6*p_b[1] - R*d) / (1 + R*tand(b))             # pb[1] = Initial cap position
if node_select == "max_subsidence"
    #q_Fc = vcat(Vector{}(0:1:27), 27.01:0.01:28.0)   # p_b[1], 
    #q_Fc = vcat(Vector{}(0:1:35), 35.01:0.01:36.8)   # p_b[1], b = 58
    #q_Fc = vcat(Vector{}(0:1:29), 29.02:.02:29.16)   
    #q_Fc = vcat(Vector{}(0:1:8), 8.02:.02:9.66)     # normal regime  
    q_Fc = vcat(Vector{}(0:1:19), 19.02:.02:19.4)    # strike-slip regime
elseif node_select == "max_rebound"
    #q_Fc = vcat(Vector{}(0:1:34), 34.1:.1:34.9)         # p_b[1], b = 52
    #q_Fc = vcat(Vector{}(0:1:40), 40.02:.02:40.7)         # p_b[1],  b= 58
    q_Fc = vcat(Vector{}(0:1:12), 12.02:.02:13.04)
end
p_Fc = p_a .+ sqrt.( abs.((R*(d + p_a*tand(b)))^2 .- (R.*q_Fc).^2) )


# Make plot stress path
f = Figure(size=(800, 350))
ax = Axis(f[1, 1], xlabel = L"p^\text{eff} = -\frac{\sigma_{kk}^\text{eff}}{3} \; \text{ [MPa]}", ylabel = L"q = \sqrt{3 J_2} \; \text{ [MPa]}", xlabelsize = 16, ylabelsize = 16)
ln1 = lines!(ax, p_Fs, q_Fs, color = RGBAf(0, 0, 0, 1), linewidth = 1.5, label=L"F_\text{s},\; \beta=55.0^\circ")
ln1 = lines!(ax, p_Fc, q_Fc, color = RGBAf(0.5, 0.5, 0.5, 1), linewidth = 1, linestyle = :solid, label=L"F_\text{c}^0")
sc1 = scatter!(ax, p_eff/1e6, q/1e6, markersize = 8, color=t, colormap=:rainbow, label=L"\text{Ranger element}")
if regime == "reverse"
    xlims!(ax, 0, 45)
    ylims!(ax, 0, 45)
elseif regime == "normal"
    xlims!(ax, 0, 16)
    ylims!(ax, 0, 16)
elseif regime == "strike-slip"
    xlims!(ax, 0, 30)
    ylims!(ax, 0, 30)
end
axislegend(position = :lt, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8), labelsize=14) 
cbar = Colorbar(f[1, 2], sc1, label=L"\text{Year}", labelsize=12, width=12, ticklabelsize=12)

#ax = Axis(f[1, 2], xlabel = L"e^\text{in}_\text{v} \; \text{ [-]}", ylabel = L"p_\text{b} \; \text{ [MPa]}", xlabelsize = 16, ylabelsize = 16)
#ln1 = lines!(ax, e_p, p_b*1e-6, color = RGBAf(0, 0, 0, 1), linewidth = 1.5)
#xlims!(ax, 0, 0.2)
#ylims!(ax, 0, 75)
if regime == "reverse"
    #ax = Axis(f[1, 3], xlabel = L"\sigma_{33} \; \text{ [MPa]}", ylabel = L"u_z \; \text{ [m]}", xlabelsize = 16, ylabelsize = 16, xticks=20:2:30)
    ax = Axis(f[1, 3], xlabel = L"p^\text{eff} = -\frac{\sigma_{kk}^\text{eff}}{3} \; \text{ [MPa]}", ylabel = L"u_z \; \text{ [m]}", xlabelsize = 16, ylabelsize = 16, xticks=20:2:30)
    #sc1 = scatter!(ax, s_33/1e6, u_z, markersize = 8, color=1:length(u_z), colormap=:rainbow, label=L"\max(u_z) \; \text{ (top Ranger)}")
    sc1 = scatter!(ax, p_eff/1e6, u_z, markersize = 8, color=1:length(u_z), colormap=:rainbow, label=L"\max(u_z) \text{, top Ranger}")
    if node_select == "max_subsidence"
        #xlims!(ax, -10, -20)
        xlims!(ax, 20, 30)
        ylims!(ax, -8, 0)
    elseif node_select == "max_rebound"
        xlims!(ax, -10, -20)
        ylims!(ax, -1, 5)
    end
elseif regime == "normal"
    ax = Axis(f[1, 3], xlabel = L"p^\text{eff} = -\frac{\sigma_{kk}^\text{eff}}{3} \; \text{ [MPa]}", ylabel = L"u_z \; \text{ [m]}", xlabelsize = 16, ylabelsize = 16, xticks=0:2:16, yticks=-10:2:0)
    sc1 = scatter!(ax, p_eff/1e6, u_z, markersize = 8, color=1:length(u_z), colormap=:rainbow, label=L"\max(u_z) \text{, top Ranger}")
    if node_select == "max_subsidence"
        xlims!(ax, 0, 16)
        ylims!(ax, -10, 0)
    elseif node_select == "max_rebound"
        xlims!(ax, 5, 15)
        ylims!(ax, -1, 5)
    end
elseif regime == "strike-slip"
    ax = Axis(f[1, 3], xlabel = L"p^\text{eff} = -\frac{\sigma_{kk}^\text{eff}}{3} \; \text{ [MPa]}", ylabel = L"u_z \; \text{ [m]}", xlabelsize = 16, ylabelsize = 16, xticks=10:2:24, yticks=-10:2:0)
    sc1 = scatter!(ax, p_eff/1e6, u_z, markersize = 8, color=1:length(u_z), colormap=:rainbow, label=L"\max(u_z) \text{, top Ranger}")
    if node_select == "max_subsidence"
        xlims!(ax, 10, 24)
        ylims!(ax, -10, 0)
    elseif node_select == "max_rebound"
        xlims!(ax, 5, 15)
        ylims!(ax, -1, 5)
    end
end
axislegend(position = :rb, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8), labelsize=14) 
f
#save(join([figs_dir, "stres_path_ranger_", node_select, "_", sim_name, ".pdf"]), f, pt_per_unit = 1)


# 
# Plot 2: max(u_z) vs t, comparison with subsidence data
#
data_u3max = [1947 1951 1954 1958 1965 1975; -2.4 -4 -6.1 -7.9 -8.8 -7.4]        # Kosloff et al., Tectonophysics (1980)
#data_u3max = [1975; 0.45]  
f = Figure(size=(400, 300))
#ax = Axis(f[1, 1], xlabel = L"t \; \text{ [y]}", ylabel = L"\max \; u_z \; \text{ [m]}", xlabelsize = 18, ylabelsize = 18, xticks=1935:5:1965, yticks =-10:2:0)
ax = Axis(f[1, 1], xlabel = L"t \; \text{ [y]}", ylabel = L"\max \; u_z \; \text{ [m]}", xlabelsize = 18, ylabelsize = 18, xticks=1935:5:1975, yticks =-10:2:0)
sc1 = scatter!(ax, data_u3max[1, :], data_u3max[2, :], markersize=10, color=RGBAf(1,1,1,0), marker=:rect, strokewidth=1, strokecolor=:black, label=L"\text{ Field data}")
ln1 = lines!(ax, t, u_z_top, color=:black, linewidth=1, label=L"\text{ Model}")
xlims!(ax, 1935, 1977)
ylims!(ax, -10, 0)
#ylims!(ax, 0, 10)
axislegend(position = :lb, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8)) 
f
#save(join([figs_dir, "u3_surface_plastic-89.pdf"]), f, pt_per_unit = 1)