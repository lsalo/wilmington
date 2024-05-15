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
includet(joinpath(lib_dir, "modules", "abaqus.jl"))
using .Abaqus
includet(joinpath(lib_dir, "modules", "fileio.jl"))
using .Fileio
includet(joinpath(lib_dir, "modules", "arrayOps.jl"))
using .ArrayOps


# Directories and options
abaqus_dir_el = "C:/Users/lsalo/abaqus/simulation_data/wilmington/output/plastic/extracted/S/ENTIRE_DOMAIN/"
abaqus_fn_el = readdir(abaqus_dir_el)
abaqus_fn_el = sort(abaqus_fn_el, by=sort_steps)

abaqus_dir_nod = "C:/Users/lsalo/abaqus/simulation_data/wilmington/output/plastic/extracted/U/TOPRANGERSURFACE/"
abaqus_fn_nod = readdir(abaqus_dir_nod)
abaqus_fn_nod = sort(abaqus_fn_nod, by=sort_steps)

abaqus_dir_nod2 = "C:/Users/lsalo/abaqus/simulation_data/wilmington/output/plastic/extracted/S/TOPRANGERSURFACE/"
abaqus_fn_nod2 = readdir(abaqus_dir_nod2)
abaqus_fn_nod2 = sort(abaqus_fn_nod2, by=sort_steps)

dir_inp = "C:/Users/lsalo/Wilmington/AbaqusGrids/"
fn_inp  = "WilmingtonMesh_AllFaults_Rotated_Coarse.inp"

dir_in = "C:/Users/lsalo/julia/output/wilmington/"
figs_dir = "C:/Users/lsalo/julia/graphics/figs/wilmington/test/" 


# Get node of max subsidence, and uz vs s33 evolution for that node
u_topranger_fin = readdlm(join([abaqus_dir_nod, abaqus_fn_nod[end]]), ',', skipstart=1)
id_umax = argmin(u_topranger_fin[:, end])
nod_umax = Int(u_topranger_fin[id_umax, 1])

n = length(abaqus_fn_nod)
u_z = Vector{Float32}(undef, n)
for i in 1:n
    file = readdlm(join([abaqus_dir_nod, abaqus_fn_nod[i]]), ',', skipstart=1)
    u_z[i] = file[id_umax, end]
end

s_33 = Vector{Int32}(undef, n)
for i in 1:n
    file = h5open(join([abaqus_dir_nod2, abaqus_fn_nod2[i]]), "r")
    abq_out = [read(file[name]) for name in ["NODE_LABEL", "Sxx", "Syy", "Szz", "Sxy", "Sxz", "Syz"]]
    abq_out = Int.(hcat(abq_out...))
    s_nod = reduce_matrix_with_repeated_ids(abq_out, "mean") # currently very slow, needs to be revised / run in parallel
    s_33[i] = s_nod[id_umax, 4]
    print("$i/$n completed\n")
end

# Get Ranger element with closest centroid to node above
# Note that we assume that Ranger element ids start at 1, and data are sorted
nodes_inp, elements_inp, elset_names_inp, nset_inp, nset_names_inp = get_connectivity_from_inp(dir_inp, fn_inp, "tetrahedral")
id_elements = elements_inp[1][any(elements_inp[1][:, 2:5] .== nod_umax, dims=2), 1]

centroids = readdlm(join([dir_in, "element_centroids.dat"]), ',')
centroids_els = centroids[id_elements, 2:end]
coord_nod_umax = nodes_inp[findfirst(nodes_inp[:, 1] .== nod_umax), 2:end]
dist = mapslices(row -> norm(row - coord_nod_umax), centroids_els, dims=2)
id_element = id_elements[argmin(dist)]


# Get geostatic and final stresses
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

# Plot parameters
# Cap model parameters
d = 1e-3        # MPa
b = 52.41       # deg
R = 0.2         # 1.227
W = 0.27        # [-]
D = 0.7143      # 1/kbar
e_p = [0, 0.0025, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2]

# Hardening curve [Pa]
# p_eff_el below should be very close to p_eff[1] above
s_ini_el = readdlm(join([dir_in, "initial_stress_rotated_criticallyStressed_mu0.6.inp"]), ',')
p_eff_ini_el = -(1/3)*sum(s_ini_el[id_element, 2:4])
p_b = 1.18*p_eff_ini_el .+ 0.25*((1/3)*(1/D)*log.(e_p./W .+ 1)*1e8) 

# Shear envelope and cap
e_p0 = 0
p_Fs = Vector{}(0:10:50)
q_Fs = p_Fs.*tand(b) .+ d

p_a = (1e-6*p_b[1] - R*d) / (1 + R*tand(b))
q_Fc = vcat(Vector{}(0:1:28))
p_Fc = p_a .+ sqrt.( abs.((R*(d + p_a*tand(b)))^2 .- (R.*q_Fc).^2) )


# Make plot stress path
f = Figure(size=(800, 400))
ax = Axis(f[1, 1], xlabel = L"p^\text{eff} = -\frac{\sigma_{kk}^\text{eff}}{3} \; \text{ [MPa]}", ylabel = L"q = \sqrt{3 J_2} \; \text{ [MPa]}", xlabelsize = 16, ylabelsize = 16)
ln1 = lines!(ax, p_Fs, q_Fs, color = RGBAf(0, 0, 0, 1), linewidth = 1.5, label=L"F_\text{s}")
ln1 = lines!(ax, p_Fc, q_Fc, color = RGBAf(0.5, 0.5, 0.5, 1), linewidth = 1, linestyle = :dashdot, label=L"F_\text{c}")
sc1 = scatter!(ax, p_eff/1e6, q/1e6, markersize = 8, color=1:length(p_eff), colormap=:viridis, label=L"\text{Ranger element}")
xlims!(ax, 0, 50)
ylims!(ax, 0, 50)
axislegend(position = :lt, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8)) 

ax = Axis(f[1, 2], xlabel = L"e^\text{in}_\text{v} \; \text{ [-]}", ylabel = L"p_\text{b} \; \text{ [MPa]}", xlabelsize = 16, ylabelsize = 16)
ln1 = lines!(ax, e_p, p_b*1e-6, color = RGBAf(0, 0, 0, 1), linewidth = 1.5)
xlims!(ax, 0, 0.2)
ylims!(ax, 0, 75)
f


# Make plot displacement
f = Figure(size=(400, 400))
ax = Axis(f[1, 1], xlabel = L"\sigma_{33} \; \text{ [MPa]}", ylabel = L"u_3 \; \text{ [m]}", xlabelsize = 16, ylabelsize = 16)
sc1 = scatter!(ax, s_33/1e6, u_z, markersize = 8, color=:black, label=L"\text{Top Ranger } \max(u_3)")
xlims!(ax, -10, -20)
f