#=
Script to compute fault tractions at wilmington
=#

# Cd into current directory and activate the Julia environment in this directory
cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__)

# startup
using Revise
using GLMakie
GLMakie.activate!()
using Colors
using ColorSchemes
using LaTeXStrings
using DelimitedFiles
using Statistics
using LinearAlgebra

lib_dir = "C:/Users/lsalo/julia/jl-tools/"
#lib_dir = "/Users/lluis/Documents/julia/jl-tools/"
includet(joinpath(lib_dir, "modules", "faultOperations.jl"))
using .FaultOperations
includet(joinpath(lib_dir, "modules", "abaqus.jl"))
using .Abaqus
includet(joinpath(lib_dir, "modules", "arrayOps.jl"))
using .ArrayOps
includet(joinpath(lib_dir, "modules", "fileio.jl"))
using .Fileio
includet(joinpath(lib_dir, "modules", "makiePolyPlot.jl"))
using .MakiePolyPlot

# Directories
dir_inp = "C:/Users/lsalo/Wilmington/AbaqusGrids/"
dir_abaqus = "C:/Users/lsalo/abaqus/simulation_data/wilmington/plastic-88/extracted/S/"
dir_out = "C:/Users/lsalo/julia/output/wilmington/fault_quantities/plastic-88/"
figs_dir = "C:/Users/lsalo/julia/graphics/figs/wilmington/test/"
#dir_inp = "/Users/lluis/Dropbox (MIT)/home_copy/abaqus/wilmington/grids/"
#dir_abaqus = "/Users/lluis/Dropbox (MIT)/home_copy/abaqus/wilmington/output/elastic/extracted/S/"
#figs_dir = "/Users/lluis/Documents/julia/graphics/figs/wilmington/" 

# Get connectivity fron .inp grid
fn_inp  = "WilmingtonMesh_AllFaults_Rotated_ENTIRE_DOMAIN_Coarse.inp"
nodes_inp, elements_inp, elset_names_inp, nset_inp, nset_names_inp = get_connectivity_from_inp(dir_inp, fn_inp, "tetrahedral")

# Get nodeset faces - for now tetrahedra only
#nsets_to_export = ["Tear_Fault_2", "Tear_Fault_4", "Tear_Fault_5", "Tear_Fault_6", "Tear_Fault_8a", "Tear_Fault_8b", "Tear_Fault_9", "Tear_Fault_10a", "Tear_Fault_10b", "Tear_Fault_11", "Tear_Fault_13a", "Tear_Fault_13b", "Tear_Fault_15", "Tear_Fault_19", "Tear_Fault_21a", "Tear_Fault_21b", "Tear_Fault_26", "Tear_Fault_27"]
nsets_to_export = ["Wilmington"]
face_list = get_face_connectivity_nodesets(nsets_to_export, nset_names_inp, nset_inp, elements_inp[1], "tetrahedral")

# Face centroids and normals
face_unit_normals, face_centroids = face_normals_and_centroids(face_list, nodes_inp)

# Calculate shear, normal tractions at each node, and dcff
mu = 0.6
nodal_tractions = [Vector{Any}() for _ in 1:length(nsets_to_export)]
#dcff = Vector{Matrix{Union{Int64,Float64}}}(undef, length(nsets_to_export))
dcff = [Vector{Any}() for _ in 1:length(nsets_to_export)]
compute_tractions = false
for nset in eachindex(nsets_to_export)
    abaqus_fn_nod = readdir(join([dir_abaqus, nsets_to_export[nset]]))
    abaqus_fn_nod = sort(abaqus_fn_nod, by=sort_steps)
    for (i, fn) in enumerate(abaqus_fn_nod)
        if compute_tractions
            nodal_stress = readdlm(join([dir_abaqus, nsets_to_export[nset], "/", fn]), ',', Int64, '\n', skipstart=1)
            dir_path = joinpath([dir_out, nsets_to_export[nset]])
            fname = join(["tractions_", fn])
            push!(nodal_tractions[nset], tractions_from_sij_and_n(nodal_stress, nodes_inp, face_unit_normals[nset], face_centroids[nset], dir_path, fname))           

            print(join([fn " completed.\n"]))

        else # we already have them, so we just read
            data = readdlm(joinpath([dir_out, nsets_to_export[nset], join(["tractions_", fn])]), '\t', Float64, '\n')
            push!(nodal_tractions[nset], data)
        end

        if i > 1
            dtau = nodal_tractions[nset][i][:, end] .- nodal_tractions[nset][1][:, end]
            dsn = nodal_tractions[nset][i][:, 2] .- nodal_tractions[nset][1][:, 2]
            cff0 = nodal_tractions[nset][1][:, end] .+ mu*nodal_tractions[nset][1][:, 2]
            cff = nodal_tractions[nset][i][:, end] .+ mu*nodal_tractions[nset][i][:, 2] 
            push!(dcff[nset], hcat(nodal_tractions[nset][1][:,1], cff .- cff0, dsn, dtau))
        end
    end
    
    print(join(["Fault tractions for nodeset ", nsets_to_export[nset], " completed.\n"]))
end

# Averages for each fault
n_steps = length(dcff[1])
dcff_avg = [Vector{Float64}(undef, n_steps) for _ in 1:length(nsets_to_export)]
dcff_10 = [Vector{Float64}(undef, n_steps) for _ in 1:length(nsets_to_export)]
dcff_90 = [Vector{Float64}(undef, n_steps) for _ in 1:length(nsets_to_export)]
for nset in eachindex(nsets_to_export)
    for i in 1:n_steps
        dcff_avg[nset][i] = mean(dcff[nset][i][:, 2])
        dcff_10[nset][i] = quantile(dcff[nset][i][:, 2], 0.1)
        dcff_90[nset][i] = quantile(dcff[nset][i][:, 2], 0.9)
    end
end

# Get intersection of TopRangerSurface and BottomRangerSurface with the fault
nsets_to_intersect = ["TopRangerSurface", "BottomRangerSurface", "TopBasementSurface"]
face_list_intersect = get_face_connectivity_nodesets(nsets_to_intersect, nset_names_inp, nset_inp, elements_inp[1], "tetrahedral")
nodes_coords_plot = Vector{}(undef, length(nsets_to_intersect))
for nset in eachindex(nsets_to_intersect)
    unique_nodes = unique(face_list_intersect[nset])
    common_nodes = intersect(unique_nodes, face_list[1]) # update index of face_list according to nset for intersections
    nodes_coords_plot[nset]  = vcat([nodes_inp[nodes_inp[:, 1] .== node_id, 2:4] for node_id in common_nodes]...)
end


# Plots
# 1. Plot average dcff on each fault vs t
using CairoMakie
CairoMakie.activate!()

t = collect(1937:1961)  # Make Auto 
#cmap = cgrad(:nipy_spectral, length(nsets_to_export))
cmap = cgrad(:roma, length(nsets_to_export))

f = Figure(size=(900, 350))
ax = Axis(f[1, 1], xlabel = L"t \; \text{ [y]}", ylabel = L"\overline{\text{DCFF}} \; \text{ [bar]}", xlabelsize = 16, ylabelsize = 16)
for (i, data) in enumerate(dcff_avg)
    normalized_i = i / length(nsets_to_export)
    lbl = replace(nsets_to_export[i], "Tear_Fault_" => "TF")
    lbl = latexstring("\\text{" * lbl * "}")
    #lbl = latexstring("\\text{" * nsets_to_export[i] * "}")
    lines!(ax, t, data./1e5, color = cmap[normalized_i], linewidth = 1.5, label=lbl)
end
xlims!(ax, 1936, 1960)
ylims!(ax, -10, 1)
axislegend(position = :lb, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8), nbanks=2) 

#ax = Axis(f[1, 2], ylabel = L"\min{\text{ DCFF}} \; \text{ [bar]}", xlabelsize = 16, ylabelsize = 16)
ax = Axis(f[1, 2], ylabel = L"\text{P}{10} \text{ DCFF} \; \text{ [bar]}", xlabelsize = 16, ylabelsize = 16)
for (i, data) in enumerate(dcff_10)
    normalized_i = i / length(nsets_to_export)
    lines!(ax, t, data./1e5, color = cmap[normalized_i], linewidth = 1.5)
end
xlims!(ax, 1936, 1960)
ylims!(ax, -25, 1)

#ax = Axis(f[1, 3], ylabel = L"\max{\text{ DCFF}} \; \text{ [bar]}", xlabelsize = 16, ylabelsize = 16)
ax = Axis(f[1, 3], ylabel = L"\text{P}{90} \text{ DCFF} \; \text{ [bar]}", xlabelsize = 16, ylabelsize = 16)
for (i, data) in enumerate(dcff_90)
    normalized_i = i / length(nsets_to_export)
    lines!(ax, t, data./1e5, color = cmap[normalized_i], linewidth = 1.5)
end
xlims!(ax, 1936, 1960)
ylims!(ax, -1, 3)
f
save(join([figs_dir, "dcff_vs_t_TearFaults.pdf"]), f, pt_per_unit = 1)


# 2. Plot Map
# Better like this than using a loop: unique 1 to n numbering and plot without lines showing fault segment boundaries
# we assume nodes_inp to be monotonically increasing
using GLMakie
GLMakie.activate!()

unique_ids = sort(unique(vcat(face_list...)))   # nodes may still be repeated when >1 nodesets share nodes
coords_map = Vector{}(undef, 1)
coords_map[1] = nodes_inp[in.(nodes_inp[:, 1], Ref(unique_ids)), 2:end]

# vals nodes
#tstep = "step_264_nodes.dat" # t=1958, before injectors turn on
tstep = "step_468_nodes.dat" # t=1975, end of repressurization
abaqus_fn_nod = readdir(join([dir_abaqus, nsets_to_export[1]]))
abaqus_fn_nod = sort(abaqus_fn_nod, by=sort_steps)
idx = findfirst(==(tstep), abaqus_fn_nod)
dcff_plot = Vector{Matrix{Union{Int64,Float64}}}(undef, length(nsets_to_export))
for nset in eachindex(nsets_to_export)
    dcff_plot[nset] = dcff[nset][idx-1]
end
tmp = sortslices(vcat(dcff_plot...), dims=1, by=x->x[1])
values_plot = zeros(Float64, length(unique_ids), 3)
for (i, node) in enumerate(unique_ids) 
    id_rows = tmp[:,1] .== node
    values_plot[i, :] = mean(tmp[id_rows, 2:end], dims=1)
end
vals_nodes = Vector{}(undef, 1)
vals_nodes[1] = values_plot[:, 3] ./ 1e5         # bar, column index 1 = dcff, 2 = dsn, 3 = dtau
#labl = L"$\sigma_\text{n}^\prime$ [MPa]"       # "$|\tau|$ [MPa]"
#labl = L"\text{DCFF [bar]}"
#labl = L"\Delta \sigma_\text{n}^{\text{eff}} \; \text{ [bar]}"
labl = L"\Delta |\tau| \; \text{ [bar]}"

# face connectivity
face_list_map = Vector{}(undef, 1)
face_list_map[1], _ = renumber_array_1_to_n(vcat(face_list...))

#mycmap = cgrad(:oleron)
#mycmap = cgrad(:lajolla)
mycmap = cgrad(:curl)
fig_descriptor = (cmap=mycmap, cbar_label=labl, fig_size = (1200, 900), plot_height = 700, cbar_height = 300, plot_azi = -0.2*pi, plot_elev = 0.3*pi, plot_linew=0.1, xlms=(nothing, nothing), ylms=(nothing, nothing), save_fig = true, save_dir = figs_dir, title = "wilmington_dtau_plastic-88_tstep468", target_dpi=400, target_size_in=(12,9))
plotNodalValues(face_list_map, coords_map, vals_nodes, fig_descriptor, nodes_coords_plot)


# test tractions
#=
S = [15 0 0; 0 12 0; 0 0 10]
coord = [2 2 1; 6 3 0; 4 1 1]
f = coord[2,:]-coord[1,:]
g = coord[3,:]-coord[1,:]
n = cross(f, g) 
for i in 3:-1:1 # start with z component, we want the outward normal vector
    if n[i] != 0
        sgn = n[i] >= 0 ? 1 : -1
        n = n .* sgn / norm(n)
        break
    end
end

t = S*n
tn = dot(n, t)
ts = sqrt(norm(t)^2 - tn^2)
=#
