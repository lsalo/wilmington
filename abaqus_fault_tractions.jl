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

#lib_dir = "C:/Users/lsalo/julia/jl-tools/"
lib_dir = "/Users/lluis/Documents/julia/jl-tools/"
includet(joinpath(lib_dir, "modules", "faultOperations.jl"))
using .FaultOperations
includet(joinpath(lib_dir, "modules", "gridOperations.jl"))
using .GridOperations
includet(joinpath(lib_dir, "modules", "plotFromMRST.jl"))
using .PlotFromMRST

# Directories
# dir_inp = "C:/Users/lsalo/Wilmington/AbaqusGrids/"
# dir_abaqus = "C:/Users/lsalo/abaqus/simulation_data/wilmington/output/elastic/extracted/S/"
# figs_dir = "C:/Users/lsalo/julia/graphics/figs/wilmington/test/" 
dir_inp = "/Users/lluis/Dropbox (MIT)/home_copy/abaqus/wilmington/grids/"
dir_abaqus = "/Users/lluis/Dropbox (MIT)/home_copy/abaqus/wilmington/output/elastic/extracted/S/"
figs_dir = "/Users/lluis/Documents/julia/graphics/figs/wilmington/" 

# Get connectivity fron .inp grid
fn_inp  = "WilmingtonMesh_AllFaults_Rotated_ENTIRE_DOMAIN_Coarse.inp"
nodes_inp, elements_inp, elset_names_inp, nset_inp, nset_names_inp = get_connectivity_from_inp(dir_inp, fn_inp, "tetrahedral")

# Get nodeset faces - for now tetrahedra only
nsets_to_export = ["ComptonLosAlamitos_a", "ComptonLosAlamitos_b"]
face_list = get_face_connectivity_nodesets(nsets_to_export, nset_names_inp, nset_inp, elements_inp[1], "tetrahedral")

# Face centroids and normals
face_unit_normals, face_centroids = face_normals_and_centroids(face_list, nodes_inp)

# Assign normal vectors to each node (where we export the stress from Abaqus)
# Calculate shear, normal tractions at each node
repeated_nodes = true
nodal_tractions = [zeros(Union{Int64, Float64}, 0, 3) for _ in 1:length(nsets_to_export)]
for nset in eachindex(nsets_to_export)
    tstep = "Geostatic"
    nodal_stress = readdlm(join([dir_abaqus, nsets_to_export[nset], "/", tstep, ".dat"]), ',', Int64, '\n', skipstart=1)
    fname = join([dir_abaqus, nsets_to_export[nset], "/tractions_", tstep, ".dat"])
    nodal_tractions[nset] = tractions_from_sij_and_n(nodal_stress, repeated_nodes, face_list[nset], face_unit_normals[nset], face_centroids[nset], nodes_inp, fname)   
end


# Test plot
# Better like this than using a loop: unique 1 to n numbering and plot without lines showing fault segment boundaries
# we assume nodes_inp to be monotonically increasing
unique_ids = sort(unique(vcat(face_list...)))   # nodes may still be repeated when >1 nodesets share nodes
coords_map = Vector{}(undef, 1)
coords_map[1] = nodes_inp[in.(nodes_inp[:, 1], Ref(unique_ids)), 2:end]

# vals nodes
vals_nodes = Vector{}(undef, 1)
tmp = sortslices(vcat(nodal_tractions...), dims=1, by=x->x[1])
nodal_tractions_plot = zeros(Float64, length(unique_ids), 2)
for (i, node) in enumerate(unique_ids) 
    id_rows = tmp[:,1] .== node
    nodal_tractions_plot[i, :] = mean(tmp[id_rows, 2:end], dims=1)
end
vals_nodes[1] = nodal_tractions_plot[:,1] ./ 1e6 # MPa

# face connectivity
face_list_map = Vector{}(undef, 1)
face_list_map[1], _ = renumber_array_1_to_n(vcat(face_list...))

mycmap = cgrad(:nipy_spectral)
fig_descriptor = (cmap=mycmap, cbar_label=L"$t_\text{n}$ [MPa]", fig_size = (800, 600), plot_height = 900, cbar_height = 200, plot_azi = 1.8*pi, plot_elev = 0.3*pi, plot_linew=0.1, xlms=(nothing, nothing), ylms=(nothing, nothing), save_fig = false, save_dir = figs_dir, title = "t_n", target_dpi=350, target_size_in=(8,6))
plotCellValues(face_list_map, coords_map, vals_nodes, fig_descriptor)

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
