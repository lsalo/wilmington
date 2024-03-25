"""

"""
# Cd into current directory and activate the Julia environment in this directory
cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__)

# startup
using LinearAlgebra
using Revise
using GLMakie
GLMakie.activate!()
using Colors
using ColorSchemes
using LaTeXStrings

lib_dir = "C:/Users/lsalo/julia/jl-tools/"
includet(joinpath(lib_dir, "modules", "gridOperations.jl"))
using .GridOperations
includet(joinpath(lib_dir, "modules", "plotFromMRST.jl"))
using .PlotFromMRST

# Directories
dir_inp = "C:/Users/lsalo/Wilmington/AbaqusGrids/"
dir_abaqus = "C:/Users/lsalo/abaqus/simulation_data/wilmington/output/elastic/extracted/S/"
figs_dir = "C:/Users/lsalo/julia/graphics/figs/wilmington/test/"  

# Get connectivity fron .inp grid
fn_inp  = "WilmingtonMesh_AllFaults_Rotated_ENTIRE_DOMAIN_Coarse.inp"
nodes_inp, elements_inp, elset_names_inp, nset_inp, nset_names_inp = get_connectivity_from_inp(dir_inp, fn_inp, "tetrahedral")

# Get nodeset faces - for now tetrahedra only
nsets_to_export = ["ComptonLosAlamitos_a", "ComptonLosAlamitos_b"]
face_list = [zeros(Int64, 0, 3) for _ in 1:length(nsets_to_export)]
for nset in eachindex(nsets_to_export)
    nset_id = findfirst(nset_names_inp .== nsets_to_export[nset])
    nset_list =  nset_inp[nset_id]

    for row in eachrow(elements_inp[1][:, 2:end])
        in_nset = in.(row, Ref(nset_list))
        if sum(in_nset) == 3                                             # element contains a face in the nodeset
            face_list[nset] = vcat(face_list[nset], sort(row[in_nset])')
        elseif sum(in_nset) > 3
            error("Multiple fault faces in a single element? Check")
        end
    end
    face_list[nset] = unique(face_list[nset], dims=1)                     # remove duplicate faces
end

# Face centroids and normals, and assign normal vectors to each node (where we export the stress from Abaqus)
# Get coordinates from nodes_inp, we'll have a loop for each face
unit_normals = [zeros(Float64, 0, 3) for _ in 1:length(nsets_to_export)]
for nset in eachindex(nsets_to_export)
    faces = face_list[nset]
    for face_node_ids in eachrow(faces)
        # Coords of face nodes
        coords = nodes_inp[in.(nodes_inp[:,1], Ref(face_node_ids)), 2:end]

        # unit upward-pointing normal from two edges
        f = coords[2,:]-coords[1,:]
        g = coords[3,:]-coords[1,:]
        n = cross(f, g) 
        for i in 3:-1:1 # ensure that we get upward normal vector, start with z, and normalize (unit vector)
            if n[i] != 0
                sgn = n[i] >= 0 ? 1 : -1
                n = n .* sgn / norm(n)  
                unit_normals[nset] = vcat(unit_normals[nset], n')
                break
            end
        end

    end
end

# Calculate shear, normal tractions at each node
# Export a file with node, ts, tn, 

# Test plot
# Better like this than using a loop to have unique 1 to n numbering and plot without lines showing fault segment boundaries
coords_map = Vector{}(undef, 1)
vals_nodes = Vector{}(undef, 1)
face_list_map = Vector{}(undef, 1)
unique_ids = sort(unique(vcat(face_list...)))
coords_map[1] = nodes_inp[in.(nodes_inp[:, 1], Ref(unique_ids)), 2:end]
vals_nodes[1] = coords[:, 3] ./ 1e3
face_list_map[1], _ = renumber_array_1_to_n(vcat(face_list...))

mycmap = reverse(cgrad(:turbid))
fig_descriptor = (cmap=mycmap, cbar_label=L"$z$ [km]", fig_size = (800, 600), plot_height = 900, cbar_height = 200, plot_azi = 1.8*pi, plot_elev = 0.3*pi, plot_linew=0.1, xlms=(nothing, nothing), ylms=(nothing, nothing), save_fig = false, save_dir = figs_dir, title = "Depth", target_dpi=350, target_size_in=(8,6))
plotCellValues(face_list_map, coords_map, vals_nodes, fig_descriptor)

# test tractions
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