#=
Script to compute fault tractions at wilmington
=#

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
using DelimitedFiles
using Statistics

lib_dir = "C:/Users/lsalo/julia/jl-tools/"
includet(joinpath(lib_dir, "modules", "faultOperations.jl"))
using .FaultOperations
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
face_list = get_face_connectivity_nodesets(nsets_to_export, nset_names_inp, nset_inp, elements_inp[1], "tetrahedral")

# Face centroids and normals
face_unit_normals = [zeros(Float64, 0, 3) for _ in 1:length(nsets_to_export)]
face_centroids = [zeros(Float64, 0, 3) for _ in 1:length(nsets_to_export)]
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
                face_unit_normals[nset] = vcat(face_unit_normals[nset], n')
                break
            end
        end

        # centroids
        x = (coords[1,1] + coords[2,1] + coords[3,1]) / 3
        y = (coords[1,2] + coords[2,2] + coords[3,2]) / 3 
        z = (coords[1,3] + coords[2,3] + coords[3,3]) / 3
        face_centroids[nset] = vcat(face_centroids[nset], [x, y, z]')   
    end

end

# Assign normal vectors to each node (where we export the stress from Abaqus)
# Calculate shear, normal tractions at each node
repeated_nodes = true
nodal_tractions = [zeros(Union{Int64, Float64}, 0, 3) for _ in 1:length(nsets_to_export)]
for nset in eachindex(nsets_to_export)
    tstep = "Geostatic"
    nodal_stress = readdlm(join([dir_abaqus, nsets_to_export[nset], "/", tstep, ".dat"]), ',', Int64, '\n', skipstart=1)
    nodal_stress = sortslices(nodal_stress, dims=1, by=x->x[1])
    
    # Assign normals to nodes based on distance of node to face centroids
    # Average nodal stresses in case of repeated nodes (extrapolation of stress to the nodes)
    unique_ids = sort(unique(face_list[nset]))
    tmp = zeros(Int, length(unique_ids), 7)
    node_unit_normals = zeros(Union{Int64, Float64}, length(unique_ids), 4)
    for (i, node) in enumerate(unique_ids)
        if repeated_nodes
            tmp[i, 1] = node
            id_rows = nodal_stress[:,1] .== node
            tmp[i, 2:end] = Int.(round.(mean(nodal_stress[id_rows, 2:end], dims=1)))
        end

        node_unit_normals[i, 1] = node
        node_coord = nodes_inp[nodes_inp[:,1] .== node, 2:end]
        dist = sqrt.(sum((face_centroids[nset] .- node_coord) .^2, dims=2))
        node_unit_normals[i, 2:end] = face_unit_normals[nset][argmin(dist)[1], :]
    end
    if repeated_nodes
        nodal_stress = tmp
    end

    # tractions
    n = node_unit_normals[:, 2:end]
    s = nodal_stress[:, 2:end]
    sxx = s[:,1]
    syy = s[:,2]
    szz = s[:,3]
    sxy = s[:,4]
    sxz = s[:,5]
    syz = s[:,6]

    tx = sxx.*n[:,1] + sxy.*n[:,2] + sxz.*n[:,3]
    ty = sxy.*n[:,1] + syy.*n[:,2] + syz.*n[:,3]
    tz = sxz.*n[:,1] + syz.*n[:,2] + szz.*n[:,3]

    t = hcat(tx, ty, tz)
    tn = dot.(eachrow(n), eachrow(t))
    ts = sqrt.(norm.(eachrow(t)).^2 - tn.^2)
    
    tmp = Matrix{Union{Int64, Float64}}(undef, length(unique_ids), 3)
    tmp[:, 1] = unique_ids
    tmp[:, 2:end] = hcat(tn, ts)
    nodal_tractions[nset] = tmp

    # Export a file with node, ts, tn, 
    open(join([dir_abaqus, nsets_to_export[nset], "/tractions_", tstep, ".dat"]), "w") do f
        for row in eachrow(nodal_tractions[nset])
            println(f, join(row, "\t"))
        end
    end
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
