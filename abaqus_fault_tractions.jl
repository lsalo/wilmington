"""

"""

# startup
using LinearAlgebra
using Revise

lib_dir = "C:/Users/lsalo/julia/jl-tools/"
includet(joinpath(lib_dir, "modules", "gridOperations.jl"))
using .GridOperations

# Directories
dir_inp = "C:/Users/lsalo/Wilmington/AbaqusGrids/"
dir_abaqus = "C:/Users/lsalo/abaqus/simulation_data/wilmington/output/elastic/extracted/S/"

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

# Face centroids and normals, and assign normal vectors to each (where we obtain the stress from Abaqus)
# Get coordinates from nodes_inp, we'll have a loop for each face
# First: PLOT FACES with dummy values to ensure faces are correct.
face_node_ids = face_list[nset][face, :]
face_node_coords = nodes_inp[in.(nodes_inp[:,1], Ref(face_node_ids)), :]

# Calculate shear, normal tractions at each node
# Export a file with node, ts, tn, 

# test
S = [15 0 0; 0 12 0; 0 0 10]
coord = [2 2 1; 6 3 0; 4 1 1]
f = coord[2,:]-coord[1,:]
g = coord[3,:]-coord[1,:]
n = cross(f, g) 
for i in 3:-1:1 # start with z component, we want the normal vector positive upward
    if n[i] != 0
        sgn = n[i] >= 0 ? 1 : -1
        n = n .* sgn / norm(n)
        break
    end
end

t = S*n
tn = dot(n, t)
ts = sqrt(norm(t)^2 - tn^2)
