#=
Script to extract nodeset info from .inp file
=#

# Packages
cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__)

using DelimitedFiles
using Printf

using Revise
lib_dir = "C:/Users/lsalo/julia/jl-tools/"
includet(joinpath(lib_dir, "modules", "arrayOps.jl"))
using .ArrayOps
includet(joinpath(lib_dir, "modules", "abaqus.jl"))
using .Abaqus

# Directories and inputs
dir_inp = "C:/Users/lsalo/Wilmington/AbaqusGrids/"
fn_inp  = "WilmingtonMesh_AllFaults_Coarse.inp"
dir_out = "C:/Users/lsalo/julia/output/wilmington/topRanger/"

# Get info
nsets_to_export = ["TopRangerSurface"]
nodes_inp, elements_inp, elset_names_inp, nset_inp, nset_names_inp = get_connectivity_from_inp(dir_inp, fn_inp, "tetrahedral")
face_list = get_face_connectivity_nodesets(nsets_to_export, nset_names_inp, nset_inp, elements_inp[1], "tetrahedral")
nodes_coords_topRanger = nodes_inp[in.(nodes_inp[:, 1], Ref(nset_inp[2])), :]

# Format and write to file
formatted_nodes_inp = [@sprintf("%10d\t%10.1f\t%10d\t%12.4f", row[1], row[2], row[3], row[4]) for row in eachrow(nodes_coords_topRanger)]
formatted_face_list = [@sprintf("%10d\t%10d\t%10d\t", row[1], row[2], row[3]) for row in eachrow(face_list[1])]

#write(joinpath(dir_out, "nodes_coordinates_topRanger.dat"), join(formatted_nodes_inp, "\n"))
#write(joinpath(dir_out, "nodes_topRanger.dat"), join(nset_inp[2], "\n"))
#write(joinpath(dir_out, "facelist_topRanger.dat"), join(formatted_face_list, "\n"))


