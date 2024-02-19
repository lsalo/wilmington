# startup 
using MAT
using DelimitedFiles
using Revise
using StatsBase

lib_dir = "C:/Users/lsalo/julia/jl-tools/"
includet(joinpath(lib_dir, "modules", "gridOperations.jl"))
using .GridOperations

# Directories
dir_inp = "C:/Users/lsalo/Wilmington/AbaqusGrids/"
fn_inp  = "WilmingtonMesh_AllFaults_Rotated_Coarse.inp"
dir_G   = "C:/Users/lsalo/matlab/mrst-dev/mrst-ls/ls-proj/wilmington/grid/"
fn_G    = "wilmington_ranger_rotated_tetrahedral"

# Get connectivity
nodes_inp, elements_inp, elset_names_inp, nset_inp, nset_names_inp = get_connectivity_from_inp(dir_inp, fn_inp, "tetrahedral")

# Write to text files to generate MRST grid
elsets_to_export = ["Ranger"] # "all" or vector of strings 
nsets_to_export = ["Tear_Fault_2"]
nset_info = [nset_inp, nset_names_inp, nsets_to_export]
writeTetrahedralGridForMRST(nodes_inp, elements_inp, elset_names_inp, elsets_to_export, nset_info, dir_G, fn_G)

# Find fault faces for transmissibility multipliers

