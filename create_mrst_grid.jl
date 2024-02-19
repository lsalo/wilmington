# startup 
using MAT
using DelimitedFiles
using Revise
using StatsBase

lib_dir = "C:/Users/lsalo/julia/jl-tools/"
includet(joinpath(lib_dir, "modules", "gridOperations.jl"))
using .GridOperations

# Directories
dir_inp     = "C:/Users/lsalo/Wilmington/AbaqusGrids/"
fn_inp      = "WilmingtonMesh_AllFaults_Rotated_Coarse.inp"
dir_grid    = "C:/Users/lsalo/matlab/mrst-dev/mrst-ls/ls-proj/wilmington/grid/"
fn_grid     = "wilmington_ranger_rotated_tetrahedral"
fn_G        = "wilmington_ranger_G.mat"

# Get connectivity
nodes_inp, elements_inp, elset_names_inp, nset_inp, nset_names_inp = get_connectivity_from_inp(dir_inp, fn_inp, "tetrahedral")

# Write to text files to generate MRST grid
elsets_to_export = ["Ranger"] # "all" or vector of strings 
nsets_to_export = ["Tear_Fault_2", "Tear_Fault_4", "Tear_Fault_5", "Tear_Fault_6", "Tear_Fault_9", "Tear_Fault_11", "Tear_Fault_15", "Tear_Fault_19", "Tear_Fault_26", "Tear_Fault_27", "Tear_Fault_8a", "Tear_Fault_10a", "Tear_Fault_13a", "Tear_Fault_21a", "Tear_Fault_8b", "Tear_Fault_10b", "Tear_Fault_13b", "Tear_Fault_21b"]
nset_info = [nset_inp, nset_names_inp, nsets_to_export]
nset_mrst_ids = writeTetrahedralGridForMRST(nodes_inp, elements_inp, elset_names_inp, elsets_to_export, nset_info, dir_grid, fn_grid)

# Find fault faces for transmissibility multipliers
G = matread(join([dir_G, fn_G]))["G"]
faces_nodes, nset_mrst_face_ids = findFaultFacesMRSTGrid(G, nset_mrst_ids, nsets_to_export, dir_grid)

# Computation of transmissibility multipliers
# TBD if needed. For now, we just want 1s and 0s (for impermeable fault faces) and assign directly in MRST.
