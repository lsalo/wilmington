# Code do compute pore pressure changes at nodal points (for Abaqus) based on cell-based pore pressure values from reservoir simulation (MRST)

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
fn_G    = "wilmington_G.mat"
dir_p   = "C:/Users/lsalo/matlab/mrst-dev/mrst-ls/ls-proj/wilmington/results/"
fn_p    = "pressures.mat"

# Import
f2cn        = readdlm(join([dir_G, "f2cn.dat"]), '\t', Int, '\n')[:,1]  # from MRST, see Grid structure in textbook
n2fn        = readdlm(join([dir_G, "n2fn.dat"]), '\t', Int, '\n')[:,1] # "
G           = matread(join([dir_G, fn_G]))["G"]
nodes_inp, elements_inp, elset_names_inp = get_nodes_from_inp(dir_inp, fn_inp, "tetrahedral")

# Node-cell connections map
nodes_to_cells = map_nodes_to_cells(G, n2fn, f2cn)

# Get nodal pressures and compute Dp
# This needs to be done at each timestep, and then compute Dp
p_test  = rand(1:10^6, maximum(f2cn), 1) 
nodal_p = Dict((node_id, mean(p_test[cell_ids])) for (node_id, cell_ids) in nodes_to_cells)
# TBD: Compute Dp

# Map MRST nodes to global (.inp) nodes using nodal coordinates
nodes_MRST = Matrix{Union{Int, Float64}}(undef, size(nodes_inp, 1), 4)
nodes_MRST[:,1] = 1:size(nodes_inp, 1)
nodes_MRST[:,2:4] = G["nodes"]["coords"] # Attention: Here in the G structure we have ALL nodes, not just those of the cells used in reservoir simulation
nodes_MRST[:,4] = -1*nodes_MRST[:,4]  #
#for node in 1:size(nodes_MRST, 1)
#    coords = nodes_MRST[node,2:4]
#    findall(coords' .== nodes_inp[:,2:4]) # MRST coordinates don't overlap...TBD
#end

# Assign nodal Dp to global node ids (focus on reservoir cells only)
