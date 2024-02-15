# Code do compute pore pressure changes at nodal points (for Abaqus) based on cell-based pore pressure values from reservoir simulation (MRST)

# startup 
using MAT
using DelimitedFiles
using Revise

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
f2cn        = readdlm(join([dir_G, "f2cn.dat"]), '\t', Int, '\n')[:,1]
G           = matread(join([dir_G, fn_G]))["G"]
nodes_inp   = get_nodes_from_inp(dir_inp, fn_inp)
