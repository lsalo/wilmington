# Code do compute pore pressure changes at nodal points (for Abaqus) based on cell-based pore pressure values from reservoir simulation (MRST)

# startup 
using MAT
using DelimitedFiles

# Directories
dir_inp = "C:/Users/lsalo/Wilmington/AbaqusGrids/"
fn_inp  = "WilmingtonMesh_AllFaults_Rotated_Coarse.inp"
dir_G   = "C:/Users/lsalo/matlab/mrst-dev/mrst-ls/ls-proj/wilmington/grid/"
fn_G    = "wilmington_G.mat"
dir_p   = "C:/Users/lsalo/matlab/mrst-dev/mrst-ls/ls-proj/wilmington/results/"
fn_p    = "pressures.mat"

# Import
f2cn = readdlm(join([dir_G, "f2cn.dat"]), '\t', Int, '\n')[:,1]
G = matread(join([dir_G, fn_G]))["G"]
# Import node info from .inp file
node_data = Any[]
open(join([dir_inp, fn_inp])) do f
    start_reading = false

    for line in eachline(f) 
        if occursin("*NODE, NSET=ALLNODES", line)   
            start_reading = true
            continue 
        end
        
        if start_reading && occursin("**", line)   
            break 
        end
        
        if start_reading
            node_line_str = split(strip(line), ',')
            node_line = Any[parse(Int, node_line_str[1])]
            append!(node_line, parse.(Float64, node_line_str[2:end]))
            push!(node_data, node_line)
        end
    end

end
nodes_inp = hcat(node_data...)' 
