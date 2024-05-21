# Code do compute pore pressure changes at nodal points (for Abaqus) based on cell-based pore pressure values from reservoir simulation (MRST)

# startup 
cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__)

using MAT
using DelimitedFiles
using Revise
using StatsBase
using Printf

using GLMakie
GLMakie.activate!()
using Colors
using ColorSchemes
using LaTeXStrings

lib_dir = "C:/Users/lsalo/julia/jl-tools/"
includet(joinpath(lib_dir, "modules", "abaqus.jl"))
using .Abaqus
includet(joinpath(lib_dir, "modules", "mrst.jl"))
using .Mrst
includet(joinpath(lib_dir, "modules", "makiePolyPlot.jl"))
using .MakiePolyPlot

# Directories
dir_inp = "C:/Users/lsalo/Wilmington/AbaqusGrids/"
fn_inp  = "WilmingtonMesh_AllFaults_Rotated_Coarse.inp"
dir_G   = "C:/Users/lsalo/matlab/mrst-dev/mrst-ls/ls-proj/wilmington/grid/"
fn_G    = "wilmington_ranger_G_withGeometry.mat"
dir_p   = "C:/Users/lsalo/matlab/sim_data/mrst/wilmington/liquidRate_JS_LS_lowRsprops_owc1800/result_data/1936_1975/"
dir_p_abq = "C:/Users/lsalo/matlab/sim_data/mrst/wilmington/liquidRate_JS_LS_lowRsprops_owc1800/result_data/dp_for_abaqus_1936_1975/"
fn_p    = "states_p_Pa.mat"
fn_t    = "states_t_s.mat"
figs_dir = "C:/Users/lsalo/julia/graphics/figs/wilmington/liquidRate_JS_LS_lowRsprops_owc1800/"  

# Connectivity from .inp file (Abaqus grid)
nodes_inp, _, _, _, _ = get_connectivity_from_inp(dir_inp, fn_inp, "tetrahedral")

# Node-cell connections map
G           = matread(join([dir_G, fn_G]))["G"]
f2cn        = Int.(G["cells"]["f2n"])
n2fn        = Int.(G["faces"]["n2fn"])
nodes_to_cells, faces_nodes = map_nodes_to_cells(G, n2fn, f2cn)

# Get nodal pressures and compute Dp
# This needs to be done at each timestep, and then compute Dp
t = 1936 .+ matread(join([dir_p fn_t]))["t"]./31556952  # convert to year
p = matread(join([dir_p fn_p]))["p"]
p0 = matread(join([dir_p "state0_p_Pa.mat"]))["p0"]
#dp = p[:, 1:end] - p0  # dp from previous tstep (not from hydrostatic pressure)
dp = p[:, 1:end] .- p[:, 1]     # bad convergence with p0, dp from first tstep (no inj yet)
nodal_dp = Vector{}(undef, size(dp, 2))
for i in 1:length(nodal_dp)
    nodal_dp[i] = Dict((node_id, mean(dp[cell_ids, i])) for (node_id, cell_ids) in nodes_to_cells)
end

# Update MRST nodes to global (.inp) nodes using map
# Note: Given that we are using the "not rotated" grid in MRST, and the rotated grid for Abaqus, this works because the global node ids (in .inp file) of both the rotated and not rotated grids remain the same. Otherwise, we would have to use the rotated grid in MRST as well. We decided not to use the rotated grid in MRST because the well coordinates (to assign well cell ids) are reported in the "not rotated" frame, although you could move from one to the other by matching cell volumes + depth. See create_mrst_inputs.jl.
map_inp_to_mrst_node = readdlm(join([dir_G, "map_inp_to_mrst_node.dat"]), '\t', Int64, '\n')
nodal_dp_inp = Matrix{Union{Int64, Float64}}(undef, Int(G["nodes"]["num"]), size(dp, 2)+1)
nodal_dp_inp[:,1] = map_inp_to_mrst_node[:,1]
for i in 1:length(nodal_dp)
    nodal_dp_inp[:, i+1] = [nodal_dp[i][row[2]] for row in eachrow(map_inp_to_mrst_node)]
end
nodal_dp_inp_sorted = sortslices(nodal_dp_inp, dims=1)  # sort based on node id

# For each tstep, export nodal values which are != 0
fn = "nodal_dp_"
for i in 1:length(nodal_dp)
    data = hcat(nodal_dp_inp_sorted[:,1], nodal_dp_inp_sorted[:,i+1])
    open(join([dir_p_abq, fn, string(round(t[i], digits=2)), ".dat"]), "w") do f
        for row in eachrow(data)
            if abs(row[2]) > 10    # we consider only values > 10 Pa 
                #writedlm(f, row', delim=",")
                v = map(x -> @sprintf("%.2f", x), row[2]) # NOTE: convert fromat because notation like 1.315e7 may not work in Abaqus
                str=string(row[1])*","*v*"\n"
                write(f, str)
            end
        end
    end
end

# Plot nodal dp using poly
# get dp with respect to initial 
nodal_dp_inp_plot = Matrix{Union{Int64, Float64}}(undef, size(nodal_dp_inp)[1],  size(nodal_dp_inp)[2])
nodal_dp_inp_plot[:,1] = map_inp_to_mrst_node[:,2]
nodal_dp_inp_plot[:,2:end] = nodal_dp_inp[:,2:end]
nodal_dp_inp_plot = sortslices(nodal_dp_inp_plot, dims=1)  

# put variables in vector containers
t_id = 309
coords_map = Vector{}(undef, 1)
cells_nodes_map_plot = Vector{}(undef, 1)
cells_values_nodes = Vector{}(undef, 1)
coords_map[1] = G["nodes"]["coords"]
cells_nodes_map_plot[1] = faces_nodes
cells_values_nodes[1] = nodal_dp_inp_plot[:,t_id]./1e5

# plot
mycmap=cgrad(:curl)

fig_descriptor = (cmap=mycmap, cbar_label=L"$\Delta p$ [bar]", fig_size = (800, 600), plot_height = 900, cbar_height = 200, plot_azi = 1.57*pi, plot_elev = 0.47*pi, plot_linew=0.1, xlms=(378449, 400339), ylms=(3.728e6, 3.752e6), save_fig = true, save_dir = figs_dir, title = "dp_reservoir_1961", target_dpi=350, target_size_in=(8,6))

plotNodalValues(cells_nodes_map_plot, coords_map, cells_values_nodes, fig_descriptor)