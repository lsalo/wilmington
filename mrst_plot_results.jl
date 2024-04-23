#
#
# Plot cell-based values for a 3D tetrahedral grid (fault cells are 2D cells)
#
# We define vertices and connectivity separately for each element, so that we can plot all vertices of each element with the corresponding cell value (plotting function poly below)
#
# To add packages, check status of project, etc:
# In Julia REPL: ] activate . or using Pkg; Pkg.activate("path/to/project') (to activate environment in project folder)

# 0. Load modules
# REMEMBER: Update path in include file below (windows vs mac)
# In Julia REPL: include("startup_mrst_plot_results.jl") # Run in julia REPL (once at the beginning of session).

# 1. Load data from MRST grid
in_dir           = "C:\\Users\\lsalo\\matlab\\mrst-dev\\mrst-ls\\ls-proj\\wilmington\\grid\\"             #"/Users/lluis/Documents/MATLAB/mrst-dev/mrst-ls/ls-proj/wilmington/grid/"
#figs_dir         = "C:\\Users\\lsalo\\julia\\graphics\\figs\\wilmington\\JSsetup\\initialization\\"       #"/Users/lluis/Documents/julia/graphics/figs/wilmington/"
figs_dir         = "C:\\Users\\lsalo\\julia\\graphics\\figs\\wilmington\\revised_20y\\initialization\\"  
#figs_dir         = "C:\\Users\\lsalo\\julia\\graphics\\figs\\wilmington\\JSsetup\\results\\"       #"/Users/lluis/Documents/julia/graphics/figs/wilmington/"
files_dir        = "C:\\Users\\lsalo\\julia\\wilmington\\grid\\poly_cell_connectivity\\"                  #"/Users/lluis/Documents/julia/wilmington/grid/poly_cell_connectivity/"
sim_data_ini_dir     = "C:\\Users\\lsalo\\matlab\\sim_data\\mrst\\wilmington\\liquidRate_JS_LSprops_owc1800_test\\ini_data\\"
# "C:\\Users\\lsalo\\matlab\\sim_data\\mrst\\wilmington\\runJSsetup\\ini_data\\"         #"/Users/lluis/Dropbox (MIT)/home_copy/matlab/sim_data/wilmington/test/"
sim_data_fin_dir     = "C:\\Users\\lsalo\\matlab\\sim_data\\mrst\\wilmington\\runJSsetup\\result_data\\"  

write_connectivity = false
centroids       = readdlm(join([in_dir, "centroids.dat"]), '\t', Float64, '\n')
nodes_coords    = readdlm(join([in_dir, "coords.dat"]), '\t', Float64, '\n')
n2fn            = readdlm(join([in_dir, "n2fn.dat"]), '\t', Int, '\n')[:,1]
faces_nodes     = readdlm(join([in_dir, "nodes.dat"]), '\t', Int, '\n')[:,1]
cells_faces     = readdlm(join([in_dir, "faces.dat"]), '\t', Int, '\n')[:,1]
f2cn            = readdlm(join([in_dir, "f2cn.dat"]), '\t', Int, '\n')[:,1]
materialId      = readdlm(join([in_dir, "materialId.dat"]), '\t', Int, '\n')[:,1]
isFaultCell     = readdlm(join([in_dir, "isFaultCell.dat"]), '\t', Int, '\n')[:,1]

# state 0
state0          = matread(join([sim_data_ini_dir, "state0.mat"]))["data"]

# results (for larger files, we save in .mat binary format, it is much faster to read than text files)
#p               = matread(join([sim_data_fin_dir, "states_p_Pa.mat"]))["p"]
#dp_fin_MPa      = (p[:,end]-p[:,1])/1e6 
#sw              = matread(join([sim_data_fin_dir, "states_sw.mat"]))["sw"]
#so              = matread(join([sim_data_fin_dir, "states_so.mat"]))["so"]
#sg              = matread(join([sim_data_fin_dir, "states_sg.mat"]))["sg"]



#Adjust cell values
#cells_values = centroids[:,3]/1e3
#cells_values = state0["s"][:,3]
cells_values = state0["pressure"]/1e6
#cells_values = state0["rsRatio"]
#cells_values = state0["rs"]
#cells_values = state0["poro"]
#cells_values = state0["perm"][:,1]/9.869233e-16 #mD

#sw = state0["s"][:,1]
#cells_values = ones(Int, size(sw, 1))
#cells_values[state0["s"][:,2] .> 1e-4] .= 2
#cells_values[state0["s"][:,3] .> 1e-4] .= 3

#cells_values = dp_fin_MPa
#cells_values = ones(Int, size(sw, 1))
#cells_values[so[:,1] .> 1e-4] .= 2
#cells_values[sg[:,1] .> 1e-4] .= 3 

# 2. Reshape faces_nodes
faces_nodes = get_faces_nodes(n2fn, faces_nodes)

# 3. Get connectivity for plot
if write_connectivity
    write_tetrahedral_cell_connectivity(f2cn, faces_nodes, cells_faces, nodes_coords, materialId, isFaultCell, files_dir)
end
cells_nodes_map, cells_materialId_map, cells_isFaultCell_map, coords_map, cells_values_map = load_cell_connectivity(f2cn, files_dir)

# 4. Figure
# get values at nodes and select a subset of cells based on materialId or isFaultCell
# within_reservoir_faults_id = [102, 104, 105, 106, 109, 111, 115, 119, 126, 127, 208, 210, 213, 221, 310, 321]
# non_reservoir_faults_id = [50, 51, 53, 252, 308, 313, 352, 452]
material_ids = [1]
cell_subset = true
cells_nodes_map_new, cells_values_nodes = getValuesAtNodes(material_ids, cells_materialId_map, cell_subset, f2cn, cells_nodes_map, cells_values, cells_values_map)

# Plot
#mycmap = reverse(cgrad(:winter))
#mycmap = cgrad(:winter)
#mycmap = cgrad(:gist_stern)
#mycmap = reverse(cgrad(:turbid))
#mycmap = ColorScheme([colorant"blue", colorant"green", colorant"red"])
mycmap=cgrad(:nipy_spectral)

# Full mesh
#fig_descriptor = (cmap=mycmap, cbar_label="So [-]", fig_size = (1000, 800), plot_height = 800, cbar_height = 200, plot_azi = 1.8*pi, plot_elev = 0.3*pi, plot_linew = 0.1, xlms=(nothing, nothing), ylms=(nothing, nothing), save_fig = false, save_dir = figs_dir, title = "s0", target_dpi=450, target_size_in=(8,6.4))

# Reservoir 
fig_descriptor = (cmap=mycmap, cbar_label="Pressure [MPa]", fig_size = (800, 600), plot_height = 800, cbar_height = 200, plot_azi = 1.8*pi, plot_elev = 0.25*pi, plot_linew=0.1, xlms=(378546, 400263), ylms=(3.728e6, 3.751e6), save_fig = false, save_dir = figs_dir, title = "s0_reservoir", target_dpi=450, target_size_in=(8,6))

# Wilmington fault
#fig_descriptor = (cmap=mycmap, cbar_label="Sg [-]", fig_size = (1600, 1200), plot_height = 1600, cbar_height = 200, plot_azi = 1.7*pi, plot_elev = 0.1*pi, plot_linew=0.05, xlms=(378546, 400263), ylms=(3.728e6, 3.751e6), save_fig = false, save_dir = figs_dir, title = "So_wilmingtonFault", target_dpi=450, target_size_in=(8,6))

# Plot
plotCellValues(cells_nodes_map_new, coords_map, cells_values_nodes, fig_descriptor)





















