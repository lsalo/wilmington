#
#
# Plot cell-based values for a 3D tetrahedral grid (fault cells are 2D cells)
#
# We define vertices and connectivity separately for each element, so that we can plot all vertices of each element with the corresponding cell value (plotting function poly below)
#
# To add packages, check status of project, etc:
# In Julia REPL: activate . (to activate environment in project folder)

# 0. Load modules 
# In Julia REPL: include("startup_plot_from_mrst.jl") # Run in julia REPL (once at the beginning of session).

# 1. Load data from MRST grid
# dir             = "/Users/lluis/Documents/MATLAB/mrst-dev/mrst-ls/ls-tests/simpleTetrahedralGrid/output/"
in_dir           = "/Users/lluis/Documents/MATLAB/mrst-dev/mrst-ls/ls-proj/wilmington/grid/"
figs_dir         = "/Users/lluis/Documents/julia/graphics/figs/wilmington/"
files_dir        = "/Users/lluis/Documents/julia/wilmington/grid/poly_cell_connectivity/"
sim_data_dir     = "/Users/lluis/Dropbox (MIT)/home_copy/matlab/sim_data/wilmington/test/"
write_connectivity = false
centroids       = readdlm(join([in_dir, "centroids.dat"]), '\t', Float64, '\n')
nodes_coords    = readdlm(join([in_dir, "coords.dat"]), '\t', Float64, '\n')
n2fn            = readdlm(join([in_dir, "n2fn.dat"]), '\t', Int, '\n')[:,1]
faces_nodes     = readdlm(join([in_dir, "nodes.dat"]), '\t', Int, '\n')[:,1]
cells_faces     = readdlm(join([in_dir, "faces.dat"]), '\t', Int, '\n')[:,1]
f2cn            = readdlm(join([in_dir, "f2cn.dat"]), '\t', Int, '\n')[:,1]
materialId      = readdlm(join([in_dir, "materialId.dat"]), '\t', Int, '\n')[:,1]
isFaultCell     = readdlm(join([in_dir, "isFaultCell.dat"]), '\t', Int, '\n')[:,1]
f2cn            = readdlm(join([in_dir, "f2cn.dat"]), '\t', Int, '\n')[:,1]
state0          = matread(join([sim_data_dir, "state0.mat"]))["data"]

#Adjust cell values
#cells_values = centroids[:,3]
#cells_values = state0["s"][:,2]
#cells_values = state0["pressure"]/1e6
#cells_values = state0["rsRatio"]
#cells_values = state0["poro"]
cells_values = state0["perm"][:,1]/9.869233e-16 #mD

# 2. Reshape faces_nodes
faces_nodes = get_faces_nodes(n2fn, faces_nodes)

# 3. Get connectivity for plot
if write_connectivity
    write_tetrahedral_cell_connectivity(f2cn, faces_nodes, cells_faces, nodes_coords, materialId, isFaultCell, files_dir)
end
cells_nodes_map, cells_materialId_map, cells_isFaultCell_map, coords_map, cells_values_map = load_cell_connectivity(f2cn, files_dir)

# 4. Figure
# get values at nodes and select a subset of cells based on materialId or isFaultCell?
material_ids = [1]
cell_subset = false
cells_nodes_map_new, cells_values_nodes = getValuesAtNodes(material_ids, cell_subset, f2cn, cells_nodes_map, cells_values, cells_values_map)

# Plot
#cmap = reverse(cgrad(:winter))
#cmap = cgrad(:winter)
#cmap = cgrad(:gist_stern)
fig_descriptor = (cmap=reverse(cgrad(:turbid)), cbar_label="Permeability [mD]", fig_size = (700, 800), plot_height = 1100, angles= (-0.1, -1.57, 0), save_fig = true, save_dir = figs_dir, title = "permeability")
plotCellValues(cells_nodes_map_new, coords_map, cells_values_nodes, fig_descriptor)





















