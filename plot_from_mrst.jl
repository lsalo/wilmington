#
#
# Plot cell-based values for a 3D tetrahedral grid (fault cells are 2D cells)
#
# We define vertices and connectivity separately for each element, so that we can plot all vertices of each element with the corresponding cell value (plotting function poly below)
#

# 0. Load modules                           
# include("startup_plot_from_mrst.jl") # Run in julia REPL (once at the beginning of session).

# 1. Load data from MRST grid
# dir             = "/Users/lluis/Documents/MATLAB/mrst-dev/mrst-ls/ls-tests/simpleTetrahedralGrid/output/"
in_dir           = "/Users/lluis/Documents/MATLAB/mrst-dev/mrst-ls/ls-proj/wilmington/grid/"
figs_dir         = "/Users/lluis/Documents/julia/graphics/figs/wilmington/"
files_dir        = "/Users/lluis/Documents/julia/wilmington/grid/poly_cell_connectivity/"
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

#Adjust cell values
#minx, miny = minimum(centroids[:,1]), minimum(centroids[:,2])
#cells_values = mean([centroids[:,1].-minx centroids[:,2].-miny], dims=2)
# cells_values = 1:112573
cells_values = centroids[:,3]

# 2. Reshape faces_nodes
faces_nodes = get_faces_nodes(n2fn, faces_nodes)

# 3. Get connectivity for plot
if write_connectivity
    write_tetrahedral_cell_connectivity(f2cn, faces_nodes, cells_faces, nodes_coords, materialId, isFaultCell, files_dir)
end
cells_nodes_map, cells_materialId_map, cells_isFaultCell_map, coords_map, cells_values_map = load_cell_connectivity(f2cn, files_dir)

# 4. Figure
# Select a subset of cells based on materialId or isFaultCell?
material_ids = [1]
cell_subset = false
n_faces_per_cell = countmap(f2cn)
unique_nc        = sort(unique(values(n_faces_per_cell)))
if cell_subset
    cells_nodes_map_new  = Matrix{Int64}[]
else
    cells_nodes_map_new = cells_nodes_map
end
cells_values_nodes  = Vector{Float64}[]

for f in 1:length(unique_nc)
    if cell_subset
        cells_nodes_id = findall(in(material_ids), cells_materialId_map[f])
        push!(cells_nodes_map_new, cells_nodes_map[f][cells_nodes_id, :])
    end
    push!(cells_values_nodes, cells_values[cells_values_map[f]])
end


# Plot
# NOTE: In the plot loop, the order with which each cell face is plotted will lead to some differences for interior (shared) faces. Here we are not plotting based on global cell index since we have data organized based on number of faces of each cell, so it will not be exactly equal to MRST view.

minv, maxv = minimum(flatten(cells_values_nodes)), maximum(flatten(cells_values_nodes))
#minv, maxv = 1, 67733
#cmap = reverse(cgrad(:turbid))
cmap = cgrad(:gist_stern)
size_in_inches = (8, 6)
dpi = 300
size_in_pixels = size_in_inches .* dpi
fig = Figure(resolution=size_in_pixels)
#ax1 = Axis3(fig[1, 1], aspect=(1,1,1/5), perspectiveness=0.5)
#ax1.zreversed = true
lscene = LScene(fig[1, 1], show_axis=false)
plotcbar = true
for n in eachindex(cells_nodes_map_new)
    if length(cells_nodes_map_new[n])>0
        coords_plot = coords_map[n]
        if coords_plot[1,3] > 0
            coords_plot[:,end] = coords_plot[:,end]*-1
        end
        pl = poly!(lscene, coords_plot, cells_nodes_map_new[n], color=cells_values_nodes[n], colormap=cmap, colorrange=(minv, maxv), strokecolor=:black, strokewidth=1, shading=false)
        if plotcbar == true
            Colorbar(fig[1, 2], pl, label="values", height=Relative(0.5))
            global plotcbar = false
        end
    end
end
colgap!(fig.layout, 5)
fig

savefig = false
if savefig
    # For interactive plots with GLMakie, you can rotate, zoom, etc and it will save the view on the fig window
    save(join([output_dir, "test.png"]), fig)
end




















