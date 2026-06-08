module MakiePolyPlot
"""
This module contains functionality for making "poly" plots with Makie based on simulation results at cell centers or nodes.
"""

using StatsBase
using DelimitedFiles
using StatsBase, Base.Iterators
using Colors
using GLMakie
GLMakie.activate!()

export write_tetrahedral_cell_connectivity, load_cell_connectivity, getValuesAtNodes, plotNodalValues

# Generate ploy-specific connectivity
function write_tetrahedral_cell_connectivity(f2cn, faces_nodes, cells_faces, nodes_coords, materialId, isFaultCell, dir=[])
    """
    Write connectivity to plot results from MRST (values obtained at cell centers, but poly uses values at the nodes).

    This function is intended for usage with the poly function in Makie.jl. The main steps are:
    1. Get nodes for each cell face (cells_nodes)
    2. Regenerate cell-specific connectivity for poly (coloring based on nodal values, so we have to duplicate nodes to have distinct nodes for each cell)

    Not intended to be general. It is specific for Cubit tetrahedral grids where the fault cells are 2D entities. The inputs to this function, which describe the grid connectivity, are in MRST format. Note that because fault cells are 2D, fault faces are 1D edges and therefore the third column in faces_nodes is dummy.

    This function is currently slow due to the inner loop over all cells, but this accommodates variability in cell_nodes. We typically run it once, save the connectivity, and then simply load it.
    """
    n_faces_per_cell = countmap(f2cn)
    unique_nc        = sort(unique(values(n_faces_per_cell)))
    cells_nodes_map  = Matrix{Int64}[]
    coords_map       = Matrix{Float64}[]
    cells_values_map = Vector{Int64}[]
    cells_materialId_map = Vector{Int64}[]
    cells_isFaultCell_map = Vector{Int64}[]

    for n in 1:length(unique_nc)
        cells_with_n_faces = sort(collect(filter(x -> n_faces_per_cell[x] == n, keys(n_faces_per_cell))))
        n_cells = length(cells_with_n_faces)
        faces_of_cells_with_n_faces = cells_faces[findall(in(cells_with_n_faces), f2cn)]
        cells_nodes = faces_nodes[faces_of_cells_with_n_faces, :]
        materialId_of_cells_with_n_faces = materialId[cells_with_n_faces]
        isFaultCell_of_cells_with_n_faces = isFaultCell[cells_with_n_faces]
        f2cn_cells_new = repeat(1:n_cells, inner=n)
        cells_nodes_itmap = zeros(Int64, size(cells_nodes))
        coords_itmap = zeros(Float64, 0, 3)
        cells_values_itmap = Vector{Int64}[]
        cells_materialId_itmap = Vector{Int64}[]
        cells_isFaultCell_itmap = Vector{Int64}[]
        to_add = 0

        for cell in 1:n_cells
            cell_nodes = cells_nodes[f2cn_cells_new .== cell, :]
            unique_cell_nodes = sort(unique(cell_nodes))
            cell_nodes_map = Dict(unique_cell_nodes[i] => i+to_add for i in 1:lastindex(unique_cell_nodes))
            cell_nodes_new = [cell_nodes_map[node] for node in cell_nodes]
            if length(unique_cell_nodes)==3 && n>1      # this will not be general
                cells_nodes_itmap[f2cn_cells_new .== cell, :] .= unique(cell_nodes_new')' 
            else
                cells_nodes_itmap[f2cn_cells_new .== cell, :] .= cell_nodes_new
            end
            cells_materialId_itmap = [cells_materialId_itmap; [materialId_of_cells_with_n_faces[cell] for i=1:n]]
            cells_isFaultCell_itmap = [cells_isFaultCell_itmap; [isFaultCell_of_cells_with_n_faces[cell] for i=1:n]]
            coords_itmap = vcat(coords_itmap, nodes_coords[unique_cell_nodes, :])
            cells_values_itmap = vcat(cells_values_itmap, repeat([cells_with_n_faces[cell]], inner=length(unique_cell_nodes)))
            
            # Prepare next timestep
            to_add = maximum(cells_nodes_itmap)
        end

        # Delete repeated faces
        unique_id = unique(i -> cells_nodes_itmap[i,:], 1:size(cells_nodes_itmap, 1))
        push!(cells_nodes_map, cells_nodes_itmap[unique_id, :])
        push!(cells_materialId_map, cells_materialId_itmap[unique_id])
        push!(cells_isFaultCell_map, cells_isFaultCell_itmap[unique_id])
        push!(coords_map, coords_itmap)
        push!(cells_values_map, cells_values_itmap) 

        if !isempty(dir)
            variables = Dict{String, Any}()
            variables["cells_nodes_map_$n"] = cells_nodes_map[n]
            variables["cells_materialId_map_$n"] = cells_materialId_map[n]
            variables["cells_isFaultCell_map_$n"] = cells_isFaultCell_map[n]
            variables["coords_map_$n"] = coords_map[n]
            variables["cells_values_map_$n"] = cells_values_map[n]
    
            for vname in keys(variables)
                open(join([dir, "$vname.dat"]), "w") do f
                    for row in eachrow(variables[vname])
                        println(f, join(row, ' '))
                    end
                end
            end   
        
        end
    
    end

    if isempty(dir)
        return cells_nodes_map, cells_materialId_map, cells_isFaultCell_map, coords_map, cells_values_map
    end

end

function load_cell_connectivity(f2cn, dir)
    """
    TBD
    """
    n_faces_per_cell = countmap(f2cn)
    unique_nc        = sort(unique(values(n_faces_per_cell)))
    cells_nodes_map  = Matrix{Int64}[]
    cells_materialId_map  = Vector{Int64}[]
    cells_isFaultCell_map  = Vector{Int64}[]
    coords_map       = Matrix{Float64}[]
    cells_values_map = Vector{Int64}[]
    
    for n in 1:length(unique_nc)
        push!(cells_nodes_map, readdlm(join([dir, "cells_nodes_map_$n.dat"]), Int))
        push!(cells_materialId_map, readdlm(join([dir, "cells_materialId_map_$n.dat"]), Int)[:,1])
        push!(cells_isFaultCell_map, readdlm(join([dir, "cells_isFaultCell_map_$n.dat"]), Int)[:,1])
        push!(coords_map, readdlm(join([dir, "coords_map_$n.dat"]), Float64))
        push!(cells_values_map, readdlm(join([dir, "cells_values_map_$n.dat"]), Int64)[:,1])
    end

    return cells_nodes_map, cells_materialId_map, cells_isFaultCell_map, coords_map, cells_values_map
end


# Get nodal values, if quantity of interest is obtained at the cell-level
function getValuesAtNodes(unit_ids, cells_materialId_map, cell_subset, f2cn, cells_nodes_map, cells_values, cells_values_map)
    """
    MRST-specific.

    assign cell values to cell nodes (i.e., all faces of each cell, defined by the node connectivities). If desired, select a subset of cell values to plot.
    """
    n_faces_per_cell = countmap(f2cn)
    unique_nc        = sort(unique(values(n_faces_per_cell)))
    if cell_subset
        cells_nodes_map_plot  = Matrix{Int64}[]
    else
        cells_nodes_map_plot = cells_nodes_map
    end
    cells_values_nodes  = Vector{Float64}[]

    for f in 1:length(unique_nc)
        if cell_subset
            cells_nodes_id = findall(in(unit_ids), cells_materialId_map[f])
            push!(cells_nodes_map_plot, cells_nodes_map[f][cells_nodes_id, :])
        end
        push!(cells_values_nodes, cells_values[cells_values_map[f]])
    end

    return cells_nodes_map_plot, cells_values_nodes
end


# Poly Plot
function plotNodalValues(cells_nodes_map_plot, coords_map, cells_values_nodes, fig_descriptor;  
                         node_coords_plot=nothing, line_data=nothing)
    """
    plot values at the nodes (e.g., permeability, pressure, displacement, stress...) at a given timestep, and save .png to desired folder.

    NOTE: In the plot loop, the order with which each cell face is plotted will lead to some differences for interior (shared) faces. Here we are not plotting based on global cell index since we have data organized based on number of faces of each cell, so it will not be exactly equal to MRST view.
    """
    val = Vector{Float64}[]
    for f in 1:lastindex(cells_nodes_map_plot)   # get values for colormap limits in case of cell subset
        idx_val = unique(collect(flatten(cells_nodes_map_plot[f])))
        push!(val, collect(flatten(cells_values_nodes[f]))[idx_val])
    end
    val = collect(flatten(val))

    #minv, maxv = minimum(val), maximum(val)
    if all(val .>= 0) || all(val .<= 0)
        minv, maxv = minimum(val), maximum(val)
    else
        v_max_abs = maximum(abs.(val))
        minv, maxv = -v_max_abs, v_max_abs
    end
    #minv, maxv = 0, maximum(val)
    #minv, maxv = -60, 60
    
    fig = Figure(size=fig_descriptor.fig_size)
    ax1 = Axis3(fig[1, 1], aspect=(1.0, 1.0, 0.5), height=fig_descriptor.plot_height, azimuth=fig_descriptor.plot_azi, elevation=fig_descriptor.plot_elev, viewmode=:fitzoom) #azimuth default = 1.275pi, elevation default = pi/8
    #ax1.zreversed = true
    #lscene = LScene(fig[1, 1], show_axis=false, height=fig_descriptor.plot_height, valign=:center, halign=0)
    plotcbar = true
    for n in eachindex(cells_nodes_map_plot)
        if length(cells_nodes_map_plot[n])>0
            coords_plot = coords_map[n]
            if coords_plot[1,3] > 0
                coords_plot[:,end] = coords_plot[:,end]*-1
            end
            pl = poly!(ax1, coords_plot, cells_nodes_map_plot[n], color=cells_values_nodes[n], colormap=fig_descriptor.cmap, colorrange=(minv, maxv), strokecolor=Makie.RGBA(0, 0, 0, 0.2), strokewidth=fig_descriptor.plot_linew, alpha=1.0)
            if plotcbar == true
                cbar = Colorbar(fig[1, 2], pl, label=fig_descriptor.cbar_label, labelsize=18, height=fig_descriptor.cbar_height)
                #cbar.ticks = range(round(minv,digits=0), stop=round(maxv,digits=0), length=6)
                #Colorbar(fig[1, 2], pl, label=fig_descriptor.cbar_label, height=fig_descriptor.cbar_height, ticks=([1.0,2.0,3.0], ["Water", "Oil", "Gas"]))
                plotcbar = false
            end
        end
    end

    # eg to highlight basement-sedimentary interface on faults represented as nodesets in mesh
    if node_coords_plot !== nothing 
        mycolor = RGBA(0, 0, 0, 1)
        adjust_coords = true
        for line in eachindex(node_coords_plot)
            if adjust_coords == true
                vx = node_coords_plot[line][:, 1].*1.000001
                vy = node_coords_plot[line][:, 2].*1.000001
            else
                vx = node_coords_plot[line][:, 1]
                vy = node_coords_plot[line][:, 2]
            end
            scatter!(ax1, vx, vy, node_coords_plot[line][:, 3], color=mycolor, markersize=4, marker=:circle)
        end
    end

    # Plot lines (eg coastline, boundaries, other fault traces, etc)
    if line_data !== nothing
        for line in line_data
            if size(line, 2) == 3
                lines!(ax1, line[:, 1], line[:, 2], line[:, 3], color=:black, linewidth=1, linestyle=:solid)
            else
                @warn "Skipping line - expected Nx3 matrix with x, y, z coordinates"
            end
        end
    end
    colgap!(fig.layout, 0)
    xlims!(ax1, fig_descriptor.xlms)
    ylims!(ax1, fig_descriptor.ylms) 
    hidedecorations!(ax1)  # hides ticks, grid and labels
    hidespines!(ax1) 
    display(fig)   

    if fig_descriptor.save_fig
        ppu = fig_descriptor.target_size_in[1]*fig_descriptor.target_dpi/fig_descriptor.fig_size[1]
        # For interactive plots with GLMakie: you can also rotate, zoom, etc on the GLMakie window, and then save, and it will save the view on the fig window.
        save(join([fig_descriptor.save_dir, fig_descriptor.title, ".png"]), fig, px_per_unit=ppu)
    end

end


end