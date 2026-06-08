module Mrst
"""
This module performs grid-related operations for MRST (open-source FV simulator)
"""

using StatsBase
using Statistics
using RegularExpressions
using Distances
using DelimitedFiles

include("./arrayOps.jl")
using .ArrayOps


export find_fault_faces_mrst_grid, write_tetrahedral_grid_for_mrst, map_nodes_to_cells, get_faces_nodes


function find_fault_faces_mrst_grid(G::Dict, nset_mrst_ids::Vector, nsets_to_export::Vector{String}, dir_grid::String)
    """
    Finds face ids of MRST grid corresponding to nodesets in nset_mrst_ids. The nodesets must be numbered according to the MRST grid, as obtained from writeTetrahedralGridForMRST.
    """
    n2fn = Int.(G["faces"]["n2fn"])
    faces_nodes = get_faces_nodes(n2fn, Int.(G["faces"]["nodes"]))
    
    nset_mrst_face_ids = Vector{}(undef, length(nsets_to_export))
    for nset in eachindex(nset_mrst_ids)
        if !isempty(nset_mrst_ids[nset])          # nodeset, or a part of it, is present in mrst grid
            nset_mrst_face_ids[nset] = findall([all(row) for row in eachrow(in.(faces_nodes, [nset_mrst_ids[nset]]))])
            if !isempty(nset_mrst_face_ids[nset]) # at least one face defined by nodeset is present in mrst grid
                open(join([dir_grid, "faces_", nsets_to_export[nset], ".dat"]), "w") do f
                    for row in eachrow(nset_mrst_face_ids[nset])
                        println(f, join(row, " "))
                    end
                end
            else
                println(join(["This nodeset does not have any faces in MRST grid: ", nsets_to_export[nset]]))
            end
        else
            nset_mrst_face_ids[nset] = []
            println(join(["Nodeset not present MRST grid: ", nsets_to_export[nset]]))
        end
    end

    return faces_nodes, nset_mrst_face_ids

end


function write_tetrahedral_grid_for_mrst(nodes_inp, elements_inp, elset_names_inp, elsets_to_export, nset_info, dir_grid, fname_grid)
    """
    write nodal coordinates and tessellation to generate grid in MRST using:
    G = tetrahedralGrid(coord, t)

    Note: nodes_inp[:,1] (ids) can be numbered nonconsecutively, but they must be monotonically increasing.
    Note: may need to be adapted for .inp files without blocks or with a single block
    """
    # Select subset?
    if elsets_to_export == "all"
        T = vcat(elements_inp...)
        T = T[:, 2:end]
    else
        T = Matrix{Int64}[]
        for elset in eachindex(elsets_to_export)
            elset_id = findfirst(elset_names_inp .== elsets_to_export[elset])
            push!(T, elements_inp[elset_id][:, 2:end])
        end
        T = vcat(T...)
    end

    # For MRST, we need to renumber the node ids consecutively from 1 to n
    # This is slow for large grids
    T_mrst, map_inpNode_to_mrstNode = renumber_array_1_to_n(T)

    # Select correct nodes if only a subset is used
    unique_ids = sort(unique(T))
    id_nodes = in.(nodes_inp[:,1], [unique_ids])
    nodes_coord = nodes_inp[id_nodes, 2:end]

    # Write files
    open(join([dir_grid, fname_grid, "_coords.dat"]), "w") do f
        for row in eachrow(nodes_coord)
            println(f, join(row, " "))
        end
    end
    open(join([dir_grid, fname_grid, "_t.dat"]), "w") do f
        for row in eachrow(T_mrst)
            println(f, join(row, " "))
        end
    end
    open(join([dir_grid, "map_inp_to_mrst_node.dat"]), "w") do f
        for (k, v) in map_inpNode_to_mrstNode
            println(f, "$k\t$v")
        end
    end

    # Nodesets 
    if !isempty(nset_info)
        nset_inp = nset_info[1]
        nset_names_inp = nset_info[2]
        nsets_to_export = nset_info[3]
        
        # We need to renumber the nodes according to the map
        nset_mrst_ids = Vector{Vector{Int64}}(undef, length(nsets_to_export))
        for nset in eachindex(nsets_to_export)
            nset_id = findfirst(nset_names_inp .== nsets_to_export[nset])
            nset_node_ids = nset_inp[nset_id]
            nset_mrst_ids[nset] = [map_inpNode_to_mrstNode[k] for k in nset_node_ids if haskey(map_inpNode_to_mrstNode, k)]
            
            open(join([dir_grid, nsets_to_export[nset], ".dat"]), "w") do f
                for row in eachrow(nset_mrst_ids[nset])
                    println(f, join(row, " "))
                end
            end

        end

        return nset_mrst_ids

    end

end


function map_nodes_to_cells(G::Dict, n2fn::Union{Vector{Int64}, Matrix{Int64}}, f2cn::Union{Vector{Int64}, Matrix{Int64}})
    """
    Map nodes -> cells (i.e., all cells connected to each node) in a MRST grid
    Assumes that cells may have different number of faces (e.g., for grids with (n-1)-dimensional faults)
    """
    cells_faces = Int.(G["cells"]["faces"])
    faces_nodes = Int.(G["faces"]["nodes"])
    faces_nodes = get_faces_nodes(n2fn, faces_nodes)
    n_faces_per_cell = countmap(f2cn)
    unique_nc        = sort(unique(values(n_faces_per_cell)))

    # Create a map (Dictionary) for each group of cells with n faces
    maps = [Dict{Int,Vector{Int}}() for _ in 1:length(unique_nc)]
    for (i, nc) in enumerate(unique_nc)
        cells_with_n_faces = sort(collect(filter(x -> n_faces_per_cell[x] == nc, keys(n_faces_per_cell))))
        faces_of_cells_with_n_faces = cells_faces[findall(in(cells_with_n_faces), f2cn)]
        cells_nodes = faces_nodes[faces_of_cells_with_n_faces, :]
        cells_with_n_faces = repeat(cells_with_n_faces, inner=nc)
        for val in unique(cells_nodes) # Slow
            maps[i][val] = unique(cells_with_n_faces[findall(any(cells_nodes .== val, dims=2))])
        end        
    end

    # Collapse all maps into a single map that provides all cells for a given node
    nodes_to_cells = Dict{Int,Vector{Int}}()
    for map in maps
        for (k, v) in map
            nodes_to_cells[k] = union(get(nodes_to_cells, k, []), v)
        end
    end

    return nodes_to_cells, faces_nodes

end


function get_faces_nodes(n2fn::Union{Vector{Int64}, Matrix{Int64}}, faces_nodes::Union{Vector{Int64}, Matrix{Int64}})
    """
    Reshape MRST faces_nodes (G.faces.nodes) to (n_faces x n_nodes_per_face) 
    """
    #n_nodes = size(nodes_coords)[1]
    n_nodes_per_face = countmap(n2fn)
    n_faces          = maximum(keys(n_nodes_per_face))
    unique_nf        = unique(values(n_nodes_per_face))
    if length(unique_nf) == 1
        faces_nodes = reshape(faces_nodes, unique_nf[1], n_faces)'
    end
    
    return faces_nodes
end


end