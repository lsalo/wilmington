module Abaqus
    """
    Process abaqus .inp files
    """

    using StatsBase
    using Statistics
    using RegularExpressions
    using Distances
    using DelimitedFiles

    export get_element_centroids_inp, get_face_connectivity_elsets, get_face_connectivity_nodesets_C3D4, get_face_connectivity_nodesets, get_connectivity_from_inp


    function get_element_centroids_inp(elements_all::Matrix, nodes_inp::AbstractMatrix)
        """
        Get element centroids of Abaqus .inp file

        Note: because the in. function is very slow, in the loop we changed:
            row_ids = in.(nodes_inp[:, 1], Ref(nodes))  
            nodes_coords = nodes_inp[row_ids, 2:4]
            
            by:
            nodes_coords = [nodes_inp_dict[node] for node in nodes]
        """
        n_elements = size(elements_all)[1]
        centroids = Matrix{Union{Int64, Float64}}(undef, n_elements, 4)
        centroids[:, 1] = 1:n_elements
        nodes_inp_dict = Dict(zip(nodes_inp[:, 1], eachrow(nodes_inp[:, 2:4])))     # for speed
        
        for element in 1:n_elements
            nodes = Set(@view elements_all[element, 2:end])                         # Set and @view for speed
            nodes_coords = [nodes_inp_dict[node] for node in nodes]                 
            centroids[element, 2:4] = mean(hcat(nodes_coords...)', dims=1)
        end

        return centroids
    end


    function get_face_connectivity_elsets(elsets_to_export::Vector{String}, elset_names_inp::Vector{String}, elements_inp::AbstractArray, element_type::String)
        """
        Uses outputs from get_connectivity_from_inp and provides face connectivity for specific elementsets.
        """

        if element_type == "tetrahedral"
            nnf = 3                 # Number of nodes per face
            nfe = 4                 # Number of faces per element
        end


        face_connectivity_elsets = Vector{}(undef, length(elsets_to_export))
        for elset in eachindex(elsets_to_export)
            elset_id = findfirst(elset_names_inp .== elsets_to_export[elset])
            node_list =  elements_inp[elset_id][:, 2:end]
            
            face_conn_elem = Matrix{Int64}(undef, 4, nnf) 
            face_connectivity_elset = Matrix{Int64}(undef, size(node_list)[1]*4, nnf) 
            for (i, row) in enumerate(eachrow(node_list))
                row = sort(row)
                if nfe == 4
                    face_conn_elem = [row[1:3]'; row[1:2]' row[4]; row[1] row[3:4]'; row[2:4]']  
                end   
                face_connectivity_elset[(nfe*i-(nfe-1)):(nfe*i), :] = face_conn_elem                  
            end 
            
            face_connectivity_elsets[elset] = unique(face_connectivity_elset, dims=1)    
        end

        return face_connectivity_elsets

    end


    function get_face_connectivity_nodesets_C3D4(nsets_to_export::Vector{String}, nset_names_inp::Vector{String}, nset_inp::Vector, elements_inp::AbstractMatrix)
        """
        Uses outputs from get_connectivity_from_inp and provides face connectivity for specific nodesets.
        The face connectivity is defined as a matrix where each row contains the nodes defining a face of the nodeset.
        Separate function because this info is not found directly in the .inp file, and may not be needed.
        
        Note:
        elements_inp is a matrix with the element connectivity. All faces (in all nodesets) must be contained in the elements provided.
        """
        nnf = 3 # C3D4 faces are triangles

        face_connectivity_nodesets = [zeros(Int64, 0, nnf) for _ in 1:length(nsets_to_export)]
        for nset in eachindex(nsets_to_export)
            nset_id = findfirst(nset_names_inp .== nsets_to_export[nset])
            nset_list = nset_inp[nset_id]

            for row in eachrow(elements_inp[:, 2:end])  # row has 4 nodes for C3D4
                in_nset = in.(row, Ref(nset_list))
                k = sum(in_nset)

                if k == nnf
                    face_connectivity_nodesets[nset] = vcat(face_connectivity_nodesets[nset], sort(row[in_nset])')
                elseif k == 4
                    a, b, c, d = row
                    face_connectivity_nodesets[nset] = vcat(face_connectivity_nodesets[nset], sort!([a, b, c])', sort!([a, b, d])', sort!([a, c, d])', sort!([b, c, d])')
                end
            end

            face_connectivity_nodesets[nset] = unique(face_connectivity_nodesets[nset], dims=1)
            print("Face connectivity for nodeset $(nsets_to_export[nset]) completed.\n")
        end

        return face_connectivity_nodesets
    end


    function get_face_connectivity_nodesets(nsets_to_export::Vector{String}, nset_names_inp::Vector{String}, nset_inp::Vector, elements_inp::AbstractMatrix, element_type::String)
        """
        Uses outputs from get_connectivity_from_inp and provides face connectivity for specific nodesets.
        The face connectivity is defined as a matrix where each row contains the nodes defining a face of the nodeset.
        Separate function because this info is not found directly in the .inp file, and may not be needed.
        
        Note:
        elements_inp is a matrix with the element connectivity. All faces (in all nodesets) must be contained in the elements provided.
        """
        if element_type == "tetrahedral"
            nnf = 3                 # Number of nodes per face
        end
        
        face_connectivity_nodesets = [zeros(Int64, 0, nnf) for _ in 1:length(nsets_to_export)]
        for nset in eachindex(nsets_to_export)
            nset_id = findfirst(nset_names_inp .== nsets_to_export[nset])
            nset_list =  nset_inp[nset_id]
        
            for row in eachrow(elements_inp[:, 2:end])
                in_nset = in.(row, Ref(nset_list))
                if sum(in_nset) == nnf                                              # element contains a face in the nodeset
                    face_connectivity_nodesets[nset] = vcat(face_connectivity_nodesets[nset], sort(row[in_nset])')
                elseif sum(in_nset) > nnf
                    #error("Multiple nodeset/surface faces in a single element? Check")
                    k = count(in_nset)
                    @warn "Element has $k nodes in nodeset (nnf=$nnf). Nodeset probably includes entire element." element_nodes=row nset=nsets_to_export[nset]
                end
            end
            face_connectivity_nodesets[nset] = unique(face_connectivity_nodesets[nset], dims=1)     # remove duplicate faces
        
            print(join(["Face connectivity for nodeset ", nsets_to_export[nset], " completed.\n"]))
        end
        
        return face_connectivity_nodesets
        
    end


    function get_connectivity_from_inp(dir_inp::String, fn_inp::String, element_type::String, parse_elsets=true)
        """
        Reads an Abaqus .inp file and imports:
        - nodes_inp: n_node x 4 matrix, with node ids and x, y, z coordinates
        - elements_inp: n_elset x (n_elements x (n_point+1), with element ids and n-point tessellation (e.g., 4 nodes for each tetrahedral element)
        - nsets_inp: n_nodesets x (n_nodes_in_nodeset x 1)

        Currently limited to element_type = "tetrahedral"
        """

        node_data = Any[]
        elements_inp = Any[]
        elset_names_inp = String[]
        nset_names_inp = String[]
        nset_inp = Any[]
        open(join([dir_inp, fn_inp])) do f
            # Nodes
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

            # Elements
            seekstart(f)
            start_reading = false
            elset_count = 0
            elements = Any[]
            if element_type == "tetrahedral"
                for line in eachline(f)
                    if occursin("*ELEMENT, TYPE=C3D4, ELSET=", line)   
                        start_reading = true
                        m = match(r"ELSET=(.*)", line)
                        push!(elset_names_inp, m[1])
                        elset_count += 1
                        if elset_count > 1
                            push!(elements_inp, hcat(elements...)')
                            elements = Any[]
                        end
                        continue
                    end
                        
                    if start_reading && occursin("**", line)
                        push!(elements_inp, hcat(elements...)')
                        break 
                    end
                        
                    if start_reading
                        element_line_str = split(strip(line), ',')
                        push!(elements, parse.(Int, element_line_str))
                    end
                end
            end

            # Elsets (lists of element ids, if present we will add to elements_inp)
            if parse_elsets
                seekstart(f)
                elements_all = vcat(elements_inp...)
                id = sortperm(elements_all[:, 1])
                elements_all = elements_all[id, :] # sorted by element id
                
                start_reading = false
                elements = Any[]
                for line in eachline(f)
                    if occursin("*ELSET, ELSET=", line)  
                        start_reading = true
                        m = match(r"ELSET=(.*)", line)
                        push!(elset_names_inp, m[1])

                        if !isempty(elements)
                            elements_set = Set(elements) # for speed
                            id_rows = [x in elements_set for x in elements_all[:, 1]]
                            rows_elset = elements_all[id_rows, :]
                            push!(elements_inp, rows_elset)
                            elements = Any[]
                            println("Elset $m[1] completed.")
                        end
                        continue
                    end
                        
                    if start_reading && occursin("**", line)
                        push!(elements_inp, hcat(elements...)')
                        break 
                    end
                        
                    if start_reading
                        # Use filter to remove any empty strings from element_line_str before parsing:
                        element_line_str = filter(x -> !isempty(x), split(strip(line), ','))
                        append!(elements, parse.(Int, element_line_str))
                    end
                end
            end

            # Nodesets
            seekstart(f)
            start_reading = false 
            nset_count = 0
            nodes = Any[]
            for line in eachline(f)
                if occursin("*NSET, NSET=", line)   
                    start_reading = true
                    m = match(r"NSET=(.*)", line)
                    #println(join(["Nset name: ", m[1], "\n"]))
                    push!(nset_names_inp, m[1])
                    nset_count += 1
                    if nset_count > 1
                        push!(nset_inp, vcat(nodes...))
                        nodes = Any[]
                    end
                    continue
                end
                    
                if start_reading && occursin("**", line)
                    push!(nset_inp, vcat(nodes...))
                    break 
                end
                    
                if start_reading
                    nodes_line_str = split(strip(line), ',')
                    nodes_line_str = filter(!isempty, nodes_line_str) # get rid of empty strings after comma
                    push!(nodes, parse.(Int, nodes_line_str))
                end
            end
        end
        nodes_inp = hcat(node_data...)' 
        
        return nodes_inp, elements_inp, elset_names_inp, nset_inp, nset_names_inp
    end


end