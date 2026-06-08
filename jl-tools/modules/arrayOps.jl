module ArrayOps
    """
    Common operations on arrays
    """

    using StatsBase
    using Statistics
    using RegularExpressions
    using Distances
    using DelimitedFiles

    export resample_matrix, fill_zero_rows, reduce_matrix_with_repeated_ids, find_ids_from_pairwise_distance, renumber_array_1_to_n

    function resample_matrix(data::Matrix{<:Number}, freq::Int)
        """
        Resample a 2D matrix to a new size using mean value of the original data.
        The new size is determined by the frequency parameter, which specifies how many columns to average over.
        """
        
        new_data = Matrix{Float64}(undef, size(data, 1), ceil(Int, size(data, 2) / freq))
        i = 1
        for tstep in 1:size(new_data, 2) 
            new_data[:, tstep] = mean(data[:, i:i+(freq-1)], dims=2)
            i += freq
        end

        return new_data
    end

    function fill_zero_rows(data::Matrix{Float64}; method::String="nearest_average")
        """
        Fill zero rows in a 2D matrix. Rows with (0.0, 0.0) are replaced by:
            - Method="nearest_average": the average of the nearest nonzero rows above and below. If only one is found, that one is used.
        """
        new_data = copy(data)  # work on a copy to avoid touching original
        N = size(new_data, 1)
        for i in 1:N
            if new_data[i, 1] == 0.0 && new_data[i, 2] == 0.0
                if method=="nearest_average"
                    # Find nearest nonzero row above
                    above = nothing
                    for j in (i-1):-1:1
                        if new_data[j, 1] != 0.0 || new_data[j, 2] != 0.0
                            above = new_data[j, :]
                            break
                        end
                    end

                    # Find nearest nonzero row below
                    below = nothing
                    for j in (i+1):N
                        if new_data[j, 1] != 0.0 || new_data[j, 2] != 0.0
                            below = new_data[j, :]
                            break
                        end
                    end

                    # Average available neighbors
                    if !isnothing(above) && !isnothing(below)
                        new_data[i, :] .= (above .+ below) ./ 2
                    elseif !isnothing(above)
                        new_data[i, :] .= above
                    elseif !isnothing(below)
                        new_data[i, :] .= below
                    else
                        @warn "No valid neighbors found for row $i. Leaving as (0.0, 0.0)"
                    end
                
                else
                    error("Unsupported method: $method")
                end
            end
        end

        return new_data
    end
    
    function reduce_matrix_with_repeated_ids(data_matrix::AbstractMatrix, reduction_method::String; out_matrix_type="Int")
        """
        Note: Currently limited to mean value of repeated entries.

        Example of data matrix:
            stress_vals: n by 7 matrix with first column being [repeated] identifier (e.g. node number), s11, s22, s33, s12, s13, s23
        """

        unique_ids = sort(unique(data_matrix[:, 1]))
        if out_matrix_type == "Int"
            reduced_matrix = zeros(Int, length(unique_ids), size(data_matrix)[2])
        elseif out_matrix_type == "Int-Float"
            reduced_matrix = Matrix{Union{Int, Float64}}(undef, length(unique_ids), size(data_matrix)[2])
        else
            error("Unsupported output matrix type: $out_matrix_type")
        end

        for (i, id) in enumerate(unique_ids)
            reduced_matrix[i, 1] = id
            id_rows = data_matrix[:, 1] .== id
            
            if reduction_method == "mean"
                reduced_vals = mean(data_matrix[id_rows, 2:end], dims=1)
            elseif reduction_method == "sum"
                reduced_vals = sum(data_matrix[id_rows, 2:end], dims=1)
            end

            if all(isinteger, data_matrix)
                reduced_matrix[i, 2:end] = Int.(round.(reduced_vals))
            else
                reduced_matrix[i, 2:end] = reduced_vals
            end
        end

        return reduced_matrix

    end

    function find_ids_from_pairwise_distance(A::Matrix, B::Matrix; info_to_save=[])
        """
        Finds closest ids in B to every row in A using pairwise Euclidean distances
        Assumes A and B have dimensions of (n x 2/3) and (m x 2/3), respectively
        """
        if size(A, 2) == size(B, 2)
            r = pairwise(Euclidean(), A', B')
            ids = [cell_id[2] for cell_id in argmin(r, dims=2)]

            if !isempty(info_to_save)
                writedlm(join([info_to_save[1], info_to_save[2], ".dat"]), ids)
            end

            return ids

        else
            error("Dimension mismatch. $name1 and $name2 must have same number of columns.")
        end
    end


    function renumber_array_1_to_n(T::Union{Matrix, Vector})
        """
        Used to renumber nodes, which may not be numbered starting from 1 in the .inp grid, but need to start from 1 for MRST and plotting.
        """
        unique_ids = sort(unique(T))
        new_ids = 1:length(unique_ids)
        T_new = similar(T)
        map_ini_to_new = Dict{Int, Int}()
        for n in new_ids
            k, v = unique_ids[n], new_ids[n]
            T_new[T .== k] .= v
            map_ini_to_new[k] = v
        end

        return T_new, map_ini_to_new

    end


end
