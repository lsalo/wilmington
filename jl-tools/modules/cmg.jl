module Cmg
"""
Performs operations to generate an input file for CMG
"""

export write_grid_file_for_cmg, process_indices_blockgroup

function write_grid_file_for_cmg(nx, ny, nz, dx, dy, dz, dtop, output_dir, fn_grid, fn_bg)
    output_file = open(joinpath(output_dir, fn_grid), "w")
    str_bg = join(["INCLUDE ", fn_bg])

    println(output_file, "GRID VARI $nx $ny $nz")
    println(output_file, "KDIR DOWN")
    println(output_file, "DI IVAR")
    println(output_file, " $nx*$dx")
    println(output_file, "DJ JVAR")
    println(output_file, " $ny*$dy")
    println(output_file, "DK KVAR")
    println(output_file, " $nz*$dz")
    println(output_file, "DTOP")
    println(output_file, " $nx*$dtop")
    println(output_file, " ")
    println(output_file, str_bg)
    close(output_file)
end

function write_range(current_range_i, current_range_k, name, firstln, output_file)
    bg = "BLOCKGROUP '$name' IJK "
    if !firstln
        bg = " " ^ length(bg)
    end
    if length(current_range_i) == length(current_range_k)
        println(output_file, bg, join(current_range_i), " 1 ", join(current_range_k), "    1")
    elseif length(current_range_i) > length(current_range_k)
        println(output_file, bg, join(current_range_i[1]), ":", join(current_range_i[end]), " 1 ", join(current_range_k), "    1")
     elseif length(current_range_i) < length(current_range_k)
        println(output_file, bg, join(current_range_i), " 1 ", join(current_range_k[1]), ":", join(current_range_k[end]), "    1")
    end
end

function process_indices_blockgroup(name, indices, output_file)
    """
    TBW
    """
    current_range_i, current_range_k = Int[], Int[]
    inrange_i, inrange_k = false, false
    firstln = true
    for (n, entry) in enumerate(indices)
        i, k = entry
        if n < length(indices)
            next_i, next_k = indices[n+1]
        end
        if n==length(indices) || (inrange_i && i+1 != next_i || inrange_k && k+1 != next_k || i+1 == next_i && k+1 == next_k || i+1 != next_i && k+1 != next_k)
            if inrange_i
                push!(current_range_i, i)
            elseif inrange_k
                push!(current_range_k, k)
            else
                push!(current_range_i, i)
                push!(current_range_k, k)
            end
            write_range(current_range_i, current_range_k, name, firstln, output_file)
            current_range_i = Int[]
            current_range_k = Int[]
            inrange_i, inrange_k = false, false
            firstln = false
        elseif i+1 == next_i && k == next_k
            push!(current_range_i, i)
            if !inrange_i
                push!(current_range_k, k)
                inrange_i = true
            end
        elseif i == next_i && k+1 == next_k
            push!(current_range_k, k)
            if !inrange_k
                push!(current_range_i, i)
                inrange_k = true
            end
        end
    end
    println(output_file, " ")  # add blank line after blockgroup
end

end