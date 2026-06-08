module Gmsh
"""
This module performs operations related to Gmsh: https://gmsh.info/
"""

export get_unit_names_and_element_ids_from_msh

function get_unit_names_and_element_ids_from_msh(file_path)
    """
    Open .msh file (gmsh mesh), and get physical names (of facies or units) and element ids corresponding to each unit.

    NOTE: We currently assume that the Physical Name of units is "Facies 1", "Facies 2", etc. Function needs updating if that is not the case.
    """
    # Read file and get array with group of each element
    file = open(file_path, "r")
    inside_elements = false
    inside_names = false
    groups = Int64[]
    names = String[]
    for line in eachline(file)
        # Get facies names
        if occursin(r"Facies \d+", line)
            name = match(r"Facies \d+", line).match
            push!(names, name) 
        end

        # Get element groups
        if occursin("\$Elements", line)
            inside_elements = true
            continue  
        end
        if inside_elements
            parts = split(line)
            if length(parts) >= 2
                push!(groups, parse(Float64, parts[4]))
            end
        end
    end
    close(file)

    return names, groups
end

end