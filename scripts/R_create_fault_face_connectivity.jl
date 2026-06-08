#=
Script to create tear fault input files for R-script
=#

cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__)

# startup
using Revise
using DelimitedFiles
using Statistics
using LinearAlgebra
using Printf

lib_dir = "C:/Users/lsalo/julia/jl-tools/"
#lib_dir = "/Users/lluis/Documents/julia/jl-tools/"
includet(joinpath(lib_dir, "modules", "abaqus.jl"))
using .Abaqus
includet(joinpath(lib_dir, "modules", "arrayOps.jl"))
using .ArrayOps

# Directories
dir_inp = "C:/Users/lsalo/Wilmington/AbaqusGrids/"
dir_out = "C:/Users/lsalo/julia/output/wilmington/R/"

# Get connectivity fron .inp grid
fn_inp  = "WilmingtonMesh_AllFaults_Coarse.inp"
nodes_inp, elements_inp, elset_names_inp, nset_inp, nset_names_inp = get_connectivity_from_inp(dir_inp, fn_inp, "tetrahedral")
nsets_to_export = ["Tear_Fault_2", "Tear_Fault_4", "Tear_Fault_5", "Tear_Fault_6", "Tear_Fault_8a", "Tear_Fault_8b", "Tear_Fault_9", "Tear_Fault_10a", "Tear_Fault_10b", "Tear_Fault_11", "Tear_Fault_13a", "Tear_Fault_13b", "Tear_Fault_15", "Tear_Fault_19", "Tear_Fault_21a", "Tear_Fault_21b", "Tear_Fault_26", "Tear_Fault_27"]
face_list = get_face_connectivity_nodesets(nsets_to_export, nset_names_inp, nset_inp, vcat(elements_inp...), "tetrahedral")

# Prepare quantities and write file
for n in eachindex(nsets_to_export)
    nset_id = findfirst(nset_names_inp .== nsets_to_export[n])

    # Nodes
    id_nodes = findall(in(nset_inp[nset_id]), nodes_inp[:,1])
    nset_nodes = nodes_inp[id_nodes, :]
    n_nodes = size(nset_nodes)[1]
    nset_nodes[:, 1] = 1:n_nodes

    # Face connectivity
    n_faces = size(face_list[n])[1]
    face_list_out = Matrix{Union{Int64, Float64}}(undef, n_faces, 4)
    face_list_out[:, 1] = 1:n_faces
    face_list_out[:, 2:4], _ = renumber_array_1_to_n(face_list[n])   
    
    # Write file
    open(join([dir_out, nsets_to_export[n], ".dat"]), "w") do f
        writedlm(f, [n_nodes n_faces], " ")

        for row in eachrow(nset_nodes) # Format: 1st col integers, 2nd to 4th floats (avoid scientific notation)
            formatted_row = [@sprintf("%d", Int(row[1]))]
            append!(formatted_row, [@sprintf("%.1f", x) for x in row[2]])
            append!(formatted_row, [@sprintf("%.0f", x) for x in row[3]])
            append!(formatted_row, [@sprintf("%.2f", x) for x in row[4]])
            writedlm(f, [formatted_row], " ")
        end

        writedlm(f, face_list_out, " ")
    end

end