#=
Script to compute fault tractions at wilmington
=#

# Cd into current directory and activate the Julia environment in this directory
cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__)

# startup
using Revise
using GLMakie
GLMakie.activate!()
using Colors
using ColorSchemes
using LaTeXStrings
using DelimitedFiles
using Statistics

lib_dir = "C:/Users/lsalo/julia/jl-tools/"
#lib_dir = "/Users/lluis/Documents/julia/jl-tools/"
includet(joinpath(lib_dir, "modules", "abaqus.jl"))
using .Abaqus
includet(joinpath(lib_dir, "modules", "arrayOps.jl"))
using .ArrayOps
includet(joinpath(lib_dir, "modules", "MakiePolyPlot.jl"))
using .MakiePolyPlot

# Directories
dir_inp     = "C:/Users/lsalo/Wilmington/AbaqusGrids/"
dir_abaqus  = "C:/Users/lsalo/abaqus/simulation_data/wilmington_old/output/plastic/extracted/U/"
figs_dir    = "C:/Users/lsalo/julia/graphics/figs/wilmington/2024-05-ERL/" 
#dir_inp = "/Users/lluis/Dropbox (MIT)/home_copy/abaqus/wilmington/grids/"
#dir_abaqus = "/Users/lluis/Dropbox (MIT)/home_copy/abaqus/wilmington/output/elastic/extracted/S/"
#figs_dir = "/Users/lluis/Documents/julia/graphics/figs/wilmington/"


# 1. Nodal data (connectivity and displacements)
# Get connectivity fron .inp grid: We need to get the *face* connectivity of node sets
fn_inp  = "WilmingtonMesh_AllFaults_Rotated_ENTIRE_DOMAIN_Coarse.inp"
#nsets_to_export = ["TopBC", "TopRangerSurface"]
#nsets_to_export = ["ComptonLosAlamitos_a", "ComptonLosAlamitos_b", "ComptonLosAlamitos_c", "Newport", "PalosVerdes", "Tear_Fault_2", "Tear_Fault_4", "Tear_Fault_5", "Tear_Fault_6", "Tear_Fault_8a", "Tear_Fault_8b", "Tear_Fault_9", "Tear_Fault_10a", "Tear_Fault_10b", "Tear_Fault_11", "Tear_Fault_13a", "Tear_Fault_13b", "Tear_Fault_15", "Tear_Fault_19", "Tear_Fault_21a", "Tear_Fault_21b", "Tear_Fault_26", "Tear_Fault_27", "Wilmington"]
nsets_to_export = ["Newport"]
nodes_inp, elements_inp, elset_names_inp, nset_inp, nset_names_inp = get_connectivity_from_inp(dir_inp, fn_inp, "tetrahedral")
face_list = get_face_connectivity_nodesets(nsets_to_export, nset_names_inp, nset_inp, elements_inp[1], "tetrahedral")

# Simulation data: get values at the nodes from abaqus
nodal_disp = Vector{}(undef, length(nsets_to_export))
tstep = "step_300"
for nset in eachindex(nsets_to_export)
    data = readdlm(join([dir_abaqus, nsets_to_export[nset], "/", tstep, "_nodes.dat"]), ',', Float64, '\n', skipstart=1)
    nodal_disp[nset] = Matrix{Union{Int,Float64}}(data)
    nodal_disp[nset][:, 1] .= Int.(nodal_disp[nset][:, 1])
end

# 2. Element data (connectivity and displacements)
fn_inp  = "WilmingtonMesh_AllFaults_Rotated_Coarse.inp"
elsets_to_export = ["Ranger"]
_, elements_inp, elset_names_inp, nset_inp, nset_names_inp = get_connectivity_from_inp(dir_inp, fn_inp, "tetrahedral")
face_list_el = get_face_connectivity_elsets(elsets_to_export, elset_names_inp, elements_inp, "tetrahedral")

nodal_disp_el = Vector{}(undef, length(elsets_to_export))
for elset in eachindex(elsets_to_export)
    data = readdlm(join([dir_abaqus, elsets_to_export[elset], "/", tstep, "_nodes.dat"]), ',', Float64, '\n', skipstart=1)
    nodal_disp_el[elset] = Matrix{Union{Int,Float64}}(data)
    nodal_disp_el[elset][:, 1] .= Int.(nodal_disp_el[elset][:, 1])
end

for elset in eachindex(elsets_to_export)
    push!(face_list, face_list_el[elset])
    push!(nodal_disp, nodal_disp_el[elset])
end

# 3. Plot
# For poly we need: 1) nodal coords, 2) nodal values, 3) face connectivity, starting at 1
nsets = length(nsets_to_export)+length(elsets_to_export)
coords_map = Vector{}(undef, nsets)
vals_nodes = Vector{}(undef, nsets)
face_list_map = Vector{}(undef, nsets)
for set in 1:nsets
    unique_ids = sort(unique(face_list[set]))      
    coords_map[set] = nodes_inp[in.(nodes_inp[:, 1], Ref(unique_ids)), 2:end]
    nodal_disp_plot = sortslices(nodal_disp[set], dims=1, by=x->x[1])  
    vals_nodes[set] = nodal_disp_plot[:, end] # U3 [m]   
    face_list_map[set], _ = renumber_array_1_to_n(face_list[set])
end

# Make plot
labl = L"$u_z$ [m]"
mycmap = cgrad(:lapaz)
fig_descriptor = (cmap=mycmap, cbar_label=labl, fig_size = (1500, 1200), plot_height = 1500, cbar_height = 400, plot_azi = 1*pi, plot_elev = 0.08*pi, plot_linew=0.1, xlms=(nothing, nothing), ylms=(nothing, nothing), save_fig = false, save_dir = figs_dir, title = "u_3_ranger_Newport", target_dpi=400, target_size_in=(15,12))
plotNodalValues(face_list_map, coords_map, vals_nodes, fig_descriptor)