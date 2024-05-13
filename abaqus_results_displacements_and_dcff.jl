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
dir_abaqus  = "C:/Users/lsalo/abaqus/simulation_data/wilmington/output/plastic/extracted/U/"
figs_dir    = "C:/Users/lsalo/julia/graphics/figs/wilmington/test/" 
#dir_inp = "/Users/lluis/Dropbox (MIT)/home_copy/abaqus/wilmington/grids/"
#dir_abaqus = "/Users/lluis/Dropbox (MIT)/home_copy/abaqus/wilmington/output/elastic/extracted/S/"
#figs_dir = "/Users/lluis/Documents/julia/graphics/figs/wilmington/" 


# 1. Get connectivity fron .inp grid
# We need to get the *face* connectivity of element sets or node sets
fn_inp  = "WilmingtonMesh_AllFaults_Rotated_ENTIRE_DOMAIN_Coarse.inp"
nsets_to_export = ["TopBC", "TopRangerSurface"]
nodes_inp, elements_inp, elset_names_inp, nset_inp, nset_names_inp = get_connectivity_from_inp(dir_inp, fn_inp, "tetrahedral")
face_list = get_face_connectivity_nodesets(nsets_to_export, nset_names_inp, nset_inp, elements_inp[1], "tetrahedral")


# 2. Simulation data: get values at the nodes from abaqus
nodal_disp = Vector{}(undef, length(nsets_to_export))
tstep = "step_300"
for nset in eachindex(nsets_to_export)
    data = readdlm(join([dir_abaqus, nsets_to_export[nset], "/", tstep, "_nodes.dat"]), ',', Float64, '\n', skipstart=1)
    nodal_disp[nset] = Matrix{Union{Int,Float64}}(data)
    nodal_disp[nset][:, 1] .= Int.(nodal_disp[nset][:, 1])
end

# 3. Plot
# For poly we need: 1) nodal coords, 2) nodal values, 3) face connectivity, starting at 1
coords_map = Vector{}(undef, length(nsets_to_export))
vals_nodes = Vector{}(undef, length(nsets_to_export))
face_list_map = Vector{}(undef, length(nsets_to_export))
for nset in eachindex(nsets_to_export)
    unique_ids = sort(unique(face_list[nset]))      
    coords_map[nset] = nodes_inp[in.(nodes_inp[:, 1], Ref(unique_ids)), 2:end]
    nodal_disp_plot = sortslices(nodal_disp[nset], dims=1, by=x->x[1])  
    vals_nodes[nset] = nodal_disp_plot[:, end] # U3 [m]   
    face_list_map[nset], _ = renumber_array_1_to_n(face_list[nset])
end

# Make plot
labl = L"$u_z$ [m]"
mycmap = cgrad(:vik)
fig_descriptor = (cmap=mycmap, cbar_label=labl, fig_size = (800, 600), plot_height = 700, cbar_height = 200, plot_azi = 1.35*pi, plot_elev = 0.2*pi, plot_linew=0.1, xlms=(nothing, nothing), ylms=(nothing, nothing), save_fig = false, save_dir = figs_dir, title = "u_3", target_dpi=350, target_size_in=(8,6))
plotNodalValues(face_list_map, coords_map, vals_nodes, fig_descriptor)