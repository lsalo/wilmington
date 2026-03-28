#=
Script to visualize displacements on grid at wilmington
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

lib_dir = "../jl-tools/"
includet(joinpath(lib_dir, "modules", "abaqus.jl"))
using .Abaqus
includet(joinpath(lib_dir, "modules", "arrayOps.jl"))
using .ArrayOps
includet(joinpath(lib_dir, "modules", "MakiePolyPlot.jl"))
using .MakiePolyPlot

# Directories
sim_name = "plastic-120"
platform = Sys.ARCH
if platform == :x86_64 
    dir_inp  = "C:/Users/lsalo/Wilmington/AbaqusGrids/"
    sim_dir  = "C:/Users/lsalo/abaqus/simulation_data/wilmington/"
    figs_dir = "C:/Users/lsalo/julia/graphics/figs/wilmington/test/"
elseif platform == :aarch64
    dir_inp  = "/Users/lluis/Dropbox (MIT)/home_copy/abaqus/wilmington/AbaqusGrids/"
    sim_dir  = "/Users/lluis//Dropbox (MIT)/home_copy/abaqus/simulation_data/wilmington/"
    figs_dir = "/Users/lluis/Dropbox (MIT)/home_copy/julia/graphics/wilmington/test/"
    data_dir = "/Users/lluis/Dropbox (MIT)/home_copy/julia/data/wilmington/"
    dir_out = "/Users/lluis/Dropbox (MIT)/home_copy/julia/output/wilmington/output_data/displacement/" 
end
sim_dir = sim_dir * sim_name * "/extracted/U/"

# 1. Nodal data (connectivity and displacements)
# Get connectivity from .inp grid: We need to get the *face* connectivity of node sets
fn_inp  = "WilmingtonMesh_AllFaults_Rotated_ENTIRE_DOMAIN_Coarse.inp"
#nsets_to_export = ["TopBC", "TopRangerSurface"]
#nsets_to_export = ["Wilmington"]
#nsets_to_export = ["ComptonLosAlamitos_a", "ComptonLosAlamitos_b", "ComptonLosAlamitos_c", "Newport", "PalosVerdes", "Tear_Fault_2", "Tear_Fault_4", "Tear_Fault_5", "Tear_Fault_6", "Tear_Fault_8a", "Tear_Fault_8b", "Tear_Fault_9", "Tear_Fault_10a", "Tear_Fault_10b", "Tear_Fault_11", "Tear_Fault_13a", "Tear_Fault_13b", "Tear_Fault_15", "Tear_Fault_19", "Tear_Fault_21a", "Tear_Fault_21b", "Tear_Fault_26", "Tear_Fault_27", "Wilmington"]
#nsets_to_export = ["Tear_Fault_2", "Tear_Fault_4", "Tear_Fault_5", "Tear_Fault_6", "Tear_Fault_8a", "Tear_Fault_8b", "Tear_Fault_9", "Tear_Fault_10a", "Tear_Fault_10b", "Tear_Fault_11", "Tear_Fault_13a", "Tear_Fault_13b", "Tear_Fault_15", "Tear_Fault_19", "Tear_Fault_21a", "Tear_Fault_21b", "Tear_Fault_26", "Tear_Fault_27"]
nsets_to_export = []
if ~isempty(nsets_to_export)
    nodes_inp, elements_inp, elset_names_inp, nset_inp, nset_names_inp = get_connectivity_from_inp(dir_inp, fn_inp, "tetrahedral", false)
    face_list = get_face_connectivity_nodesets(nsets_to_export, nset_names_inp, nset_inp, elements_inp[1], "tetrahedral")
else
    face_list = []
end

# Simulation data: get values at the nodes from abaqus
nodal_disp = Vector{}(undef, length(nsets_to_export))
#tstep = "step_468"
tstep = "step_89"
#tstep = "step_117"
#tstep = "step_157"
#tstep = "step_337"
for nset in eachindex(nsets_to_export)
    data = readdlm(join([sim_dir, nsets_to_export[nset], "/", tstep, "_nodes.dat"]), ',', Float64, '\n', skipstart=1)
    nodal_disp[nset] = Matrix{Union{Int,Float64}}(data)
    nodal_disp[nset][:, 1] .= Int.(nodal_disp[nset][:, 1])
end

# 2. Element data (connectivity and displacements)
fn_inp  = "WilmingtonMesh_AllFaults_Rotated_Coarse.inp"
elsets_to_export = ["OverburdenPico"] #OverburdenPico
_, elements_inp, elset_names_inp, nset_inp, nset_names_inp = get_connectivity_from_inp(dir_inp, fn_inp, "tetrahedral")
face_list_el = get_face_connectivity_elsets(elsets_to_export, elset_names_inp, elements_inp, "tetrahedral")

nodal_disp_el = Vector{}(undef, length(elsets_to_export))
for elset in eachindex(elsets_to_export)
    data = readdlm(join([sim_dir, elsets_to_export[elset], "/", tstep, "_nodes.dat"]), ',', Float64, '\n', skipstart=1)
    nodal_disp_el[elset] = Matrix{Union{Int,Float64}}(data)
    nodal_disp_el[elset][:, 1] .= Int.(nodal_disp_el[elset][:, 1])
end

for elset in eachindex(elsets_to_export)
    push!(face_list, face_list_el[elset])
    push!(nodal_disp, nodal_disp_el[elset])
end

# Not rotated
fn_inp2  = "WilmingtonMesh_AllFaults_Coarse.inp" # Not rotated grid, with actual UTM coords (zone 11)
nodes_inp2, _, _, _, _ = get_connectivity_from_inp(dir_inp, fn_inp2, "tetrahedral", false)

# 3. Plot
# For poly we need: 1) nodal coords, 2) nodal values, 3) face connectivity, starting at 1
nsets = length(nsets_to_export)+length(elsets_to_export)
coords_map = Vector{}(undef, nsets)
vals_nodes = Vector{}(undef, nsets)
face_list_map = Vector{}(undef, nsets)
for set in 1:nsets
    unique_ids = sort(unique(face_list[set]))      
    coords_map[set] = nodes_inp2[in.(nodes_inp2[:, 1], Ref(unique_ids)), 2:end]
    nodal_disp_plot = sortslices(nodal_disp[set], dims=1, by=x->x[1])  
    #vals_nodes[set] = nodal_disp_plot[:, end] # U3 [m]   
    vals_nodes[set] = sqrt.(nodal_disp_plot[:, 2].^2 + nodal_disp_plot[:, 3].^2) # U_hzntal [m]  
    face_list_map[set], _ = renumber_array_1_to_n(face_list[set])
end

# Make plot
#labl = L"$u_\text{z}$ [m]"
labl = L"$||\mathbf{u}_h||$ [m]"
mycmap = cgrad(:bone_1, rev=true)

# Fig descriptor
fig_descriptor = (cmap=mycmap, cbar_label=labl, fig_size = (1500, 1200), plot_height = 1000, cbar_height = 400, plot_azi = -0.5*pi, plot_elev = 0.5*pi, plot_linew=0.1, xlms=(minimum(coords_map[1][:,1]), maximum(coords_map[1][:,1])), ylms=(minimum(coords_map[1][:,2]), maximum(coords_map[1][:,2])), save_fig = true, save_dir = figs_dir, title = "u_h_ovbPico_plastic-120_t89", target_dpi=400, target_size_in=(15,12))

plotNodalValues(face_list_map, coords_map, vals_nodes, fig_descriptor)

#Line
coast_data = readdlm(join([data_dir, "cal_outline.dat"]))
cst = [coast_data[coast_data[:, 3] .== i, :] for i in 1:18] # third column is "category"
coastline = [cst[6], cst[17], cst[18]]
#coastline = [hcat(coastline[:, 2], coastline[:, 1], coastline[:, 3])]
# Plotting
plotNodalValues(face_list_map, coords_map, vals_nodes, fig_descriptor, line_data=coastline)

# Save coordinates and horizontal disp magnitudes
writedlm(dir_out * sim_name * "/nodal_coords_and_u_h_step_89.dat", hcat(coords_map[1], vals_nodes[1])) 