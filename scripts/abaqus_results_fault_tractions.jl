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
using LinearAlgebra
using HDF5

lib_dir = "../jl-tools/"
includet(joinpath(lib_dir, "modules", "faultOperations.jl"))
using .FaultOperations
includet(joinpath(lib_dir, "modules", "abaqus.jl"))
using .Abaqus
includet(joinpath(lib_dir, "modules", "arrayOps.jl"))
using .ArrayOps
includet(joinpath(lib_dir, "modules", "fileio.jl"))
using .Fileio
includet(joinpath(lib_dir, "modules", "makiePolyPlot.jl"))
using .MakiePolyPlot

# Directories
sim_name = "plastic-120"
platform = Sys.ARCH
if platform == :x86_64 
    dir_inp = "C:/Users/lsalo/Wilmington/AbaqusGrids/"
    sim_dir = "C:/Users/lsalo/abaqus/simulation_data/wilmington/"
    dir_out = "C:/Users/lsalo/julia/output/wilmington/fault_quantities/"
    figs_dir = "C:/Users/lsalo/julia/graphics/figs/wilmington/test/"
elseif platform == :aarch64  # MBP
    dir_inp = "/Users/lluis/Dropbox (MIT)/home_copy/abaqus/wilmington/AbaqusGrids/"
    sim_dir = "/Users/lluis/Dropbox (MIT)/home_copy/abaqus/simulation_data/wilmington/"
    dir_out = "/Users/lluis/Dropbox (MIT)/home_copy/julia/output/wilmington/output_data/fault_quantities/" 
    figs_dir = "/Users/lluis/Dropbox (MIT)/home_copy/julia/graphics/wilmington/test/" 
    data_dir = "/Users/lluis/Dropbox (MIT)/home_copy/julia/data/wilmington/"
end
sim_dir = sim_dir * sim_name * "/extracted/S/"
dir_out = dir_out * sim_name

# Get connectivity fron .inp grid
fn_inp  = "WilmingtonMesh_AllFaults_Rotated_ENTIRE_DOMAIN_Coarse.inp"
nodes_inp, elements_inp, elset_names_inp, nset_inp, nset_names_inp = get_connectivity_from_inp(dir_inp, fn_inp, "tetrahedral")

# Get nodeset faces - for now tetrahedra only
#nsets_to_export = ["Tear_Fault_2", "Tear_Fault_4", "Tear_Fault_5", "Tear_Fault_6", "Tear_Fault_8a", "Tear_Fault_8b", "Tear_Fault_9", "Tear_Fault_10a", "Tear_Fault_10b", "Tear_Fault_11", "Tear_Fault_13a", "Tear_Fault_13b", "Tear_Fault_15", "Tear_Fault_19", "Tear_Fault_21a", "Tear_Fault_21b", "Tear_Fault_26", "Tear_Fault_27"]
#nsets_to_export = ["Wilmington"]
nsets_to_export = ["PicoSurface"]
#nsets_to_export = ["TopRangerSurface", ]
#nsets_to_export = ["Tear_Fault_4", "Tear_Fault_9"]
#nsets_to_export = ["Tear_Fault_4", "Tear_Fault_9", "Tear_Fault_10a", "Tear_Fault_10b","Tear_Fault_11"]

#face_list = get_face_connectivity_nodesets(nsets_to_export, nset_names_inp, nset_inp, elements_inp[1], "tetrahedral")
face_list = get_face_connectivity_nodesets_C3D4(nsets_to_export, nset_names_inp, nset_inp, elements_inp[1])

# Face centroids and normals
face_unit_normals, face_centroids = face_normals_and_centroids(face_list, nodes_inp)

# Calculate shear, normal tractions at each node, and dcff
mu = 0.6
nodal_tractions = [Vector{Any}() for _ in 1:length(nsets_to_export)]
#dcff = Vector{Matrix{Union{Int64,Float64}}}(undef, length(nsets_to_export))
dcff = [Vector{Any}() for _ in 1:length(nsets_to_export)]
cff = Vector{Matrix{Float64}}(undef, length(nsets_to_export))
cff0 = [Vector{Any}() for _ in 1:length(nsets_to_export)]
compute_tractions = false
for nset in eachindex(nsets_to_export)
    abaqus_fn_nod = readdir(join([sim_dir, nsets_to_export[nset]]))
    abaqus_fn_nod = sort(abaqus_fn_nod, by=sort_steps)
    cff_cols = Vector{Vector{Float64}}()
    for (i, fn) in enumerate(abaqus_fn_nod)
        fpath = joinpath(sim_dir, nsets_to_export[nset], fn)
        ext = splitext(fn)[2]
        
        if compute_tractions
            if ext == ".dat"
                nodal_stress = readdlm(fpath, ',', Float64, '\n', skipstart=1)
            elseif ext == ".h5"
                nodal_stress = h5open(fpath, "r") do f
                                    hcat(read(f, "NODE_LABEL"), read(f, "Sxx"), read(f, "Sxy"), read(f, "Sxz"), read(f, "Syy"), read(f, "Syz"), read(f, "Szz"))
                end
                nodal_stress = reduce_matrix_with_repeated_ids(nodal_stress, "mean")
            else
                error("Unsupported file extension: $ext for file $fn")
            end

            dir_path = joinpath([dir_out, nsets_to_export[nset]])
            fname = join(["tractions_", fn])
            push!(nodal_tractions[nset], tractions_from_sij_and_n(nodal_stress, nodes_inp, face_unit_normals[nset], face_centroids[nset], dir_path, fname))           

            print(join([fn " completed.\n"]))

        else # we already have them, so we just read
            if ext == ".dat"
                data = readdlm(joinpath([dir_out, nsets_to_export[nset], join(["tractions_", fn])]), '\t', Float64, '\n')
            elseif ext == ".h5"
                #TBD
            else
                error("Unsupported file extension: $ext for file $fn")
            end
            push!(nodal_tractions[nset], data)
        end

        if i == 1
            cff0[nset] = nodal_tractions[nset][1][:, end] .+ mu*nodal_tractions[nset][1][:, 2]
        end

        if i > 1
            dtau = nodal_tractions[nset][i][:, end] .- nodal_tractions[nset][1][:, end]
            dsn = nodal_tractions[nset][i][:, 2] .- nodal_tractions[nset][1][:, 2]
            cff_i = nodal_tractions[nset][i][:, end] .+ mu*nodal_tractions[nset][i][:, 2]
            push!(cff_cols, cff_i) 
            push!(dcff[nset], hcat(nodal_tractions[nset][1][:,1], cff_i .- cff0[nset], dsn, dtau))
        end
    end
    
    cff[nset] = reduce(hcat, cff_cols)
    print(join(["Fault tractions for nodeset ", nsets_to_export[nset], " completed.\n"]))
end


# load nodal disp value file for area filtering
node_subset_for_tseries = true
if node_subset_for_tseries
    pth_fn = "/Users/lluis/Dropbox (MIT)/home_copy/julia/output/wilmington/output_data/displacement/" * sim_name * "/nodal_coords_and_u_h_step_89.dat"
    nodal_coords_uh = readdlm(pth_fn);
    id_above_uh_thresh = nodal_coords_uh[:, end] .>= 1.0
    x_min, x_max = minimum(nodal_coords_uh[id_above_uh_thresh, 1]), maximum(nodal_coords_uh[id_above_uh_thresh, 1])
    y_min, y_max = minimum(nodal_coords_uh[id_above_uh_thresh, 2]), maximum(nodal_coords_uh[id_above_uh_thresh, 2])
    id_nodes_in_thresh = vec(all(hcat(coords_map[1][:,1] .>= x_min,  coords_map[1][:,1] .<= x_max, coords_map[1][:,2] .>= y_min,  coords_map[1][:,2] .<= y_max), dims=2))
else
    # use nsets_to_export to find number of nsets and export trues for number of nodes
end

# Averages for each fault
n_steps = length(dcff[1])
dcff_avg = [Vector{Float64}(undef, n_steps) for _ in 1:length(nsets_to_export)]
dcff_95 = [Vector{Float64}(undef, n_steps) for _ in 1:length(nsets_to_export)]
cff_avg = [Vector{Float64}(undef, n_steps) for _ in 1:length(nsets_to_export)]
cff_95 = [Vector{Float64}(undef, n_steps) for _ in 1:length(nsets_to_export)]
for nset in eachindex(nsets_to_export)
    for i in 1:n_steps
        dcff_avg[nset][i] = mean(dcff[nset][i][id_nodes_in_thresh, 2])
        dcff_95[nset][i] = quantile(dcff[nset][i][id_nodes_in_thresh, 2], 0.95)
        cff_avg[nset][i] = mean(cff[nset][id_nodes_in_thresh, i])
        cff_95[nset][i] = quantile(cff[nset][id_nodes_in_thresh, i], 0.95)
    end
end

# Get intersection of TopRangerSurface and BottomRangerSurface with the fault
nsets_to_intersect = ["TopRangerSurface", "BottomRangerSurface", "TopBasementSurface"]
#nsets_to_intersect = ["TopRangerSurface", "BottomRangerSurface"]
face_list_intersect = get_face_connectivity_nodesets(nsets_to_intersect, nset_names_inp, nset_inp, elements_inp[1], "tetrahedral")
nodes_coords_plot = Vector{}(undef, length(nsets_to_intersect))
for nset in eachindex(nsets_to_intersect)
    unique_nodes = unique(face_list_intersect[nset])
    all_common_nodes = Int[]
    for nset2 in eachindex(nsets_to_export)
        common_nodes = intersect(unique_nodes, face_list[nset2]) # update index of face_list according to nset for intersections
        append!(all_common_nodes, common_nodes)
    end
    all_common_nodes = unique(all_common_nodes)
    nodes_coords_plot[nset]  = vcat([nodes_inp[nodes_inp[:, 1] .== node_id, 2:4] for node_id in all_common_nodes]...)
end


# Plots
# 1. Plot average dcff on each fault vs t
using CairoMakie
CairoMakie.activate!()

t = collect(1936:2020)  # Make Auto 

ut_orange = colorant"#FF8200"
f = Figure(size=(900, 250))
ax = Axis(f[1, 1], xlabel = L"t \; \text{ [y]}", ylabel = L"\text{DCFF} \; \text{ [bar]}", xlabelsize = 16, ylabelsize = 16, xticks=1940:20:2020)
lines!(ax, t, dcff_avg[1]./1e5, color = "black", linewidth = 1.5, label=L"\text{Mean}")
lines!(ax, t, dcff_95[1]./1e5, color = ut_orange, linewidth = 1.5, label=L"\text{P95}")
xlims!(ax, 1935, 2020)
#ylims!(ax, -10, 1)
axislegend(position = :rt, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8), nbanks=1) 

ax = Axis(f[1, 2], xlabel = L"t \; \text{ [y]}", ylabel = L"\frac{\Delta\overline{\text{CFF}}}{\Delta t} \; \left[\frac{\text{bar}}{\text{y}}\right]", xlabelsize = 16, ylabelsize = 16, xticks=1940:20:2020)
lines!(ax, t[1:end-1] .+ diff(t), diff(cff_avg[1]./1e5)./diff(t), color = "black", linewidth = 1.5, label=L"\text{Mean}")
xlims!(ax, 1935, 2020)

ax = Axis(f[1, 3], xlabel = L"t \; \text{ [y]}", ylabel = L"\frac{\Delta\text{P95 CFF}}{\Delta t} \; \left[\frac{\text{bar}}{\text{y}}\right]", xlabelsize = 16, ylabelsize = 16, xticks=1940:20:2020)
lines!(ax, t[1:end-1] .+ diff(t), diff(cff_95[1]./1e5)./diff(t), color = ut_orange, linewidth = 1.5, label=L"\text{P95}")
xlims!(ax, 1935, 2020)
#ylims!(ax, -1.25, 1.25)
f
#save(join([figs_dir, "dcff_vs_t_TopRangerSurface_1936-2020.pdf"]), f, pt_per_unit = 1)

cmap = cgrad(:roma, length(nsets_to_export), categorical=true)
#cmap = cgrad(:nipy_spectral, length(nsets_to_export))
f = Figure(size=(600, 600))
ax = Axis(f[1, 1], xlabel = L"t \; \text{ [y]}", ylabel = L"\overline{\text{DCFF}} \; \text{ [bar]}", xlabelsize = 16, ylabelsize = 16, xticks=1940:20:2020)
for (i, data) in enumerate(dcff_avg)
        #norm_i = i / length(idplot)
        lbl = replace(nsets_to_export[i], "Tear_Fault_" => "TF")
        lbl = latexstring("\\text{" * lbl * "}")
        #lbl = latexstring("\\text{" * nsets_to_export[i] * "}")
        lines!(ax, t, data./1e5, color = cmap[i], linewidth = 1.5, label=lbl)
end
xlims!(ax, 1935, 2020)
#ylims!(ax, -10, 1)
axislegend(position = :lt, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8), nbanks=2) 

#ax = Axis(f[1, 2], ylabel = L"\min{\text{ DCFF}} \; \text{ [bar]}", xlabelsize = 16, ylabelsize = 16)
#ax = Axis(f[1, 2], ylabel = L"\text{P}{10} \text{ DCFF} \; \text{ [bar]}", xlabelsize = 16, ylabelsize = 16)
#for (i, data) in enumerate(dcff_10)
#    normalized_i = i / length(nsets_to_export)
#    lines!(ax, t, data./1e5, color = cmap[normalized_i], linewidth = 1.5)
#end
#xlims!(ax, 1936, 1975)
#ylims!(ax, -25, 1)

ax = Axis(f[2, 1], ylabel = L"\text{P}{95} \text{ DCFF} \; \text{ [bar]}", xlabelsize = 16, ylabelsize = 16, xticks=1940:20:2020)
#ax = Axis(f[1, 2], ylabel = L"\max \text{ DCFF} \; \text{ [bar]}", xlabelsize = 16, ylabelsize = 16)
for (i, data) in enumerate(dcff_95)
        #normalized_i = k / length(idplot)
        lines!(ax, t, data./1e5, color = cmap[i], linewidth = 1.5)
end
xlims!(ax, 1935, 2020)
#ylims!(ax, -1, 3)
f
#save(join([figs_dir, "dcff_vs_t_TearFaults_1936-2020.pdf"]), f, pt_per_unit = 1)


# 2. Plot Map
# Better like this than using a loop: unique 1 to n numbering and plot without lines showing fault segment boundaries
# we assume nodes_inp to be monotonically increasing
using GLMakie
GLMakie.activate!()

# Not rotated
fn_inp2  = "WilmingtonMesh_AllFaults_Coarse.inp" # Not rotated grid, with actual UTM coords (zone 11)
nodes_inp2, _, _, _, _ = get_connectivity_from_inp(dir_inp, fn_inp2, "tetrahedral", false)

unique_ids = sort(unique(vcat(face_list...)))   # nodes may still be repeated when >1 nodesets share nodes
coords_map = Vector{}(undef, 1)
#coords_map[1] = nodes_inp[in.(nodes_inp[:, 1], Ref(unique_ids)), 2:end] # for node_coords_plot
coords_map[1] = nodes_inp2[in.(nodes_inp2[:, 1], Ref(unique_ids)), 2:end] # for coastline

# vals nodes
ids_nset_export = [2, 9, 10]
tstep = "step_89_nodes.h5"
#tstep = "step_89_nodes.dat" # t=1958
#tstep = "step_97_nodes.dat" # t=1960
#tstep = "step_149_nodes.dat" # t=1973
#tstep = "step_157_nodes.dat" # t=1975
#tstep = "step_217_nodes.dat" # t=1990
#tstep = "step_277_nodes.dat" # t=2005
#tstep = "step_337_nodes.dat" # t=2020
abaqus_fn_nod = readdir(join([sim_dir, nsets_to_export[1]]))
abaqus_fn_nod = sort(abaqus_fn_nod, by=sort_steps)
idx = findfirst(==(tstep), abaqus_fn_nod)
dcff_plot = Vector{Matrix{Union{Int64,Float64}}}(undef, length(nsets_to_export))
for nset in eachindex(nsets_to_export)
    dcff_plot[nset] = dcff[nset][idx-1]
end
tmp = sortslices(vcat(dcff_plot...), dims=1, by=x->x[1])
values_plot = zeros(Float64, length(unique_ids), 3)
for (i, node) in enumerate(unique_ids) 
    id_rows = tmp[:,1] .== node
    values_plot[i, :] = mean(tmp[id_rows, 2:end], dims=1)
end
vals_nodes = Vector{}(undef, 1)
col_idx = 3;
vals_nodes[1] = values_plot[:, col_idx] ./ 1e5         # bar, column index 1 = dcff, 2 = dsn, 3 = dtau
if col_idx == 1
    labl = L"\text{DCFF [bar]}"
elseif col_idx == 2
    labl = L"\Delta \sigma_\text{n}^{\text{eff}} \; \text{ [bar]}"
else
    labl = L"\Delta |\tau| \; \text{ [bar]}"
end

# face connectivity
face_list_map = Vector{}(undef, 1)
face_list_map[1], _ = renumber_array_1_to_n(vcat(face_list...))

#mycmap = cgrad(:oleron)
#mycmap = cgrad(:lajolla)
mycmap = cgrad(:balance)
#-0.8*pi, 0.2*pi
fig_descriptor = (cmap=mycmap, cbar_label=labl, fig_size = (1200, 900), plot_height = 800, cbar_height = 300, plot_azi = -0.5*pi, plot_elev = 0.5*pi, plot_linew=0.2, xlms=(minimum(coords_map[1][:,1]), maximum(coords_map[1][:,1])), ylms=(minimum(coords_map[1][:,2]), maximum(coords_map[1][:,2])), save_fig = true, save_dir = figs_dir, title = "TopRangerSurface_dtau_plastic-120_tstep89", target_dpi=400, target_size_in=(12,9))
#plotNodalValues(face_list_map, coords_map, vals_nodes, fig_descriptor, node_coords_plot = nodes_coords_plot)
plotNodalValues(face_list_map, coords_map, vals_nodes, fig_descriptor)

coast_data = readdlm(join([data_dir, "cal_outline.dat"]))
cst = [coast_data[coast_data[:, 3] .== i, :] for i in 1:18] # third column is "category"
coastline = [cst[6], cst[17], cst[18]]
#coastline = [hcat(coastline[:, 2], coastline[:, 1], coastline[:, 3])]
# Plotting
plotNodalValues(face_list_map, coords_map, vals_nodes, fig_descriptor, line_data=coastline)


# test tractions
#=
S = [15 0 0; 0 12 0; 0 0 10]
coord = [2 2 1; 6 3 0; 4 1 1]
f = coord[2,:]-coord[1,:]
g = coord[3,:]-coord[1,:]
n = cross(f, g) 
for i in 3:-1:1 # start with z component, we want the outward normal vector
    if n[i] != 0
        sgn = n[i] >= 0 ? 1 : -1
        n = n .* sgn / norm(n)
        break
    end
end

t = S*n
tn = dot(n, t)
ts = sqrt(norm(t)^2 - tn^2)
=#
