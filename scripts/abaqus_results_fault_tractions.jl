#=
Fault tractions on the Wilmington model (plastic-120). Reproduces:
  Fig 14  — Wilmington fault: Δ|τ|, Δσ_n^eff, DCFF maps at 1958/1975/1990/2020 (rel. 1936)
  Fig 15b — mean & P95 DCFF time series on tear faults (no area filter)
  Fig 15c — DCFF map on TF4 + TF9 at 1960
Map sections are 3D (GLMakie) and REQUIRE A DISPLAY; the time series (CairoMakie) is headless-OK.
=#

include(joinpath(@__DIR__, "..", "config.jl"))

using Colors, ColorSchemes, LaTeXStrings
using DelimitedFiles, Statistics
include(tool(:faultops));  using .FaultOperations   # face_normals_and_centroids, tractions_from_sij_and_n
include(tool(:abaqus));    using .Abaqus            # get_connectivity_from_inp, get_face_connectivity_*
include(tool(:arrayops));  using .ArrayOps          # reduce_matrix_with_repeated_ids, renumber_array_1_to_n
include(tool(:fileio));    using .Fileio            # sort_steps  (file fileIO.jl — case fixed via tool())
include(tool(:makieplot)); using .MakiePolyPlot     # plotNodalValues

## ---- Configuration -------------------------------------------------------------------
sim_name = "plastic-120"
mu       = 0.6

const YEAR_TO_STEP = Dict(1958=>89, 1960=>97, 1975=>157, 1990=>217, 2020=>337)  # extracted-file step index per year

const CASES = Dict(
    "fig14" => (nsets      = ["Wilmington"],
                map_nsets  = ["Wilmington"],
                map_years  = [1958, 1975, 1990, 2020],
                map_cols   = [3, 2, 1],           # 3=Δ|τ|, 2=Δσ_n^eff, 1=DCFF
                timeseries = false,
                azi = 0.3π, elev = 0.05π,
                clim = (-4, 4),
                z_aspect = 1.0, lift = 50.0),            
    "fig15" => (nsets      = ["Tear_Fault_4","Tear_Fault_9","Tear_Fault_10a","Tear_Fault_10b","Tear_Fault_11"],
                map_nsets  = ["Tear_Fault_4","Tear_Fault_9"],
                map_years  = [1960],
                map_cols   = [1],                # DCFF
                timeseries = true,
                azi = -0.3π, elev = 0.15π,
                clim = nothing,
                z_aspect = 0.5, lift = 50.0),
)
case = CASES["fig14"]      # ← select figure here

ext_S     = joinpath(abaqus_extracted(sim_name), "S")        # source stress (recompute path)
FAULT_DIR = datapath("for_fault_analysis", sim_name)
tractions_read_dir  = joinpath(FAULT_DIR, "tractions")               # deposited (read, default)
tractions_write_dir = joinpath(OUT_DIR, "fault_quantities", sim_name) # recompute cache (write)
compute_tractions   = false   # false → read deposited (fast); true → recompute from ext_S (slow)

const SAVE_TO_DEPOSIT = true
figdir = SAVE_TO_DEPOSIT ? datapath("figures", "abaqus", sim_name) : FIGS_DIR
isdir(figdir) || mkpath(figdir)

## ---- Geometry (timestep-independent) -------------------------------------------------
fn_inp = "WilmingtonMesh_AllFaults_Rotated_ENTIRE_DOMAIN_Coarse.inp"
nodes_inp, elements_inp, elset_names_inp, nset_inp, nset_names_inp = get_connectivity_from_inp(GRIDS_ABAQUS, fn_inp, "tetrahedral")
face_list = get_face_connectivity_nodesets_C3D4(case.nsets, nset_names_inp, nset_inp, elements_inp[1])
face_unit_normals, face_centroids = face_normals_and_centroids(face_list, nodes_inp)

# Non-rotated mesh for display coords (UTM zone 11) 
fn_inp2 = "WilmingtonMesh_AllFaults_Coarse.inp"
nodes_inp2, _, _, _, _ = get_connectivity_from_inp(GRIDS_ABAQUS, fn_inp2, "tetrahedral", false)

## ---- Tractions → CFF / DCFF per nodeset ----------------------------------------------
nodal_tractions = [Vector{Any}() for _ in eachindex(case.nsets)]
dcff = [Vector{Any}() for _ in eachindex(case.nsets)]
cff  = Vector{Matrix{Float64}}(undef, length(case.nsets))
cff0 = [Vector{Any}() for _ in eachindex(case.nsets)]
for nset in eachindex(case.nsets)
    src_dir = joinpath(ext_S, case.nsets[nset])
    fns = sort(readdir(src_dir), by=sort_steps)
    cff_cols = Vector{Vector{Float64}}()
    for (i, fn) in enumerate(fns)
        stem = splitext(fn)[1]
        if compute_tractions
            splitext(fn)[2] == ".dat" || error("Expected .dat stress file, got $fn")
            nodal_stress = readdlm(joinpath(src_dir, fn), ',', Float64, '\n', skipstart=1)   # NODE_LABEL, Sxx, Syy, Szz, Sxy, Sxz, Syz
            wdir = joinpath(tractions_write_dir, case.nsets[nset]); isdir(wdir) || mkpath(wdir)
            push!(nodal_tractions[nset],
                  tractions_from_sij_and_n(nodal_stress, nodes_inp, face_unit_normals[nset], face_centroids[nset], wdir, "tractions_$(stem).dat"))
        else
            # deposited tractions are always tab-delimited text — read uniformly
            data = readdlm(joinpath(tractions_read_dir, case.nsets[nset], "tractions_$(stem).dat"), '\t', Float64, '\n')
            push!(nodal_tractions[nset], data)
        end

        i == 1 && check_t0_tractions(nodal_tractions[nset][1], case.nsets[nset])   # faults reach the basement: μ = 0.6

        if i == 1
            cff0[nset] = nodal_tractions[nset][1][:, end] .+ mu*nodal_tractions[nset][1][:, 2]
        else
            dtau  = nodal_tractions[nset][i][:, end] .- nodal_tractions[nset][1][:, end]
            dsn   = nodal_tractions[nset][i][:, 2]   .- nodal_tractions[nset][1][:, 2]
            cff_i = nodal_tractions[nset][i][:, end] .+ mu*nodal_tractions[nset][i][:, 2]
            push!(cff_cols, cff_i)
            push!(dcff[nset], hcat(nodal_tractions[nset][1][:, 1], cff_i .- cff0[nset], dsn, dtau))
        end
    end
    cff[nset] = reduce(hcat, cff_cols)
    println("Tractions for $(case.nsets[nset]) completed.")
end

## ---- Reference horizons for map overlays (Producing interval, Top basement) ----------
nsets_to_intersect = ["TopRangerSurface", "BottomRangerSurface", "TopBasementSurface"]
face_list_intersect = get_face_connectivity_nodesets(nsets_to_intersect, nset_names_inp, nset_inp, elements_inp[1], "tetrahedral")
nodes_coords_plot = Vector{}(undef, length(nsets_to_intersect))
for ns in eachindex(nsets_to_intersect)
    unique_nodes = unique(face_list_intersect[ns])
    common = Int[]
    for n2 in eachindex(case.map_nsets)
        idx = findfirst(==(case.map_nsets[n2]), case.nsets)
        append!(common, intersect(unique_nodes, face_list[idx]))
    end
    common = unique(common)
    nodes_coords_plot[ns] = vcat([nodes_inp2[nodes_inp2[:, 1] .== id, 2:4] for id in common]...)
end

## ---- Fig 15b: time series (CairoMakie) ----------------------------------
if case.timeseries
    using CairoMakie
    CairoMakie.activate!()
    t = collect(1936:2020)
    n_steps = length(dcff[1])
    @assert n_steps == length(t) "DCFF steps ($n_steps) ≠ time axis ($(length(t))) — check curated extracted set"

    dcff_avg = [[mean(dcff[ns][i][:, 2]) for i in 1:n_steps] for ns in eachindex(case.nsets)]
    dcff_95  = [[quantile(dcff[ns][i][:, 2], 0.95) for i in 1:n_steps] for ns in eachindex(case.nsets)]

    cmap = cgrad(:roma, length(case.nsets), categorical=true)
    f = Figure(size=(600, 600))
    ax = Axis(f[1, 1], xlabel=L"t \;[\text{y}]", ylabel=L"\overline{\text{DCFF}}\;[\text{bar}]", xticks=1940:20:2020)
    for (i, d) in enumerate(dcff_avg)
        lbl = latexstring("\\text{" * replace(case.nsets[i], "Tear_Fault_"=>"TF") * "}")
        lines!(ax, t, d ./ 1e5, color=cmap[i], linewidth=1.5, label=lbl)
    end
    xlims!(ax, 1935, 2020); axislegend(position=:lt, framevisible=true, nbanks=2)
    ax = Axis(f[2, 1], xlabel=L"t \;[\text{y}]", ylabel=L"\text{P95 DCFF}\;[\text{bar}]", xticks=1940:20:2020)
    for (i, d) in enumerate(dcff_95)
        lines!(ax, t, d ./ 1e5, color=cmap[i], linewidth=1.5)
    end
    xlims!(ax, 1935, 2020)
    SAVE_TO_DEPOSIT && save(joinpath(figdir, "dcff_vs_t_tearfaults_$(sim_name).pdf"), f, pt_per_unit=1)
    display(f)
end

## ---- Fig 14 / 15c: maps (GLMakie, needs display) -------------------------------------
using GLMakie
GLMakie.activate!()

# Combine the mapped nodesets onto one surface
map_idx    = [findfirst(==(n), case.nsets) for n in case.map_nsets]
face_comb  = vcat([face_list[k] for k in map_idx]...)
unique_ids = sort(unique(face_comb))
coords_map = Vector{}(undef, 1)
coords_map[1] = nodes_inp2[in.(nodes_inp2[:, 1], Ref(unique_ids)), 2:end]
face_list_map = Vector{}(undef, 1)
face_list_map[1], _ = renumber_array_1_to_n(face_comb)

col_labels = Dict(1 => L"\text{DCFF [bar]}", 2 => L"\Delta\sigma_\text{n}^{\text{eff}}\;[\text{bar}]", 3 => L"\Delta|\tau|\;[\text{bar}]")
col_tags   = Dict(1 => "dcff", 2 => "dsn", 3 => "dtau")

for yr in case.map_years
    step = YEAR_TO_STEP[yr]
    fns  = sort(readdir(joinpath(ext_S, case.map_nsets[1])), by=sort_steps)
    idx  = findfirst(f -> occursin("step_$(step)_", f), fns)
    isnothing(idx) && error("year $yr (step_$(step)) not found in $(case.map_nsets[1])")

    tmp = sortslices(vcat([dcff[k][idx-1] for k in map_idx]...), dims=1, by=x->x[1])
    values_plot = zeros(Float64, length(unique_ids), 3)
    for (i, node) in enumerate(unique_ids)
        rows = tmp[:, 1] .== node
        values_plot[i, :] = mean(tmp[rows, 2:end], dims=1)
    end

    for col in case.map_cols
        vals_nodes = Vector{}(undef, 1)
        vals_nodes[1] = values_plot[:, col] ./ 1e5      # bar
        ttl = "$(join(case.map_nsets, "_"))_$(col_tags[col])_$(sim_name)_$(yr)"
        fig_descriptor = (cmap=cgrad(:balance), cbar_label=col_labels[col], clim=case.clim,
            fig_size=(1200,900), plot_height=800, cbar_height=300,
            plot_azi = case.azi, plot_elev = case.elev, plot_linew=0.2,
            z_aspect=case.z_aspect, lift=case.lift,
            xlms=(minimum(coords_map[1][:,1]), maximum(coords_map[1][:,1])),
            ylms=(minimum(coords_map[1][:,2]), maximum(coords_map[1][:,2])),
            save_fig=SAVE_TO_DEPOSIT, save_dir=figdir, title=ttl, target_dpi=400, target_size_in=(12,9))
        plotNodalValues(face_list_map, coords_map, vals_nodes, fig_descriptor, node_coords_plot = nodes_coords_plot, line_data = nothing)
    end
end