#=
Surface tractions on the Wilmington model (plastic-120). Reproduces Fig 17 b–e:
  17b — Top Pico:   Δ|τ|, DCFF maps at 1958
  17c — Top Ranger: Δσ_n^eff, Δ|τ|, DCFF maps at 1958
  17d — Top Pico:   mean & P95 DCFF + yearly rates (filtered to surface ‖u_h‖ ≥ 1 m)
  17e — Top Ranger: same as 17d
Map sections are 3D (GLMakie) and REQUIRE A DISPLAY; the time series (CairoMakie) is headless-OK.
=#
include(joinpath(@__DIR__, "..", "config.jl"))

using Colors, ColorSchemes, LaTeXStrings
using DelimitedFiles, Statistics
include(tool(:faultops));  using .FaultOperations   # face_normals_and_centroids, tractions_from_sij_and_n
include(tool(:abaqus));    using .Abaqus            # get_connectivity_from_inp, get_face_connectivity_*
include(tool(:arrayops));  using .ArrayOps          # reduce_matrix_with_repeated_ids, renumber_array_1_to_n
include(tool(:fileio));    using .Fileio            # sort_steps
include(tool(:makieplot)); using .MakiePolyPlot     # plotNodalValues

## ---- Configuration -------------------------------------------------------------------
sim_name = "plastic-120"
mu       = 0.6
const YEAR_TO_STEP = Dict(1958=>89, 1960=>97, 1975=>157, 1990=>217, 2020=>337)

const CASES = Dict(
    "pico"   => (nset="PicoSurface",      surface="Top Pico",
                 map_year=1958, map_cols=[3,1],     # Δ|τ|, DCFF        (17b)
                 azi=-0.5π, elev=0.5π, z_aspect=0.1, lift=0.0,
                 clims=Dict()),
    "ranger" => (nset="TopRangerSurface", surface="Top Ranger",
                 map_year=1958, map_cols=[2,3,1],   # Δσ_n, Δ|τ|, DCFF  (17c)
                 azi=-0.5π, elev=0.5π, z_aspect=0.1, lift=0.0,
                 clims=Dict(1 => (-35, 35), 2 => (-55, 55), 3 => (-20, 20))),
)
case = CASES["pico"]      # ← select surface here

ext_S     = joinpath(abaqus_extracted(sim_name), "S")
FAULT_DIR = datapath("for_fault_analysis", sim_name)
tractions_read_dir  = joinpath(FAULT_DIR, "tractions")
tractions_write_dir = joinpath(OUT_DIR, "fault_quantities", sim_name)
compute_tractions   = false   # false → read deposited (fast); true → recompute from ext_S (slow)

const SAVE_TO_DEPOSIT = true
figdir = SAVE_TO_DEPOSIT ? datapath("figures", "abaqus", sim_name) : FIGS_DIR
isdir(figdir) || mkpath(figdir)

## ---- Geometry (timestep-independent) -------------------------------------------------
fn_inp = "WilmingtonMesh_AllFaults_Rotated_ENTIRE_DOMAIN_Coarse.inp"
nodes_inp, elements_inp, elset_names_inp, nset_inp, nset_names_inp = get_connectivity_from_inp(GRIDS_ABAQUS, fn_inp, "tetrahedral")
face_list = get_face_connectivity_nodesets_C3D4([case.nset], nset_names_inp, nset_inp, elements_inp[1])
face_unit_normals, face_centroids = face_normals_and_centroids(face_list, nodes_inp)

fn_inp2 = "WilmingtonMesh_AllFaults_Coarse.inp"      # non-rotated, UTM zone 11 (display coords)
nodes_inp2, _, _, _, _ = get_connectivity_from_inp(GRIDS_ABAQUS, fn_inp2, "tetrahedral", false)

# Display coords for this surface (built once, ABOVE the filter that needs it)
unique_ids = sort(unique(face_list[1]))
coords_map = Vector{}(undef, 1)
coords_map[1] = nodes_inp2[in.(nodes_inp2[:, 1], Ref(unique_ids)), 2:end]
face_list_map = Vector{}(undef, 1)
face_list_map[1], _ = renumber_array_1_to_n(face_list[1])

## ---- Area filter: surface ‖u_h‖ ≥ 1 m bounding box (17d/e) ---------------------------
uh = readdlm(joinpath(FAULT_DIR, "nodal_coords_and_u_h_1958.dat"))
above = uh[:, end] .>= 1.0
xmin, xmax = extrema(uh[above, 1]); ymin, ymax = extrema(uh[above, 2])
id_in_thresh = vec(all(hcat(coords_map[1][:,1] .>= xmin, coords_map[1][:,1] .<= xmax,
                            coords_map[1][:,2] .>= ymin, coords_map[1][:,2] .<= ymax), dims=2))

## ---- Tractions → CFF / DCFF ----------------------------------------------------------
src_dir = joinpath(ext_S, case.nset)
fns = sort(readdir(src_dir), by=sort_steps)
nodal_tractions = Vector{Any}()
dcff = Vector{Any}()
cff_cols = Vector{Vector{Float64}}()
cff0 = Float64[]
for (i, fn) in enumerate(fns)
    stem = splitext(fn)[1]
    if compute_tractions
        splitext(fn)[2] == ".dat" || error("Expected .dat stress file, got $fn")
        nodal_stress = readdlm(joinpath(src_dir, fn), ',', Float64, '\n', skipstart=1)   # NODE_LABEL, Sxx, Syy, Szz, Sxy, Sxz, Syz
        wdir = joinpath(tractions_write_dir, case.nset); isdir(wdir) || mkpath(wdir)
        push!(nodal_tractions, tractions_from_sij_and_n(nodal_stress, nodes_inp, face_unit_normals[1], face_centroids[1], wdir, "tractions_$(stem).dat"))
    else
        push!(nodal_tractions, readdlm(joinpath(tractions_read_dir, case.nset, "tractions_$(stem).dat"), '\t', Float64, '\n'))
    end

    i == 1 && check_t0_tractions(nodal_tractions[1], case.nset)   

    if i == 1
        global cff0 = nodal_tractions[1][:, end] .+ mu*nodal_tractions[1][:, 2]
    else
        dtau  = nodal_tractions[i][:, end] .- nodal_tractions[1][:, end]
        dsn   = nodal_tractions[i][:, 2]   .- nodal_tractions[1][:, 2]
        cff_i = nodal_tractions[i][:, end] .+ mu*nodal_tractions[i][:, 2]
        push!(cff_cols, cff_i)
        push!(dcff, hcat(nodal_tractions[1][:, 1], cff_i .- cff0, dsn, dtau))
    end
end
cff = reduce(hcat, cff_cols)
println("Tractions for $(case.nset) completed.")


## ---- Fig 17d/e: time series with rate panels (CairoMakie) ----------------------------
using CairoMakie
CairoMakie.activate!()
t = collect(1936:2020)
n_steps = length(dcff)
@assert n_steps == length(t) "DCFF steps ($n_steps) ≠ time axis ($(length(t)))"

dcff_avg = [mean(dcff[i][id_in_thresh, 2]) for i in 1:n_steps]
dcff_95  = [quantile(dcff[i][id_in_thresh, 2], 0.95) for i in 1:n_steps]
cff_avg  = [mean(cff[id_in_thresh, i]) for i in 1:n_steps]
cff_95   = [quantile(cff[id_in_thresh, i], 0.95) for i in 1:n_steps]
ut_orange = colorant"#FF8200"
tmid = t[1:end-1] .+ diff(t)/2

f = Figure(size=(900, 250))
ax = Axis(f[1,1], xlabel=L"t \;[\text{y}]", ylabel=L"\text{DCFF}\;[\text{bar}]", xticks=1940:20:2020)
lines!(ax, t, dcff_avg ./ 1e5, color=:black, linewidth=1.5, label=L"\text{Mean}")
lines!(ax, t, dcff_95  ./ 1e5, color=ut_orange, linewidth=1.5, label=L"\text{P95}")
xlims!(ax, 1935, 2020); axislegend(position=:lt, framevisible=true, backgroundcolor=RGBAf(1,1,1,0.8))
ax = Axis(f[1,2], xlabel=L"t \;[\text{y}]", ylabel=L"\frac{\Delta\overline{\text{CFF}}}{\Delta t}\;\left[\frac{\text{bar}}{\text{y}}\right]", xticks=1940:20:2020)
lines!(ax, tmid, diff(cff_avg ./ 1e5) ./ diff(t), color=:black, linewidth=1.5)
xlims!(ax, 1935, 2020)
ax = Axis(f[1,3], xlabel=L"t \;[\text{y}]", ylabel=L"\frac{\Delta\text{P95 CFF}}{\Delta t}\;\left[\frac{\text{bar}}{\text{y}}\right]", xticks=1940:20:2020)
lines!(ax, tmid, diff(cff_95 ./ 1e5) ./ diff(t), color=ut_orange, linewidth=1.5)
xlims!(ax, 1935, 2020)
SAVE_TO_DEPOSIT && save(joinpath(figdir, "dcff_vs_t_$(case.nset)_$(sim_name).pdf"), f, pt_per_unit=1)
display(f)

## ---- Fig 17b/c: maps (GLMakie, needs display) ----------------------------------------
using GLMakie
GLMakie.activate!()
col_labels = Dict(1 => L"\text{DCFF [bar]}", 2 => L"\Delta\sigma_\text{n}^{\text{eff}}\;[\text{bar}]", 3 => L"\Delta|\tau|\;[\text{bar}]")
col_tags   = Dict(1 => "dcff", 2 => "dsn", 3 => "dtau")

step = YEAR_TO_STEP[case.map_year]
idx  = findfirst(f -> occursin("step_$(step)_", f), fns)
isnothing(idx) && error("year $(case.map_year) (step_$(step)) not found in $(case.nset)")

tmp = sortslices(dcff[idx-1], dims=1, by=x->x[1])
values_plot = zeros(Float64, length(unique_ids), 3)
for (i, node) in enumerate(unique_ids)
    rows = tmp[:, 1] .== node
    values_plot[i, :] = mean(tmp[rows, 2:end], dims=1)
end

# Load coastline & plot
coast_data = readdlm(joinpath(GEO_DIR, "cal_outline.dat"))
cst = [coast_data[coast_data[:, 3] .== i, :] for i in 1:18]   # col 3 = category
coastline = [cst[6], cst[17], cst[18]]

for col in case.map_cols
    vals_nodes = Vector{}(undef, 1)
    vals_nodes[1] = values_plot[:, col] ./ 1e5
    ttl = "$(case.nset)_$(col_tags[col])_$(sim_name)_$(case.map_year)"
    xl, yl = extrema(coords_map[1][:, 1]), extrema(coords_map[1][:, 2])
    fig_descriptor = (cmap=cgrad(:balance), cbar_label=col_labels[col], clim=get(case.clims, col, nothing),
        fig_size=(1200,900), plot_height=800, cbar_height=300,
        plot_azi=case.azi, plot_elev=case.elev, z_aspect=case.z_aspect, lift=case.lift, plot_linew=0.2,
        equal_xy=true, fontsize_annot=40, coast_linewidth=2,
        scalebar_km = col == first(case.map_cols) ? 5 : nothing,
        scalebar_xy = (xl[1] + 0.33 * (xl[2] - xl[1]), yl[1] + 0.05 * (yl[2] - yl[1])),
        north_arrow = col == first(case.map_cols),
        xlms=(minimum(coords_map[1][:,1]), maximum(coords_map[1][:,1])),
        ylms=(minimum(coords_map[1][:,2]), maximum(coords_map[1][:,2])),
        save_fig=SAVE_TO_DEPOSIT, save_dir=figdir, title=ttl, target_dpi=400, target_size_in=(12,9))
    plotNodalValues(face_list_map, coords_map, vals_nodes, fig_descriptor, line_data=coastline)
end