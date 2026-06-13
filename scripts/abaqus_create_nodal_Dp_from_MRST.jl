#=
Couple MRST → Abaqus: compute nodal pore-pressure changes (Δp) for the Abaqus
geomechanics model from cell-based pressures in the MRST flow simulation.
For each timestep, MRST cell pressures are averaged onto Abaqus mesh nodes (via the
node→cell map), differenced from the first step, and written per-timestep as the
nodal_dp_<year>.dat forcing files that drive the geomechanics.
Optional diagnostic map (GLMakie) is off by default; the script runs headless.
=#
include(joinpath(@__DIR__, "..", "config.jl"))

using MAT, DelimitedFiles, StatsBase, Printf
include(tool(:abaqus));    using .Abaqus          # get_connectivity_from_inp
include(tool(:mrst));      using .Mrst            # map_nodes_to_cells
include(tool(:makieplot));      using .MakiePolyPlot            # plot

const SAVE_TO_DEPOSIT = true      # true → (re)write nodal_dp_*.dat into the deposit
const SHOW_PLOT       = false      # true → 3D Δp diagnostic map (GLMakie, needs display)

## Filenames
fn_inp = "WilmingtonMesh_AllFaults_Rotated_Coarse.inp"
fn_G   = "wilmington_ranger_G_withGeometry.mat"
fn_p   = "states_p_Pa.mat"
fn_t   = "states_t_s.mat"

## Output: per-timestep nodal Δp forcing for Abaqus
dp_out = SAVE_TO_DEPOSIT ? FOR_ABAQUS_PRESSURE : joinpath(OUT_DIR, "from_mrst", "dp_for_abaqus")
isdir(dp_out) || mkpath(dp_out)

## Abaqus mesh node coordinates
nodes_inp, _, _, _, _ = get_connectivity_from_inp(GRIDS_ABAQUS, fn_inp, "tetrahedral")

## MRST grid + node→cell map
G    = matread(joinpath(GRIDS_MRST, fn_G))["G"]
f2cn = Int.(G["cells"]["f2n"])
n2fn = Int.(G["faces"]["n2fn"])
nodes_to_cells, faces_nodes = map_nodes_to_cells(G, n2fn, f2cn)

## Time axis (years) and cell pressures; Δp relative to the first step (pre-injection)
dty = 0.25                                                          # MRST timestep, years
t = 1936 - dty .+ matread(joinpath(MRST_RESULT_DATA, fn_t))["t"] ./ 31556952
p = matread(joinpath(MRST_RESULT_DATA, fn_p))["p"]
dp = p[:, 1:end] .- p[:, 1]                                         # Δp from first step (no injection yet)

## Average cell Δp onto MRST nodes, per timestep
nodal_dp = Vector{}(undef, size(dp, 2))
for i in eachindex(nodal_dp)
    nodal_dp[i] = Dict((node_id, mean(dp[cell_ids, i])) for (node_id, cell_ids) in nodes_to_cells)
end

## Map MRST nodes → global (.inp) node ids
# The MRST grid is unrotated and the Abaqus grid rotated, but global node ids match
# between them, so the map is a direct id↔id correspondence (see wilmington-prep / mrst_create_inputs.jl).
map_inp_to_mrst_node = readdlm(joinpath(GRIDS_MRST, "map_inp_to_mrst_node.dat"), '\t', Int64, '\n')
nodal_dp_inp = Matrix{Union{Int64,Float64}}(undef, Int(G["nodes"]["num"]), size(dp, 2) + 1)
nodal_dp_inp[:, 1] = map_inp_to_mrst_node[:, 1]
for i in eachindex(nodal_dp)
    nodal_dp_inp[:, i+1] = [nodal_dp[i][row[2]] for row in eachrow(map_inp_to_mrst_node)]
end
nodal_dp_inp_sorted = sortslices(nodal_dp_inp, dims=1)              # by global node id

## Internal verification check
# scratch = joinpath(OUT_DIR, "from_mrst", "nodal_dp_1936_to_2020_tri")
#for yr in ["1961.0", "1985.25", "2010.5"]    # random years that exist on disk
#    fn = "nodal_dp_$(yr).dat"
#    a = readdlm(joinpath(scratch, fn), ',', Float64)
#    b = readdlm(joinpath(FOR_ABAQUS_PRESSURE, fn), ',', Float64)
#    println(rpad(fn, 22), " shape ", size(a)==size(b), "  maxdiff ", size(a)==size(b) ? maximum(abs.(a .- b)) : NaN)
#end

## Write per-timestep nodal Δp (only nodes with |Δp| > 10 Pa; fixed-point to be Abaqus-safe)
for i in eachindex(nodal_dp)
    open(joinpath(dp_out, "nodal_dp_$(round(t[i], digits=2)).dat"), "w") do f
        for row in eachrow(hcat(nodal_dp_inp_sorted[:, 1], nodal_dp_inp_sorted[:, i+1]))
            if abs(row[2]) > 10
                write(f, string(Int(row[1])) * "," * @sprintf("%.2f", row[2]) * "\n")
            end
        end
    end
end

## Optional diagnostic: 3D map of nodal Δp at one timestep (GLMakie, needs display), to verify against MRST maps.
if SHOW_PLOT
    using GLMakie, Colors, ColorSchemes, LaTeXStrings
    GLMakie.activate!()
    t_id = 157 # 89=1958, 157=1975, 337=2020
    nodal_dp_inp_plot = Matrix{Union{Int64,Float64}}(undef, size(nodal_dp_inp)...)
    nodal_dp_inp_plot[:, 1] = map_inp_to_mrst_node[:, 2]
    nodal_dp_inp_plot[:, 2:end] = nodal_dp_inp[:, 2:end]
    nodal_dp_inp_plot = sortslices(nodal_dp_inp_plot, dims=1)

    coords_map = [G["nodes"]["coords"]]
    cells_nodes_map_plot = [faces_nodes]
    cells_values_nodes = [nodal_dp_inp_plot[:, t_id] ./ 1e5]        # bar
    fig_descriptor = (cmap=cgrad(:curl), cbar_label=L"$\Delta p$ [bar]", fig_size=(800,600),
        plot_height=900, cbar_height=200, plot_azi=1.57π, plot_elev=0.47π, plot_linew=0.1,
        xlms=(378449, 400339), ylms=(3.728e6, 3.752e6),
        save_fig=false, save_dir=FIGS_DIR, title="dp_reservoir_$(round(Int, t[t_id]))",
        target_dpi=350, target_size_in=(8,6))
    plotNodalValues(cells_nodes_map_plot, coords_map, cells_values_nodes, fig_descriptor)
end