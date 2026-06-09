#=
Maps of nodal displacement on the Wilmington model surface (OverburdenPico elset),
overlaid with the coastline. Regenerates the canonical paper panels for one variant:
  plastic-130 → Fig 12c (u_z, 1965 & 2020)
  plastic-120 → Fig 13c (u_z, 1965 & 2020) + Fig 17a (u_h, 1958)
  
Note: These are 3D plots — REQUIRE A DISPLAY (GLMakie).
=#
include(joinpath(@__DIR__, "..", "config.jl"))

using GLMakie
GLMakie.activate!()
using Colors, ColorSchemes, LaTeXStrings
using DelimitedFiles
include(tool(:abaqus));    using .Abaqus         # get_connectivity_from_inp, get_face_connectivity_*
include(tool(:arrayops));  using .ArrayOps       # renumber_array_1_to_n
include(tool(:makieplot)); using .MakiePolyPlot  # plotNodalValues  (file makiePolyPlot.jl)

## Canonical deposited panels
const PANELS = Dict(
    "plastic-130" => [(tstep="step_117", comp="vertical", year=1965, clim=(-8,8),  azi=-0.5π, elev=0.25π),
                      (tstep="step_337", comp="vertical", year=2020, clim=(-8,8),  azi=-0.5π, elev=0.25π)],
    "plastic-120" => [(tstep="step_117", comp="vertical", year=1965, clim=(-9,9),  azi=-0.5π, elev=0.25π),
                      (tstep="step_337", comp="vertical", year=2020, clim=(-9,9),  azi=-0.5π, elev=0.25π),
                      (tstep="step_89",  comp="horizontal", year=1958, clim=nothing, azi=-0.5π, elev=0.5π)],
)

## Select variant here; everything below follows
sim_name = "plastic-120"
panel_idx = 3 # see PANELS above 
p = PANELS[sim_name][panel_idx]
ext_U    = joinpath(abaqus_extracted(sim_name), "U")     # extracted nodal-displacement sets

const SAVE_TO_DEPOSIT = true
figdir = SAVE_TO_DEPOSIT ? datapath("figures", "abaqus", sim_name) : FIGS_DIR
isdir(figdir) || mkpath(figdir)

## ---- Geometry (timestep-independent: computed once) ----------------------------------
# Nodeset faces (inactive unless nsets_to_export is non-empty)
fn_inp_ns = "WilmingtonMesh_AllFaults_Rotated_ENTIRE_DOMAIN_Coarse.inp"   # simplified elsets, for plotting
nsets_to_export = String[]
if !isempty(nsets_to_export)
    _, elements_inp_ns, _, nset_inp, nset_names_inp = get_connectivity_from_inp(GRIDS_ABAQUS, fn_inp_ns, "tetrahedral", false)
    face_list_ns = get_face_connectivity_nodesets(nsets_to_export, nset_names_inp, nset_inp, elements_inp_ns[1], "tetrahedral")
else
    face_list_ns = []
end

# Element-set faces (active path: OverburdenPico)
fn_inp_el = "WilmingtonMesh_AllFaults_Rotated_Coarse.inp"
elsets_to_export = ["OverburdenPico"]
_, elements_inp, elset_names_inp, _, _ = get_connectivity_from_inp(GRIDS_ABAQUS, fn_inp_el, "tetrahedral")
face_list_el = get_face_connectivity_elsets(elsets_to_export, elset_names_inp, elements_inp, "tetrahedral")

# Combined sets (names + faces), built once — no push! in a loop
set_names = vcat(nsets_to_export, elsets_to_export)
face_list = vcat(face_list_ns, face_list_el)

# Non-rotated mesh for display coordinates (UTM zone 11)
fn_inp2 = "WilmingtonMesh_AllFaults_Coarse.inp"
nodes_inp2, _, _, _, _ = get_connectivity_from_inp(GRIDS_ABAQUS, fn_inp2, "tetrahedral", false)

# Per-set coords + renumbered faces (assumes nodes_inp2 is sorted by node id, as in the original)
nsets         = length(set_names)
coords_map    = Vector{}(undef, nsets)
face_list_map = Vector{}(undef, nsets)
for set in 1:nsets
    unique_ids = sort(unique(face_list[set]))
    coords_map[set] = nodes_inp2[in.(nodes_inp2[:, 1], Ref(unique_ids)), 2:end]
    face_list_map[set], _ = renumber_array_1_to_n(face_list[set])
end

# Coastline (UTM zone 11); column 3 is a category id
coast_data = readdlm(joinpath(GEO_DIR, "cal_outline.dat"))
cst = [coast_data[coast_data[:, 3] .== i, :] for i in 1:18]
coastline = [cst[6], cst[17], cst[18]]

## ---- Selected panel: read displacements, build values, plot, save --------------------
tstep, comp = p.tstep, p.comp

vals_nodes = Vector{}(undef, nsets)
for set in 1:nsets
    data = readdlm(joinpath(ext_U, set_names[set], "$(tstep)_nodes.dat"), ',', Float64, '\n', skipstart=1)
    disp = Matrix{Union{Int,Float64}}(data)
    disp[:, 1] .= Int.(disp[:, 1])
    disp = sortslices(disp, dims=1, by=x->x[1])
    vals_nodes[set] = comp == "vertical" ? disp[:, end] :
                                           sqrt.(disp[:, 2].^2 .+ disp[:, 3].^2)
end

mycmap = comp == "vertical" ? cgrad(:balance) : cgrad(:bone_1, rev=true)
labl   = comp == "vertical" ? L"$u_\text{z}$ [m]" : L"$||\mathbf{u}_h||$ [m]"
ttl    = "u_$(comp == "vertical" ? "z" : "h")_ovbPico_$(sim_name)_$(p.year)"
fig_descriptor = (cmap=mycmap, cbar_label=labl, fig_size=(1500,1200), plot_height=1000,
                  cbar_height=400, plot_azi=p.azi, plot_elev=p.elev, plot_linew=0.1,
                  clim=p.clim,
                  xlms=(minimum(coords_map[1][:,1]), maximum(coords_map[1][:,1])),
                  ylms=(minimum(coords_map[1][:,2]), maximum(coords_map[1][:,2])),
                  save_fig=SAVE_TO_DEPOSIT, save_dir=figdir, title=ttl, target_dpi=400, target_size_in=(15,12))

plotNodalValues(face_list_map, coords_map, vals_nodes, fig_descriptor, line_data=coastline)

# Scratch output: coords + displacement at this panel
out_disp = joinpath(OUT_DIR, "displacement", sim_name)
isdir(out_disp) || mkpath(out_disp)
vtag = comp == "vertical" ? "u_z" : "u_h"
writedlm(joinpath(out_disp, "nodal_coords_and_$(vtag)_$(p.year).dat"), hcat(coords_map[1], vals_nodes[1]))