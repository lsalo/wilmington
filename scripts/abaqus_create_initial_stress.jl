#=
Initialize the Wilmington geomechanical stress state and write the Abaqus initial-stress
input files used in the paper.

Two stress fields are used in the manuscript (select via `field` below):
  - "a_fit"         : reverse-faulting, critically stressed, with depth-dependent A_phi(z)
                      from a best fit to the SCEC Community Stress Model (Hardebeck et al.,
                      2024). This is the plastic-130 field; reproduces Fig. 8a.
  - "heterogeneous" : strike-slip, sediments not critically stressed (mobilized mu=0.3,
                      A_phi=1.25) over a critically stressed basement (mu=0.6, A_phi=1.75).
                      This is the plastic-120 field; reproduces Fig. 8b.
The "normal", "strike-slip", and "reverse" options are exploratory homogeneous fields that
were tested but are NOT used in the paper; they are retained as options.

Sign convention (paper): tension positive. Pore pressure p is positive; total and effective stresses at
depth are compressive and therefore negative. The saved effective stresses are (s_i + p),
and Fig. 8 shows total stress magnitudes (pore pressure not plotted).

Abaqus/Wilmington grid convention: s11 // x, s22 // y, s33 // z.

A_phi provenance: with READ_CSM_DATA=true the script fits A_phi(z) from the SCEC
CSM .csv files (third-party, NOT redistributed in this deposit) and writes the two fit
coefficients to aphi_fit_coef.dat in the deposit. With READ_CSM_DATA=false (default) it reads
those deposited coefficients, so the field is regenerable without the CSM data. To re-derive,
download the CSM archive and set SCEC_DIR below.
=#
include(joinpath(@__DIR__, "..", "config.jl"))

using DelimitedFiles, Statistics, Printf, LinearAlgebra
using CairoMakie
CairoMakie.activate!()
include(tool(:abaqus)); using .Abaqus      # get_connectivity_from_inp, get_element_centroids_inp

## ── Options ──────────────────────────────────────────────────────────────────
field = "a_fit"                  # "a_fit" (plastic-130, Fig 8a) | "heterogeneous" (plastic-120, Fig 8b)
                                 #   | "normal" | "strike-slip" | "reverse"  (exploratory, not in paper)
const READ_CSM_DATA  = false     # true → fit A_phi(z) from SCEC CSM (needs SCEC_DIR); false → read deposited coef
const SAVE_TO_DEPOSIT = false    # true → write .inp/.dat/figures into the deposit
const SHOW_PLOT      = true      # build Fig. 8 panel for the selected field

# Local path to the SCEC CSM archive (only used when READ_CSM_DATA=true).
# The CSM .csv files are third-party and are NOT included in this deposit:
# download from the SCEC Community Stress Model and set this path.
SCEC_DIR = "/Users/lluis/MIT Dropbox/Lluis Salo/Postdoc_Harvard-MIT/research/LAB/Wilmington/2023_SCEC_CSM/SCEC_CSM_archive_2023/"   # <-- set your own path

fn_inp = "WilmingtonMesh_AllFaults_Rotated_Coarse.inp"
APHI_COEF_FILE = joinpath(FOR_ABAQUS_STRESS, "aphi_fit_coef.dat")

figdir = SAVE_TO_DEPOSIT ? datapath("figures", "abaqus", field == "a_fit" ? "plastic-130" : "plastic-120") : FIGS_DIR
isdir(figdir) || mkpath(figdir)
out_dir = SAVE_TO_DEPOSIT ? FOR_ABAQUS_STRESS : OUT_DIR
isdir(out_dir) || mkpath(out_dir)

## ── A_phi(z): fit from CSM (optional) or read deposited coefficients ──────────
# aphi_of_z takes depth-below-surface in km (positive) and returns A_phi
using CSV, GLM, StatsModels
if field == "a_fit"
    if READ_CSM_DATA
        # CSM stress models. NOTE: '#' was removed from the .csv header row (row 48) in the archive.
        filenames = ["FlatMaxwell", "Hardebeck_FM", "Luttrell-2017", "YH14-K", "YHSM-2013"]
        stress_data_all = Vector{}(undef, length(filenames))
        for (i, fn) in enumerate(filenames)
            sel = fn == "FlatMaxwell" ? (1:20) : [1, 2, 3, 14]
            stress_data_all[i] = CSV.File(joinpath(SCEC_DIR, fn * ".csv"), header=48, select=sel, skipto=49)
        end
        # Wilmington-area subset
        min_lon, max_lon = -118.41947, -117.93033
        min_lat, max_lat = 33.58679, 33.93842
        stress_data = Vector{}(undef, length(filenames))
        for i in 1:length(filenames)
            rows = [min_lon <= r.LON <= max_lon && min_lat <= r.LAT <= max_lat && 0 <= r.DEP <= 16 for r in stress_data_all[i]]
            stress_data[i] = Matrix{Float64}(undef, sum(rows), length(stress_data_all[i].names))
            k = 1
            for (j, r) in enumerate(stress_data_all[i])
                rows[j] && (stress_data[i][k, :] = [getproperty(r, n) for n in stress_data_all[i].names]; k += 1)
            end
        end
        # Linear A_phi vs depth fit across all five models (col 14 / col 4 = A_phi, col 3 = depth)
        xz = vcat([stress_data[i][:, 3] for i in 1:5]...)
        ya = vcat(stress_data[1][:, 14], stress_data[2][:, 4], stress_data[3][:, 4], stress_data[4][:, 4], stress_data[5][:, 4])
        lm_aphi = lm(@formula(y ~ x), (x=xz, y=ya))
        aphi_coef = coef(lm_aphi)                              # [intercept, slope]; A_phi = c1 + c2*depth_km
        SAVE_TO_DEPOSIT && writedlm(APHI_COEF_FILE, aphi_coef)
        @info "A_phi fit" intercept=aphi_coef[1] slope=aphi_coef[2]
    else
        aphi_coef = vec(readdlm(APHI_COEF_FILE))
    end
    aphi_of_z(depth_km) = aphi_coef[1] .+ aphi_coef[2] .* depth_km
end

## ── Mesh, centroids, depths ──────────────────────────────────────────────────
nodes_inp, elements_inp, elset_names_inp, _, _ = get_connectivity_from_inp(GRIDS_ABAQUS, fn_inp, "tetrahedral")
elements_all = vcat(elements_inp[1:end-1]...)
centroids = get_element_centroids_inp(elements_all, nodes_inp)

z      = centroids[:, end]                 # element depth, m (negative down)
z_nod  = nodes_inp[:, end]                 # nodal depth, m (negative down)
p      = -1000 * 9.807 .* z                # pore pressure, Pa (positive; tension-positive convention)
sv     =  2300 * 9.807 .* z                # vertical total stress, Pa (negative at depth)
p_nod  = -1000 * 9.807 .* z_nod
sv_nod =  2300 * 9.807 .* z_nod
t      = zeros(length(z))                  # zero shear components for the .inp

mu_crit   = 0.6
ratio_csm = ((mu_crit^2 + 1)^0.5 + mu_crit)^2   # effective-stress ratio (S1-p)/(S3-p)

## ── Build the stress field for the selected `field` ──────────────────────────
# Homogeneous fields (a_fit/reverse/normal/strike-slip) computed here; the heterogeneous
# field overwrites with per-material values in its own branch. Exactly one branch runs.
if field in ("a_fit", "reverse", "normal", "strike-slip")
    if field == "a_fit"
        a_phi     = aphi_of_z(-1e-3 .* z)
        a_phi_nod = aphi_of_z(-1e-3 .* z_nod)
        s3 = sv;       s1 = ratio_csm * (s3 .+ p) .- p;           s2 = (a_phi .- 2.0) .* (s1 .- s3) .+ s3
        s3_nod = sv_nod; s1_nod = ratio_csm*(s3_nod .+ p_nod) .- p_nod; s2_nod = (a_phi_nod .- 2.0).*(s1_nod .- s3_nod) .+ s3_nod
    elseif field == "reverse"
        a_phi = 2.5
        s3 = sv;       s1 = ratio_csm*(s3 .+ p) .- p;             s2 = (a_phi - 2.0).*(s1 .- s3) .+ s3
        s3_nod = sv_nod; s1_nod = ratio_csm*(s3_nod .+ p_nod) .- p_nod; s2_nod = (a_phi - 2.0).*(s1_nod .- s3_nod) .+ s3_nod
    elseif field == "normal"
        a_phi = 0.0
        s1 = sv;       s3 = (s1 .+ p)/ratio_csm .- p;             s2 = a_phi.*(s1 .- s3) .+ s3
        s1_nod = sv_nod; s3_nod = (s1_nod .+ p_nod)/ratio_csm .- p_nod; s2_nod = a_phi.*(s1_nod .- s3_nod) .+ s3_nod
    elseif field == "strike-slip"
        a_phi = 1.5
        A = a_phi + (-a_phi + 1)/ratio_csm - 2
        s2 = sv;       s1 = ((A + 1).*(-p) .- s2)./A;             s3 = (s1 .+ p)/ratio_csm .- p
        s2_nod = sv_nod; s1_nod = ((A + 1).*(-p_nod) .- s2_nod)./A; s3_nod = (s1_nod .+ p_nod)/ratio_csm .- p_nod
    end

elseif field == "heterogeneous"
    # Per-material strike-slip: sediments (mu=0.3, A_phi=1.25) over basement (mu=0.6, A_phi=1.75).
    names_sed = ["Ranger", "BottomRanger", "Pico", "AquiferFootwall", "OverburdenPico",
                 "OverburdenCenterSouth", "OverburdenCenterNorth", "OverburdenHangingWall",
                 "OverburdenRight", "OverburdenLeft"]
    id_sed = findall(x -> x in names_sed, elset_names_inp)
    id_bas = [3, 7, 8, 9, 10]
    elements_sed = vcat(elements_inp[id_sed]...);  id_el_sed = elements_sed[:, 1]
    elements_bas = vcat(elements_inp[id_bas]...);  id_el_bas = elements_bas[:, 1]

    nodes_sed_ids = sort(unique(vcat(elements_sed[:, 2:end]...)))
    id_nod_sed = findall(x -> x in nodes_sed_ids, nodes_inp[:, 1])
    nodes_bas_ids = sort(unique(vcat(elements_bas[:, 2:end]...)))
    id_nod_bas = findall(x -> x in nodes_bas_ids, nodes_inp[:, 1])
    id_nod_bas = filter(x -> !(x in id_nod_sed), id_nod_bas)   # basement-top nodes already counted as sediment

    a_phi = [1.25, 1.75]; mu = [0.3, 0.6]
    ratio = ((mu .^ 2 .+ 1) .^ 0.5 .+ mu) .^ 2
    A = a_phi .+ (-a_phi .+ 1) ./ ratio .- 2

    s2 = sv
    s1 = Vector{Float64}(undef, length(z)); s3 = similar(s1)
    s1[id_el_sed] = ((A[1] + 1) .* (-p[id_el_sed]) .- s2[id_el_sed]) ./ A[1]
    s1[id_el_bas] = ((A[2] + 1) .* (-p[id_el_bas]) .- s2[id_el_bas]) ./ A[2]
    s3[id_el_sed] = (s1[id_el_sed] .+ p[id_el_sed]) / ratio[1] .- p[id_el_sed]
    s3[id_el_bas] = (s1[id_el_bas] .+ p[id_el_bas]) / ratio[2] .- p[id_el_bas]

    s2_nod = sv_nod
    s1_nod = Vector{Float64}(undef, length(z_nod)); s3_nod = similar(s1_nod)
    s1_nod[id_nod_sed] = ((A[1] + 1) .* (-p_nod[id_nod_sed]) .- s2_nod[id_nod_sed]) ./ A[1]
    s1_nod[id_nod_bas] = ((A[2] + 1) .* (-p_nod[id_nod_bas]) .- s2_nod[id_nod_bas]) ./ A[2]
    s3_nod[id_nod_sed] = (s1_nod[id_nod_sed] .+ p_nod[id_nod_sed]) / ratio[1] .- p_nod[id_nod_sed]
    s3_nod[id_nod_bas] = (s1_nod[id_nod_bas] .+ p_nod[id_nod_bas]) / ratio[2] .- p_nod[id_nod_bas]
else
    error("Unrecognized field: $field")
end

## ── Nodal mean effective stress (for depth-dependent hardening) ───────────────
peff_ini_nod = vec(-(1/3) * sum(hcat(s1_nod .+ p_nod, s2_nod .+ p_nod, s3_nod .+ p_nod), dims=2))

## ── Write Abaqus initial-stress files (EFFECTIVE stresses, tension positive) ──
# Column order maps principal stresses to grid axes per regime:
#   a_fit/reverse: s11=sxx=s1, s22=syy=s2, s33=szz=s3
#   heterogeneous (strike-slip): s11=sxx=s1, s22=syy=s3, s33=szz=s2
if field in ("a_fit", "reverse")
    stress_out = hcat(centroids[:, 1], s1 .+ p, s2 .+ p, s3 .+ p, t, t, t)
    aphi_tag = field == "a_fit" ? round(mean(aphi_of_z(-1e-3 .* z)), digits=2) : 2.06
    fn_stress = "initial_stress_rotated_criticallyStressed_mu0.6_aphi$(aphi_tag).inp"
    fn_peff   = "peff_ini_nod_aphi$(aphi_tag).inp"
elseif field == "normal"
    stress_out = hcat(centroids[:, 1], s2 .+ p, s3 .+ p, s1 .+ p, t, t, t)
    fn_stress = "initial_stress_rotated_criticallyStressed_mu0.6_normal.inp"
    fn_peff   = "peff_ini_nod_normal.inp"
elseif field in ("strike-slip", "heterogeneous")
    stress_out = hcat(centroids[:, 1], s1 .+ p, s3 .+ p, s2 .+ p, t, t, t)
    if field == "heterogeneous"
        fn_stress = "initial_stress_rotated_sed_mu0.3_aphi1.25_bas_mu0.6_aphi1.75.inp"
        fn_peff   = "peff_ini_nod_sed_mu0.3_aphi1.25_bas_mu0.6_aphi1.75.inp"
    else
        fn_stress = "initial_stress_rotated_strikeslip_aphi1.5.inp"
        fn_peff   = "peff_ini_nod_strikeslip_aphi1.5.inp"
    end
end

if SAVE_TO_DEPOSIT
    writedlm(joinpath(out_dir, fn_stress), stress_out, ',')
    writedlm(joinpath(out_dir, fn_peff), hcat(nodes_inp[:, 1], map(x -> @sprintf("%.2f", x), peff_ini_nod)), ",\t")
    # Geometry / index products used downstream
    writedlm(joinpath(out_dir, "element_centroids.dat"), centroids, ',')
    writedlm(joinpath(out_dir, "nodal_depth.dat"), hcat(nodes_inp[:, 1], nodes_inp[:, end]), ',')
    writedlm(joinpath(out_dir, "id_ranger_pico.dat"),          vcat(elements_inp[1][:, 1], elements_inp[4][:, 1]))
    writedlm(joinpath(out_dir, "id_ranger_pico_overPico.dat"), vcat(elements_inp[1][:, 1], elements_inp[4][:, 1], elements_inp[6][:, 1]))
    writedlm(joinpath(out_dir, "id_pico.dat"),          elements_inp[4][:, 1])
    writedlm(joinpath(out_dir, "id_overPico.dat"),      elements_inp[6][:, 1])
    writedlm(joinpath(out_dir, "id_bottomRanger.dat"),  elements_inp[2][:, 1])

    id_above_bas = findall(name -> !occursin("basement", lowercase(name)) && name != "free_elements", elset_names_inp)
    id_all_above_basement = vcat([elements_inp[i][:, 1] for i in id_above_bas]...)
    writedlm(joinpath(out_dir, "id_all_above_basement.dat"), id_all_above_basement)
end

## ── Fig. 8: initial total stress magnitudes vs depth (central column) ─────────
if SHOW_PLOT
    using Distances
    xc, yc = mean(nodes_inp[:, 2]), mean(nodes_inp[:, 3])
    target = hcat(fill(xc, 76), fill(yc, 76), LinRange(0, -15000, 76))
    D = pairwise(Euclidean(), centroids[:, 2:end], target; dims=1)
    col = vec(mapslices(argmin, D; dims=1))                      # nearest element down the central column
    zc_km = centroids[col, end] ./ -1000

    f = Figure(size=(400, 300))
    ax = Axis(f[1, 1], xlabel=L"\text{Total stress magnitude [MPa]}", ylabel=L"\text{Depth [km]}",
              xlabelsize=16, ylabelsize=16,
              xticks=[0, -100, -200, -300, -400, -500, -600, -700],
              yticks=[0, 2.5, 5, 7.5, 10, 12.5, 15])
    if field == "heterogeneous"
        # panel b: sediment/basement annotations + representative boundary
        poly!(ax, Point2f.([0, -700, -700, 0, 0], [15, 15, 2.9, 2.9, 15]), color=RGBAf(0.8, 0.8, 0.8, 0.5))
        lines!(ax, [-700, 0], [2.9, 2.9], color=RGBAf(0.5, 0.5, 0.5, 1), linewidth=0.75, linestyle=:dash)
        text!(ax, -200, 2.0, text=L"\text{Sediments, mob. } \mu=0.3 \text{, } A_\phi=1.25", fontsize=14)
        text!(ax, -200, 5.0, text=L"\text{Basement, mob. } \mu=0.6 \text{, } A_\phi=1.75", fontsize=14)
    end
    # Which principal is vertical differs by regime:
    #   a_fit/reverse:  s3 = sigma_v (min);  s2 = sigma_h;  s1 = sigma_H
    #   heterogeneous:  s2 = sigma_v (int);  s3 = sigma_h;  s1 = sigma_H
    if field == "heterogeneous"
        sH, sh, svert = s1, s3, s2
    else  # a_fit / reverse
        sH, sh, svert = s1, s2, s3
    end
    lines!(ax, sH[col]    ./ 1e6, zc_km, color=RGBAf(0, 0, 0, 1), linewidth=0.5, label=L"\sigma_\text{H}")
    lines!(ax, sh[col]    ./ 1e6, zc_km, color=RGBAf(0, 0, 1, 1), linewidth=0.5, label=L"\sigma_\text{h}")
    lines!(ax, svert[col] ./ 1e6, zc_km, color=RGBAf(1, 0, 0, 1), linewidth=0.5, label=L"\sigma_\text{v}")
    ylims!(ax, 15, 0); xlims!(ax, 0, -700)
    axislegend(position=:lb, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1, 1, 1, 0.8))
    display(f)
    if SAVE_TO_DEPOSIT
        panel = field == "a_fit" ? "fig8a_initial_stress_reverse" : "fig8b_initial_stress_heterogeneous"
        save(joinpath(figdir, "$panel.pdf"), f, pt_per_unit=1)
    end
end