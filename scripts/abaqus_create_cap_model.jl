#=
Build the Modified Drucker-Prager Cap Model hardening curves for the Wilmington
geomechanical model and reproduce Fig. 9.

For each variant the hydrostatic yield p_b is computed at a set of reference depths from
the initial effective pressure (Eq. 20 in the paper): p_b = f_p * p_eff + f_K * volumetric
hardening term. The resulting hardening table is written as an Abaqus-ready .inp include
file (the contents are pasted into the simulation deck). The same hardening curve is used
for all units above the basement ("above-basement").

  plastic-130 (reverse, CSM):          R=0.20, f_p(=a)=1.18, f_K(=b)=0.25
  plastic-120 (strike-slip, heterog.): R=0.80, f_p(=a)=1.17, f_K(=b)=0.47

Fig. 9 = the strike-slip (plastic-120) hardening curves, p'_b vs volumetric plastic strain.

Sign convention: tension positive; p_eff and p_b are compression-positive (consistent
with peff_ini_nod*.inp). Stresses read in Pa, plotted/written in MPa where noted.
=#
include(joinpath(@__DIR__, "..", "config.jl"))

using DelimitedFiles, Printf, Interpolations
using CairoMakie
CairoMakie.activate!()

## ── Options ──────────────────────────────────────────────────────────────────
variant = "plastic-130"             # "plastic-120" (Fig. 9) | "plastic-130"
SAVE_TO_DEPOSIT = true              # true → write hardening .inp + Fig. 9 into the deposit
SHOW_PLOT = true                    # build Fig. 9 (only meaningful for plastic-120)

# Per-variant cap parameters and the initial-stress files they read.
CAP_CFG = Dict(
    "plastic-130" => (R=0.20, f_p=1.18, f_K=0.25, aphi_tag="CSMfit",
        stress_fn="initial_stress_rotated_criticallyStressed_mu0.6_aphi2.07.inp",
        peff_fn  ="peff_ini_nod_aphi2.07.inp"),
    "plastic-120" => (R=0.80, f_p=1.17, f_K=0.47, aphi_tag="1.25-1.75",
        stress_fn="initial_stress_rotated_sed_mu0.3_aphi1.25_bas_mu0.6_aphi1.75.inp",
        peff_fn  ="peff_ini_nod_sed_mu0.3_aphi1.25_bas_mu0.6_aphi1.75.inp"),
)
cfg = CAP_CFG[variant]

# Shared Modified Drucker-Prager Cap parameters (variant-independent)
const D_COH = 1e-3      # cohesion d [MPa]
const BETA  = 55.0      # friction angle beta [deg]  (phi = 35)
const W_CAP = 0.27      # W [-]
const D_CAP = 0.7143    # D [1/kbar]
# Reference depths (pairs bracketing each unit) and plastic-strain values for above-basement
Z_PT = -[5, 6, 890, 900, 1420, 1430, 2890, 3000, 4890, 5000, 8890, 9000]
E_P  = [0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]

out_dir = SAVE_TO_DEPOSIT ? FOR_ABAQUS_STRESS : OUT_DIR
isdir(out_dir) || mkpath(out_dir)
figdir = SAVE_TO_DEPOSIT ? datapath("figures", "abaqus", variant) : FIGS_DIR
isdir(figdir) || mkpath(figdir)

## ── Read initial stress / effective pressure / geometry ──────────────────────
nodal_depth  = readdlm(joinpath(FOR_ABAQUS_STRESS, "nodal_depth.dat"), ',')
peff_ini_nod = readdlm(joinpath(FOR_ABAQUS_STRESS, cfg.peff_fn), ',')
# (initial stress is read for the p-q scatter in the cap-envelope figure, not needed for Fig. 9;
#  kept minimal here since Fig. 3 is drawn elsewhere.)

## ── Hardening curve p_b at each reference depth (Eq. 20) ─────────────────────
nc = Int(length(Z_PT) / 2)
peff_id = Vector{Float64}(undef, nc)
# Pre-build a sorted, unique depth→peff interpolant for depths with no exact node match
id_z      = unique(i -> nodal_depth[i, 2], eachindex(nodal_depth[:, 2]))
z_unique  = nodal_depth[id_z, 2]
peff_uniq = peff_ini_nod[id_z, 2]
sp        = sortperm(z_unique)
itp       = interpolate((z_unique[sp],), peff_uniq[sp], Gridded(Linear()))
for i in 1:nc
    id = findfirst(x -> Z_PT[2i] <= x <= Z_PT[2i-1], nodal_depth[:, 2])
    peff_id[i] = isnothing(id) ? itp(Z_PT[2i-1]) : peff_ini_nod[id, end]
end
# Hardening table: p_b as a function of E_P for each reference depth
# p_b(e) = f_p*peff + f_K * (1/3)(1/D) ln(e/W + 1) * 1e8   [Pa]
hard_pb(peff, e) = cfg.f_p * peff + cfg.f_K * ((1/3)*(1/D_CAP)*log(e/W_CAP + 1)*1e8)
p_b_table = [hard_pb(peff_id[i], e) for i in 1:nc for e in E_P]    # length nc*length(E_P), Pa

## ── Write Abaqus hardening .inp (paste into the simulation deck) ──────────────
fn_hardening = "hardening_curve_$(cfg.f_p)peff_$(cfg.R)R_$(cfg.f_K)Kosloff_above-basement_aphi$(cfg.aphi_tag).inp"
p_b_str   = map(x -> @sprintf("%.2f", x), p_b_table)
peff_col  = repeat(map(x -> @sprintf("%.2f", x), peff_id), inner=length(E_P))
to_save   = hcat(p_b_str, repeat(E_P, nc), zeros(nc*length(E_P)), peff_col)  # decimal format (no 1e7) for Abaqus
SAVE_TO_DEPOSIT && writedlm(joinpath(out_dir, fn_hardening), to_save, ",\t")

## ── Fig. 9: hardening curves for the strike-slip model (plastic-120) ─────────
# p_b at full hardening (e = E_P) for representative reference depths.
pb_curve(i) = [hard_pb(peff_id[i], e) for e in E_P] .* 1e-6      # MPa
if SHOW_PLOT && variant == "plastic-120"
    sel = [6, 5, 2, 1]                                          # 0.06, 10.5, 57.2, 140 MPa (range of p')
    colors = reverse([RGBAf(0,0,0,1), RGBAf(0,0,1,1), RGBAf(0,1,0,1), RGBAf(1,0,0,1)])
    f = Figure(size=(450, 350))
    ax = Axis(f[1,1], xlabel=L"e^\text{pl}_\text{v} \; \text{ [-]}", ylabel=L"p^\prime_\text{b} \; \text{ [MPa]}",
              xlabelsize=16, ylabelsize=16, xticks=0:0.05:0.2, yticks=0:50:200,
              title=L"A_\phi = 1.25 \; | \; 1.75", titlesize=16)
    for (k, i) in enumerate(sel)
        pcurve = pb_curve(i)
        lines!(ax, E_P, pcurve, color=colors[k], linewidth=1.5,
               label=L"p^\prime = %$(round(peff_id[i]*1e-6, digits=2)) \; \text{MPa}")
    end
    xlims!(ax, 0, 0.2); ylims!(ax, 0, 200)
    axislegend(position=:lc, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8))
    display(f)
    SAVE_TO_DEPOSIT && save(joinpath(figdir, "fig9_hardening_curves_$(variant).pdf"), f, pt_per_unit=1)
end

@info "Cap hardening built" variant fn_hardening pb_at_reference_depths_MPa=round.(peff_id .* 1e-6, digits=3)