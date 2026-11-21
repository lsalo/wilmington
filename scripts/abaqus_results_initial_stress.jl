#=
Sanity check: stresses resolved after the geostatic step should match the initial
conditions imposed (per variant). Produces a 3-panel I.C.-vs-Geostatic scatter for
each model variant; not a paper figure, deposited as evidence the initialization is correct.
=#
include(joinpath(@__DIR__, "..", "config.jl"))

using DelimitedFiles, LaTeXStrings, Printf
using CairoMakie
CairoMakie.activate!()

sim_name = "plastic-130"     # ← switch variant

# initial-stress file per variant (same mapping as the reservoir script's VARIANT_CFG)
const STRESS_INI = Dict(
    "plastic-130" => "initial_stress_rotated_criticallyStressed_mu0.6_aphi2.07.inp",
    "plastic-120" => "initial_stress_rotated_sed_mu0.3_aphi1.25_bas_mu0.6_aphi1.75.inp",
)

# Per variant: physical component → (I.C. column in stress_ini, geostatic column in abq_out, label, color)
# abq_out is read by name (Sxx,Syy,Szz) so its columns are always physical 1,2,3.
# stress_ini stores s_11,s_22,s_33 in cols 2,3,4 — physical order differs by variant.
const COMP_MAP = Dict(
    "plastic-130" => [(2, 1, L"\sigma^\text{eff}_{xx}", :black),   # s_11 = xx
                      (3, 2, L"\sigma^\text{eff}_{yy}", :red),     # s_22 = yy
                      (4, 3, L"\sigma^\text{eff}_{zz}", :blue)],   # s_33 = zz
    "plastic-120" => [(2, 1, L"\sigma^\text{eff}_{xx}", :black),   # s_11 = xx
                      (3, 2, L"\sigma^\text{eff}_{yy}", :red),     # s_33 = yy
                      (4, 3, L"\sigma^\text{eff}_{zz}", :blue)],   # s_22 = zz
)

const SAVE_TO_DEPOSIT = false
figdir = SAVE_TO_DEPOSIT ? datapath("figures", "abaqus", sim_name) : FIGS_DIR
isdir(figdir) || mkpath(figdir)

# Geostatic stress (resolved) and imposed initial stress
geo_file = joinpath(abaqus_extracted(sim_name), "S", "ENTIRE_DOMAIN", "Geostatic_elements.dat")
abq_raw  = readdlm(geo_file, ',', Float64, '\n', skipstart=1)       # ELEMENT_ID, Sxx, Syy, Szz, Sxy, Sxz, Syz
stress_ini = readdlm(joinpath(FOR_ABAQUS_STRESS, STRESS_INI[sim_name]), ',')

# align both to ascending element-id order
abq_out    = abq_raw[sortperm(abq_raw[:, 1]), 2:7]                  # columns: Sxx, Syy, Szz, Sxy, Sxz, Syz
stress_ini = stress_ini[sortperm(stress_ini[:, 1]), :]
@assert abq_raw[sortperm(abq_raw[:, 1]), 1] == stress_ini[:, 1] "element ids differ"

# 3-panel comparison: imposed (I.C.) vs resolved (Geostatic), σ11 / σ22 / σ33
f = Figure(size=(1200, 420))
cmap = COMP_MAP[sim_name]
panels = COMP_MAP[sim_name]
for (p, (col_ini, col_geo, sym, c)) in enumerate(panels)
    x = 1e-6 .* stress_ini[:, col_ini]
    y = 1e-6 .* abq_out[:, col_geo]
    lo = min(minimum(x), minimum(y)); hi = max(maximum(x), maximum(y))
    pad = 0.03 * (hi - lo)
    lims = (lo - pad, hi + pad)            # shared x/y range so the box is square and frames all data
    
    xlab = latexstring(sym.s * raw" \text{ (I.C.)} \; [\text{MPa}]")
    ylab = latexstring(sym.s * raw" \text{ (Geostatic)} \; [\text{MPa}]")
    ax = Axis(f[1, p], xlabel=xlab, ylabel=ylab, xlabelsize=16, ylabelsize=16, aspect=1, xreversed=true)
    scatter!(ax, 1e-6 .* stress_ini[:, col_ini], 1e-6 .* abq_out[:, col_geo],
             markersize=6, color=RGBAf(0,0,0,0), strokecolor=c, strokewidth=0.5, label=L"\text{ Data}")
    lines!(ax, collect(lims), collect(lims), color=RGBAf(0.8,0.8,0.8,1), linewidth=1.25, label=L" 1:1")
    xlims!(ax, 0, lims[1]); ylims!(ax, lims...)
    p == 1 && axislegend(position=:lt, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8))

   res = y .- x                                  # geostatic − I.C., MPa (both already ×1e-6)
    @printf("%-6s  bias = %+7.3f   std = %6.3f   P95|res| = %6.3f  MPa\n",
            sym.s, mean(res), std(res), quantile(abs.(res), 0.95))
end
display(f)

ppu = 9 * 300 / 900
SAVE_TO_DEPOSIT && save(joinpath(figdir, "ic_vs_geostatic_$(sim_name).png"), f, px_per_unit=ppu)