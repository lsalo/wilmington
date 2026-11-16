#=
Initial margin to Coulomb failure vs. DCFF on surfaces of the Wilmington model (plastic-120).
Sign convention as in the manuscript (tension positive; CFF = |τ| + μ σ'_n):
  M = -CFF(t0)                      M > 0 → below failure at t0
  failure reached where DCFF ≥ M    ⇔ CFF ≥ 0
Reads the deposited tractions (columns: id, σ'_n, |τ| [Pa]); t0 = first step, as in Figs 14, 15, 17.
=#
include(joinpath(@__DIR__, "..", "config.jl"))
using DelimitedFiles, Statistics, Printf
include(tool(:abaqus)); using .Abaqus      
include(tool(:fileio)); using .Fileio   
include(tool(:faultops));  using .FaultOperations   

## ---- Configuration -------------------------------------------------------------------
sim_name  = "plastic-120"
mu        = 0.6
nsets     = ["Tear_Fault_9"]   # e.g. "Wilmington", "Tear_Fault_9"
uh_filter = false                                  # Fig 17 filter (surface ‖u_h‖ ≥ 1 m); false for faults

ext_S     = joinpath(abaqus_extracted(sim_name), "S")
FAULT_DIR = datapath("for_fault_analysis", sim_name)

## ---- ‖u_h‖ ≥ 1 m bounding box (as in abaqus_results_surface_tractions.jl) -----------------
nodes, _, _, _, _ = get_connectivity_from_inp(GRIDS_ABAQUS, "WilmingtonMesh_AllFaults_Coarse.inp", "tetrahedral", false)
row_of = Dict(Int(nodes[i, 1]) => i for i in axes(nodes, 1))
uh = readdlm(joinpath(FAULT_DIR, "nodal_coords_and_u_h_1958.dat"))
above = uh[:, end] .>= 1.0
xmin, xmax = extrema(uh[above, 1]); ymin, ymax = extrema(uh[above, 2])

## ---- Margin vs. DCFF ----------------------------------------------------------------
@printf("%-18s %7s %10s %10s %10s %10s %10s\n", "nodeset", "nodes", "σn0_med", "τ0_med", "M_med", "DCFFmax", "frac_fail")
for nset in nsets
    fns = sort(readdir(joinpath(ext_S, nset)), by=sort_steps)
    tr  = [readdlm(joinpath(FAULT_DIR, "tractions", nset, "tractions_$(splitext(fn)[1]).dat"), '\t', Float64, '\n') for fn in fns]
    check_t0_tractions(tr[1], nset)
    keep = trues(size(tr[1], 1))
    if uh_filter
        xy = nodes[[row_of[Int(id)] for id in tr[1][:, 1]], 2:3]
        keep = (xmin .<= xy[:, 1] .<= xmax) .& (ymin .<= xy[:, 2] .<= ymax)
    end
    
    mu_req = vec(maximum(reduce(hcat, [x[keep, 3] ./ abs.(x[keep, 2]) for x in tr[2:end]]), dims=2))   # max_t |τ|/|σ'n|
    @printf("   %s: friction needed for slip — P50 %.3f, P95 %.3f, max %.3f\n", nset, quantile(mu_req, 0.5), quantile(mu_req, 0.95), maximum(mu_req))

    sn0, tau0 = tr[1][keep, 2], tr[1][keep, 3]
    M      = -(tau0 .+ mu .* sn0)                                                       # margin at t0
    dmax   = maximum(reduce(hcat, [(x[keep, 3] .+ mu .* x[keep, 2]) .+ M for x in tr[2:end]]), dims=2)  # max_t DCFF
    frac   = mean(vec(dmax) .>= M)                                                      # fraction reaching failure
    @printf("%-18s %7d %10.1f %10.1f %10.1f %10.1f %10.3f\n", nset, sum(keep),
            median(sn0)/1e5, median(tau0)/1e5, median(M)/1e5, maximum(dmax)/1e5, frac)
end
println("(stresses in bar; σn0 < 0 is compression)")