#=
Stress path of a Ranger reservoir element near the max-subsidence node, its
effective-stress history, and surface subsidence vs. field data (Kosloff et al., 1980).
=#
include(joinpath(@__DIR__, "..", "config.jl"))

using LinearAlgebra                          # norm
using HDF5                                   # element stress reads
using DelimitedFiles
using CairoMakie
CairoMakie.activate!()
include(tool(:fileio)); using .Fileio        # sort_steps
include(tool(:abaqus)); using .Abaqus        # get_connectivity_from_inp

## Variant + coupled options + figures
# Per-variant configuration (single source of truth for the coupled options)
const VARIANT_CFG = Dict(
    "plastic-130" => (regime = "reverse",       # CSM, critically stressed
                      stress_ini = "initial_stress_rotated_criticallyStressed_mu0.6.inp",
                      pb = (1.18, 0.25),
                      q_Fc = vcat(Vector{}(0:1:29), 29.01:0.01:30.1),
                      R = 0.2,
                      ylm = -50,
                      lbl = L"A_\phi = 2.06",
                      lbl_pos = (1980, -48)),
    "plastic-120" => (regime = "strike-slip",   # heterogeneous, sediments not critically stressed
                      stress_ini = "initial_stress_rotated_sed_mu0.3_aphi1.25_bas_mu0.6_aphi1.75.inp",
                      pb = (1.17, 0.47),
                      q_Fc = vcat(Vector{}(0:1:7), 7.01:0.01:9, 9:.001:9.355),
                      R = 0.8,
                      ylm = -30,
                      lbl = L"A_\phi = 1.25 \, | \, 1.75", 
                      lbl_pos = (1970, -29))
)

# Select variant here
sim_name = "plastic-120"            
cfg      = VARIANT_CFG[sim_name]
node_select = "max_subsidence"
fn_inp   = "WilmingtonMesh_AllFaults_Rotated_Coarse.inp"
ext = abaqus_extracted(sim_name)             # .../simulation_output/abaqus/<variant>/extracted

abaqus_dir_el  = joinpath(ext, "S", "ENTIRE_DOMAIN")       # element stresses (h5)
abaqus_dir_nod = joinpath(ext, "U", "TOPRANGERSURFACE")    # nodal disp, top Ranger
abaqus_dir_top = joinpath(ext, "U", "TOPBC")               # nodal disp, top boundary
abaqus_fn_el  = sort(readdir(abaqus_dir_el),  by=sort_steps)
abaqus_fn_nod = sort(readdir(abaqus_dir_nod), by=sort_steps)
abaqus_fn_top = sort(readdir(abaqus_dir_top), by=sort_steps)

# Set true to (re)generate the canonical figures shipped in the dataset; false → scratch FIGS_DIR.
const SAVE_TO_DEPOSIT = true
figdir = SAVE_TO_DEPOSIT ? datapath("figures", "abaqus", sim_name) : FIGS_DIR
isdir(figdir) || mkpath(figdir)
savefig(name) = joinpath(figdir, name)

## Time axis — length must equal the extracted step count (see flag)
tn = 84
t = 1936 .+ vcat(-1, 0, cumsum(repeat([1], inner=tn)))     # [1935 (geostatic), 1936, ..., 2020]
@assert length(abaqus_fn_el) == length(abaqus_fn_nod) == length(abaqus_fn_top) == length(t)  "extracted step counts must match each other and length(t)"

## Max-subsidence node at the final top-Ranger step, tracked through time
u_topranger_fin = readdlm(joinpath(abaqus_dir_nod, abaqus_fn_nod[end]), ',', skipstart=1)
id_u  = node_select == "max_subsidence" ? argmin(u_topranger_fin[:, end]) : argmax(u_topranger_fin[:, end])
nod_u = Int(u_topranger_fin[id_u, 1])

u_z = Vector{Float32}(undef, length(abaqus_fn_nod))
for i in eachindex(abaqus_fn_nod)
    u_z[i] = readdlm(joinpath(abaqus_dir_nod, abaqus_fn_nod[i]), ',', skipstart=1)[id_u, end]
end

## Min/max surface displacement (node may differ per step)
u_z_top = Vector{Float32}(undef, length(abaqus_fn_top))
for i in eachindex(abaqus_fn_top)
    col = readdlm(joinpath(abaqus_dir_top, abaqus_fn_top[i]), ',', skipstart=1)[:, end]
    u_z_top[i] = node_select == "max_subsidence" ? minimum(col) : maximum(col)
end

## Ranger element nearest the selected node (assumes Ranger element ids start at 1, sorted)
nodes_inp, elements_inp, _, _, _ = get_connectivity_from_inp(GRIDS_ABAQUS, fn_inp, "tetrahedral")
idx = vec(any(elements_inp[1][:, 2:5] .== nod_u, dims=2))
id_elements = elements_inp[1][idx, 1]

centroids     = readdlm(joinpath(FOR_ABAQUS_STRESS, "element_centroids.dat"), ',')
centroids_els = centroids[id_elements, 2:end]
coord_nod_u   = nodes_inp[findfirst(nodes_inp[:, 1] .== nod_u), 2:end]
dist          = mapslices(row -> norm(row - coord_nod_u), centroids_els, dims=2)
id_element    = id_elements[argmin(dist)]

## Effective-stress path of that element
n_el  = length(abaqus_fn_el)
s     = Vector{Vector{Float64}}(undef, n_el)
p_eff = Vector{Float32}(undef, n_el)
q     = Vector{Float32}(undef, n_el)
for i in 1:n_el
    h5open(joinpath(abaqus_dir_el, abaqus_fn_el[i]), "r") do file
        abq  = hcat([read(file[name]) for name in ("Sxx","Syy","Szz","Sxy","Sxz","Syz")]...)
        s[i] = abq[id_element, :]
    end
    p_eff[i] = -(1/3)*sum(s[i][1:3])
    q[i]     = sqrt(0.5*((s[i][1]-s[i][2])^2 + (s[i][2]-s[i][3])^2 + (s[i][3]-s[i][1])^2) + 6*(s[i][4]^2 + s[i][5]^2 + s[i][6]^2))
end


## -------------------------------------- PLOTS ---------------------------------------------------------
# Stress paths

# Plot parameters
# Cap model parameters
d = 1e-3        # MPa
b = 55          # 55, 58.6, 52.41 deg
W = 0.27        # [-]
D = 0.7143      # 1/kbar
e_p = [0, 0.0025, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2]

# Hardening curve [Pa]
s_ini_el = readdlm(joinpath(FOR_ABAQUS_STRESS, cfg.stress_ini), ',')
p_eff_ini_el = -(1/3)*sum(s_ini_el[id_element, 2:4])
p_b = cfg.pb[1]*p_eff_ini_el .+ cfg.pb[2]*((1/3)*(1/D)*log.(e_p./W .+ 1)*1e8)

# Shear envelope and cap
e_p0 = 0
p_Fs = Vector{}(0:10:50)
q_Fs = p_Fs.*tand(b) .+ d
p_a = (1e-6*p_b[1] - cfg.R*d) / (1 +cfg.R*tand(b))             # pb[1] = Initial cap position
q_Fc = cfg.q_Fc
p_Fc = p_a .+ sqrt.( abs.((cfg.R*(d + p_a*tand(b)))^2 .- (cfg.R.*q_Fc).^2) )

# Make plot stress path
f = Figure(size=(800, 350))
ax = Axis(f[1, 1], xlabel = L"p^\text{eff} = -\frac{\sigma_{kk}^\text{eff}}{3} \; \text{ [MPa]}", ylabel = L"q = \sqrt{3 J_2} \; \text{ [MPa]}", xlabelsize = 16, ylabelsize = 16)
ln1 = lines!(ax, p_Fs, q_Fs, color = RGBAf(0, 0, 0, 1), linewidth = 1.5, label=L"F_\text{s},\; \beta=55.0^\circ")
ln1 = lines!(ax, p_Fc, q_Fc, color = RGBAf(0.5, 0.5, 0.5, 1), linewidth = 1, linestyle = :solid, label=L"F_\text{c}^0")
sc1 = scatter!(ax, p_eff/1e6, q/1e6, markersize = 6, color=t, colormap=:rainbow, label=L"\text{Ranger element}")
if cfg.regime == "reverse"
    xlims!(ax, 0, 35)
    ylims!(ax, 0, 35)
elseif cfg.regime == "normal"
    xlims!(ax, 0, 16)
    ylims!(ax, 0, 16)
elseif cfg.regime == "strike-slip"
    xlims!(ax, 0, 16)
    ylims!(ax, 0, 16)
end
axislegend(position = :lt, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8), labelsize=14) 
cbar = Colorbar(f[1, 2], sc1, label=L"\text{Year}", labelsize=12, width=12, ticklabelsize=12)

if cfg.regime == "reverse"
    ax = Axis(f[1, 3], xlabel = L"p^\text{eff} = -\frac{\sigma_{kk}^\text{eff}}{3} \; \text{ [MPa]}", ylabel = L"u_z \; \text{ [m]}", xlabelsize = 16, ylabelsize = 16, xticks=18:2:30, yticks=-10:2:0)
    sc1 = scatter!(ax, p_eff/1e6, u_z, markersize = 8, color=1:length(u_z), colormap=:rainbow, label=L"\max(u_z) \text{, top Ranger}")
    if node_select == "max_subsidence"
        xlims!(ax, 18, 30)
        ylims!(ax, -10, 0)
    elseif node_select == "max_rebound"
        xlims!(ax, -10, -20)
        ylims!(ax, -1, 5)
    end
elseif cfg.regime == "normal"
    ax = Axis(f[1, 3], xlabel = L"p^\text{eff} = -\frac{\sigma_{kk}^\text{eff}}{3} \; \text{ [MPa]}", ylabel = L"u_z \; \text{ [m]}", xlabelsize = 16, ylabelsize = 16, xticks=0:2:16, yticks=-10:2:0)
    sc1 = scatter!(ax, p_eff/1e6, u_z, markersize = 6, color=1:length(u_z), colormap=:rainbow, label=L"\max(u_z) \text{, top Ranger}")
    if node_select == "max_subsidence"
        xlims!(ax, 0, 16)
        ylims!(ax, -10, 0)
    elseif node_select == "max_rebound"
        xlims!(ax, 5, 15)
        ylims!(ax, -1, 5)
    end
elseif cfg.regime == "strike-slip"
    ax = Axis(f[1, 3], xlabel = L"p^\text{eff} = -\frac{\sigma_{kk}^\text{eff}}{3} \; \text{ [MPa]}", ylabel = L"u_z \; \text{ [m]}", xlabelsize = 16, ylabelsize = 16, xticks=0:2:16, yticks=-10:2:0)
    sc1 = scatter!(ax, p_eff/1e6, u_z, markersize = 8, color=1:length(u_z), colormap=:rainbow, label=L"\max(u_z) \text{, top Ranger}")
    if node_select == "max_subsidence"
        xlims!(ax, 0, 16)
        ylims!(ax, -10, 0)
    elseif node_select == "max_rebound"
        xlims!(ax, 5, 15)
        ylims!(ax, -1, 5)
    end
end
axislegend(position = :lb, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8), labelsize=14) 
f
save(savefig("stress_path_ranger_$(node_select)_$(sim_name).pdf"), f, pt_per_unit=1)

# Stress ratios figure
s_xx = [x[1] for x in s]
s_yy = [x[2] for x in s]
s_zz = [x[3] for x in s]
s_xx_norm = s_xx ./ s[1][1]
s_yy_norm = s_yy ./ s[1][2]
s_zz_norm = s_zz ./ s[1][3]

f = Figure(size=(700, 300))
ax = Axis(f[1, 1], xlabel = L"t \; \text{ [y]}", ylabel = L"\sigma_{ii}^{\prime} \; \text{[MPa]}", xlabelsize = 18, ylabelsize = 18, xticks=1940:20:2020)
lines!(ax, [1958, 1958], [-50, 0.0], color = RGBAf(0.5, 0.5, 0.5, 1), linewidth = 1.0, linestyle = :dash)
ln1 = lines!(ax, t, s_xx ./ 1e6, color = RGBAf(0, 0, 0, 1), linewidth = 1.0, label=L"\sigma_{xx}")
ln2 = lines!(ax, t, s_yy ./ 1e6, color = RGBAf(1, 0, 0, 1), linewidth = 1.0, label=L"\sigma_{yy}")
ln3 = lines!(ax, t, s_zz ./ 1e6, color = RGBAf(0, 0, 1, 1), linewidth = 1.0, label=L"\sigma_{zz}")
text!(ax, L"\text{Depletion}", position = (1937, -5))
#text!(ax, L"\text{Repr.}", position = (1959, -5))
text!(ax, cfg.lbl, position = cfg.lbl_pos)
xlims!(ax, 1935, 2020)
ylims!(ax, cfg.ylm, 0)
axislegend(position = :rb, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8), labelsize=16)

ax = Axis(f[1, 2], xlabel = L"t \; \text{ [y]}", ylabel = L"\sigma_{ii}^{\prime}\, / \, \sigma_{ii,0}^{\prime} \; \text{[-]}", xlabelsize = 18, ylabelsize = 18, xticks=1940:20:2020)
lines!(ax, [1958, 1958], [0.5, 2.0], color = RGBAf(0.5, 0.5, 0.5, 1), linewidth = 1.0, linestyle = :dash)
ln1 = lines!(ax, t, s_xx_norm, color = RGBAf(0, 0, 0, 1), linewidth = 1.0, label=L"\overline{\sigma}_{xx}")
ln2 = lines!(ax, t, s_yy_norm, color = RGBAf(1, 0, 0, 1), linewidth = 1.0, label=L"\overline{\sigma}_{yy}")
ln3 = lines!(ax, t, s_zz_norm, color = RGBAf(0, 0, 1, 1), linewidth = 1.0, label=L"\overline{\sigma}_{zz}")
xlims!(ax, 1935, 2020)
ylims!(ax, 0.5, 2.0)
f
save(savefig("stress_ratios_ranger_$(node_select)_$(sim_name).pdf"), f, pt_per_unit=1)

# 
# Plot 2: max(u_z) vs t, comparison with subsidence data
#
data_u3max = [1947 1951 1954 1958 1965 1975; -2.4 -4 -6.1 -7.9 -8.8 -7.4]        # Kosloff et al., Tectonophysics (1980)
#data_u3max = [1975; 1.4]  # 1975-1965, 1.4 m rebound, Kosloff et al., Tectonophysics (1980)
f = Figure(size=(400, 300))
#ax = Axis(f[1, 1], xlabel = L"t \; \text{ [y]}", ylabel = L"\max \; u_z \; \text{ [m]}", xlabelsize = 18, ylabelsize = 18, xticks=1930:10:2020, yticks =0:2:10)
ax = Axis(f[1, 1], xlabel = L"t \; \text{ [y]}", ylabel = L"\max \; u_z \; \text{ [m]}", xlabelsize = 18, ylabelsize = 18, xticks=1930:10:2020, yticks =-10:2:0)
sc1 = scatter!(ax, data_u3max[1, :], data_u3max[2, :], markersize=10, color=RGBAf(1,1,1,0), marker=:rect, strokewidth=1, strokecolor=:black, label=L"\text{ Field data}")
ln1 = lines!(ax, t, u_z_top, color=:black, linewidth=1, label=L"\text{ Model}")
xlims!(ax, 1930, 2020)
ylims!(ax, -10, 0)
#ylims!(ax, 0, 10)
axislegend(position = :rt, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8)) 
f
save(savefig("u3_surface_$(sim_name)_subsidence.pdf"), f, pt_per_unit=1)