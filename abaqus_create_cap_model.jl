#=
Script to plot modified Drucker-Prager cap model and initial reservoir state.
=#

# Startup
cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__)

using DelimitedFiles
using Printf
using CairoMakie
CairoMakie.activate!()

# Directories
dir_in = "C:/Users/lsalo/julia/output/wilmington/"
figs_dir = "C:/Users/lsalo/julia/graphics/figs/wilmington/test/" 

# Read initial stress data and convert to cap stress space
id_ranger = 1:67736     # elements of Ranger elementset, grid file WilmingtonMesh_AllFaults_Rotated_Coarse.inp
stress_ini = readdlm(join([dir_in, "initial_stress_rotated_criticallyStressed_mu0.6.inp"]), ',')
stress_ini = stress_ini[id_ranger, 2:end]/1e6 # MPa
p_ini = vec(-(1/3)*sum(stress_ini[:,1:3], dims=2))
#stress_dev_ini = stress_ini
#stress_dev_ini[:, 1:3] = stress_dev_ini[:, 1:3] .+ p_ini
q_ini = sqrt.( 0.5*((stress_ini[:,1] - stress_ini[:,2]).^2 .+ (stress_ini[:,2] - stress_ini[:,3]).^2 .+ (stress_ini[:,3] - stress_ini[:,1]).^2) ) # ok for shear stress = 0

centroids = readdlm(join([dir_in, "element_centroids.dat"]), ',')
depth_ranger_elems = centroids[id_ranger, end]
nodal_depth = readdlm(join([dir_in, "nodal_depth.dat"]), ',')
peff_ini_nod = readdlm(join([dir_in, "peff_ini_nod.inp"]), ',')

# Cap model parameters 
d = 1e-3        # MPa
b = 52.41       # deg
R = 0.2         #1.227
W = 0.27        # [-]
D = 0.7143 # 1/kbar

# Hardening curve
# We define hydrostatic yield (p_b) as 110% of the initial p_eff
#p_b = maximum(p_ini) .+ ((1/3)*(1/D)*log.(e_p./W .+ 1)*1e2) # MPa
#z_pt = -[935, 940, 2835, 2840]
z_pt = -[890, 900, 1420, 1430, 2890, 3000]
nc = Int(length(z_pt)/2)
e_p = [0, 0.0025, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2]
p_b = Vector{}(undef, nc)
peff_id = Vector{}(undef, nc)
for i in 1:nc
    id = findfirst(x -> z_pt[2*i] <= x <= z_pt[2*i-1], nodal_depth[:, 2])
    peff_id[i] = peff_ini_nod[id, end]
    #p_b[i] = 1.2*peff_id[i]*exp.(e_p/W) 
    p_b[i] = 1.18*peff_id[i] .+ 0.25*((1/3)*(1/D)*log.(e_p./W .+ 1)*1e8) # Pa
end
p_b_str = map(x -> @sprintf("%.2f", x), vcat(p_b...))
peff_id_str = map(x -> @sprintf("%.2f", x), peff_id)
to_save = hcat(p_b_str, repeat(e_p, nc, 1), zeros(nc*length(e_p)), repeat(peff_id_str, inner=length(e_p))) # NOTE: number formats like 1.315e7 won't work in Abaqus
writedlm(join([dir_in, "hardening_curve_test_peff.inp"]), to_save, ",\t")

# Plot
# Cap Position
e_p0 = 0
p_Fs = Vector{}(0:10:100)
q_Fs = p_Fs.*tand(b) .+ d

p_a = (1e-6*p_b[end][1] - R*d) / (1 + R*tand(b))
#q_Fc = vcat(Vector{}(0:1:31), 31.1:.1:31.7, 31.71:.01:31.73)
q_Fc = Vector{}(0:1:78)
p_Fc = p_a .+ sqrt.( abs.((R*(d + p_a*tand(b)))^2 .- (R.*q_Fc).^2) )

# Figure
f = Figure(size=(800, 400))
ax = Axis(f[1, 1], xlabel = L"p^\text{eff} = -\frac{\sigma_{kk}^\text{eff}}{3} \; \text{ [MPa]}", ylabel = L"q = \sqrt{3 J_2} \; \text{ [MPa]}", xlabelsize = 16, ylabelsize = 16)
ln1 = lines!(ax, p_Fs, q_Fs, color = RGBAf(0, 0, 0, 1), linewidth = 1.5, label=L"F_\text{s}")
ln1 = lines!(ax, p_Fc, q_Fc, color = RGBAf(0.5, 0.5, 0.5, 1), linewidth = 1, linestyle = :dashdot, label=L"F_\text{c}")
sc1 = scatter!(ax, p_ini, q_ini, markersize = 5, color=RGBAf(0,0,0,0), strokecolor=:blue, strokewidth=0.5, label=L"\text{\sigma_0 (Ranger)}")
#text!(ax, p_b0+1, 0, text=L"p_\text{b}", fontsize=14)
xlims!(ax, 0, 100)
ylims!(ax, 0, 100)
axislegend(position = :lt, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8)) 

ax = Axis(f[1, 2], xlabel = L"e^\text{in}_\text{v} \; \text{ [-]}", ylabel = L"p_\text{b} \; \text{ [MPa]}", xlabelsize = 16, ylabelsize = 16)
ln1 = lines!(ax, e_p, p_b[1]*1e-6, color = RGBAf(0, 0, 0, 1), linewidth = 1.5)
ln1 = lines!(ax, e_p, p_b[2]*1e-6, color = RGBAf(0, 0, 0, 1), linewidth = 1.5)
ln1 = lines!(ax, e_p, p_b[3]*1e-6, color = RGBAf(0, 0, 0, 1), linewidth = 1.5)
xlims!(ax, 0, 0.2)
ylims!(ax, 0, 100)
f

#save(join([figs_dir, "cap_model.pdf"]), f, pt_per_unit = 1)