#=
Script to plot modified Drucker-Prager cap model and initial reservoir state.
=#

# Startup
cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__)

using LaTeXStrings
using DelimitedFiles
using Printf
using CairoMakie
CairoMakie.activate!()

# Directories and options
dir_in = "C:/Users/lsalo/julia/output/wilmington/"
figs_dir = "C:/Users/lsalo/julia/graphics/figs/wilmington/test/" 

cap_case = "ranger"  # entire_domain | ranger | ranger-pico | ranger-pico-overPico
model_name = "cam-clay"   # mdpcm | cam-clay (isotropic, K=1)


# Read initial stress data and convert to cap stress space
centroids       = readdlm(join([dir_in, "element_centroids.dat"]), ',')
nodal_depth     = readdlm(join([dir_in, "nodal_depth.dat"]), ',')
peff_ini_nod    = readdlm(join([dir_in, "peff_ini_nod.inp"]), ',')
stress_ini      = readdlm(join([dir_in, "initial_stress_rotated_criticallyStressed_mu0.6.inp"]), ',')
if cap_case == "ranger"
    id_ranger = 1:67736     # elements of Ranger elementset, grid file WilmingtonMesh_AllFaults_Rotated_Coarse.inp
    stress_ini = stress_ini[id_ranger, 2:end]/1e6   # MPa
    depth_ranger_elems = centroids[id_ranger, end]
elseif cap_case == "entire_domain"
    stress_ini = stress_ini[:, 2:end]/1e6           # MPa
elseif cap_case == "ranger-pico"
    id_rp = vec(readdlm(join([dir_in, "id_ranger_pico.dat"]), Int))         # grid file WilmingtonMesh_AllFaults_Rotated_Coarse.inp
    stress_ini = stress_ini[id_rp, 2:end]/1e6   # MPa
    depth_elems = centroids[id_rp, end]
elseif cap_case == "ranger-pico-overPico"
    id_rpo = vec(readdlm(join([dir_in, "id_ranger_pico_overPico.dat"]), Int))  # grid file WilmingtonMesh_AllFaults_Rotated_Coarse.inp
    stress_ini = stress_ini[id_rpo, 2:end]/1e6   # MPa
    depth_elems = centroids[id_rpo, end]
elseif cap_case == "above-basement"
    id_ab = vec(readdlm(join([dir_in, "id_all_above_basement.dat"]), Int))  # grid file WilmingtonMesh_AllFaults_Rotated_Coarse.inp
    stress_ini = stress_ini[id_ab, 2:end]/1e6   # MPa
    depth_elems = centroids[id_ab, end]
end
p_ini = vec(-(1/3)*sum(stress_ini[:,1:3], dims=2))
q_ini = sqrt.( 0.5*((stress_ini[:,1] - stress_ini[:,2]).^2 .+ (stress_ini[:,2] - stress_ini[:,3]).^2 .+ (stress_ini[:,3] - stress_ini[:,1]).^2) ) # ok for shear stress = 0


# Define Cap model 
if model_name == "mdpcm"
    # Cap model parameters 
    d = 1e-3        # MPa
    b = 52.41       # deg
    R = 0.2         #1.227
    W = 0.27        # [-]
    D = 0.7143      # 1/kbar

    # Hardening curve
    # We define hydrostatic yield (p_b) as 118% of the initial p_eff
    if cap_case == "ranger"
        z_pt = -[890, 900, 1420, 1430, 2890, 3000]
        e_p = [0, 0.0025, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2]
        fn_hardening = "hardening_curve_peff_0.1Kosloff_1.15.inp"
    elseif cap_case == "entire_domain"
        z_pt = -[0, 1, 990, 1000, 2890, 3000, 5890, 6000, 14995, 15005]
        e_p = [0, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2]
        fn_hardening = "hardening_curve_peff_all.inp"
    elseif cap_case == "ranger-pico"
        z_pt = -[400, 405, 890, 900, 1420, 1430, 2890, 3000]
        e_p = [0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]
        fn_hardening = "hardening_curve_peff_0.1Kosloff_1.15_ranger-pico.inp"
    elseif cap_case == "ranger-pico-overPico"
        z_pt = -[0, 1, 400, 405, 890, 900, 1420, 1430, 2890, 3000]
        e_p = [0, 0.0025, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2]
        fn_hardening = "hardening_curve_peff_0.25Kosloff_ranger-pico-overPico.inp"
    elseif cap_case == "above-basement"
        z_pt = -[0, 1, 890, 900, 1420, 1430, 2890, 3000, 4890, 5000, 8890, 9000]
        e_p = [0, 0.0025, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.2]
        fn_hardening = "hardening_curve_peff_0.25Kosloff_above_basement.inp"
    end
    nc = Int(length(z_pt)/2)
    p_b = Vector{}(undef, nc)
    peff_id = Vector{}(undef, nc)
    for i in 1:nc
        id = findfirst(x -> z_pt[2*i] <= x <= z_pt[2*i-1], nodal_depth[:, 2])
        peff_id[i] = peff_ini_nod[id, end]
        p_b[i] = 1.15*peff_id[i] .+ 0.1*((1/3)*(1/D)*log.(e_p./W .+ 1)*1e8) # Pa
    end

    # write to file
    p_b_str = map(x -> @sprintf("%.2f", x), vcat(p_b...))
    peff_id_str = map(x -> @sprintf("%.2f", x), peff_id)
    to_save = hcat(p_b_str, repeat(e_p, nc, 1), zeros(nc*length(e_p)), repeat(peff_id_str, inner=length(e_p))) # NOTE: number formats like 1.315e7 won't work in Abaqus
    writedlm(join([dir_in, fn_hardening]), to_save, ",\t")

elseif model_name == "cam-clay"
    # *Clay plasticity parameters
    l = 0.1
    phi = 32 # friction angle
    M = round(6*sind(phi)/(3 - sind(phi)), digits=2)
    b = 0.5
    K = 1.0
    p_t = 0 # no cohesion, must be 0 for exponential hardening as defined here

    # hardening parameter (a0)
    if cap_case == "ranger"
        z_pt = -[890, 900, 1420, 1430, 2890, 3000]
        fn_camclay = "clay_plasticity_pc_1.4peff.inp"
    end

    nc = Int(length(z_pt)/2)
    p_c = Vector{}(undef, nc)       # yield stress in hydrostatic compression
    a0 = Vector{}(undef, nc) 
    peff_id = Vector{}(undef, nc)
    for i in 1:nc
        id = findfirst(x -> z_pt[2*i] <= x <= z_pt[2*i-1], nodal_depth[:, 2])
        peff_id[i] = peff_ini_nod[id, end]
        p_c[i] = 1.4*peff_id[i]    # Pa
        a0[i] = p_c[i]/(1+b)
    end

    # write to file
    a0_str = map(x -> @sprintf("%.2f", x), a0)
    peff_id_str = map(x -> @sprintf("%.2f", x), peff_id)
    to_save = hcat( fill(l, nc), fill(M, nc), a0_str, fill(b, nc), fill(K, nc), zeros(nc, 1), peff_id_str)
    writedlm(join([dir_in, fn_camclay]), to_save, ",\t")
end


# Plot
# Cap Position
if model_name == "mdpcm"
    e_p0 = 0
    if cap_case == "ranger" || cap_case == "ranger-pico"
        p_Fs = Vector{}(0:10:100)
        q_Fc = vcat(Vector{}(0:1:79), Vector{}(79.1:.1:79.5))
    elseif cap_case == "entire_domain"
        p_Fs = Vector{}(0:10:400)
        q_Fc = vcat(Vector{}(0:10:400), Vector{}(400:1:407))
    end
    q_Fs = p_Fs.*tand(b) .+ d
    p_a = (1e-6*p_b[end][1] - R*d) / (1 + R*tand(b))
    p_Fc = p_a .+ sqrt.( abs.((R*(d + p_a*tand(b)))^2 .- (R.*q_Fc).^2) )

    f = Figure(size=(800, 400))
    ax = Axis(f[1, 1], xlabel = L"p^\text{eff} = -\frac{\sigma_{kk}^\text{eff}}{3} \; \text{ [MPa]}", ylabel = L"q = \sqrt{3 J_2} \; \text{ [MPa]}", xlabelsize = 16, ylabelsize = 16)
    ln1 = lines!(ax, p_Fs, q_Fs, color = RGBAf(0, 0, 0, 1), linewidth = 2, label=L"F_\text{s}")
    ln1 = lines!(ax, p_Fc, q_Fc, color = RGBAf(0.5, 0.5, 0.5, 1), linewidth = 1.5, linestyle = :solid, label=L"F_\text{c}")
    sc1 = scatter!(ax, p_ini, q_ini, markersize = 3, color=RGBAf(0,0,0,0), strokecolor=:blue, strokewidth=0.5, label=L"\sigma_0 \; \text{ (elements)}")
    #text!(ax, p_b0+1, 0, text=L"p_\text{b}", fontsize=14)
    if cap_case == "ranger" || cap_case == "ranger-pico"
        xlims!(ax, 0, 100)
        ylims!(ax, 0, 100)
    elseif cap_case == "entire_domain"
        xlims!(ax, 0, 400)
        ylims!(ax, 0, 450)
    end
    axislegend(position = :lt, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8)) 

    ax = Axis(f[1, 2], xlabel = L"e^\text{in}_\text{v} \; \text{ [-]}", ylabel = L"p_\text{b} \; \text{ [MPa]}", xlabelsize = 16, ylabelsize = 16)
    ln1 = lines!(ax, e_p, p_b[1]*1e-6, color = RGBAf(0, 0, 0, 1), linewidth = 1.5)
    ln1 = lines!(ax, e_p, p_b[2]*1e-6, color = RGBAf(0, 0, 0, 1), linewidth = 1.5)
    ln1 = lines!(ax, e_p, p_b[3]*1e-6, color = RGBAf(0, 0, 0, 1), linewidth = 1.5)
    xlims!(ax, 0, 0.2)
    ylims!(ax, 0, 100)
    f

elseif model_name == "cam-clay"
    if cap_case == "ranger" || cap_case == "ranger-pico"
        p_Fs = Vector{}(0:10:100)
        #q_Fc = vcat(Vector{}(0:1:25), Vector{}(25.1:.1:25.2))
        #q_Fc = vcat(Vector{}(0:1:76), Vector{}(76.1:.1:76.1))
        q_Fc = vcat(Vector{}(0:1:78), Vector{}(78.1:.1:78.6))
    end

    q_Fs = M*p_t .+ M.*p_Fs
    a = (p_c[end] - p_t) / (1 + b)
    p_Fc = (p_t .+ a.*( 1 .+ sqrt.( b^2 .*(1 .- (q_Fc.*1e6./(M*a)).^2) ) )) ./ 1e6  # MPa

    f = Figure(size=(800, 400))
    ax = Axis(f[1, 1], xlabel = L"p^\text{eff} = -\frac{\sigma_{kk}^\text{eff}}{3} \; \text{ [MPa]}", ylabel = L"q = \sqrt{3 J_2} \; \text{ [MPa]}", xlabelsize = 16, ylabelsize = 16, xticks = 0:25:100, yticks=0:25:100)
    ln1 = lines!(ax, p_Fs, q_Fs, color = RGBAf(0, 0, 0, 1), linewidth = 2, label=L"F_\text{s}")
    ln1 = lines!(ax, p_Fc, q_Fc, color = RGBAf(0.5, 0.5, 0.5, 1), linewidth = 1.5, linestyle = :solid, label=L"F_\text{c}")
    sc1 = scatter!(ax, p_ini, q_ini, markersize = 3, color=RGBAf(0,0,0,0), strokecolor=:blue, strokewidth=0.5, label=L"\sigma_0 \; \text{ (elements)}")
    #text!(ax, p_b0+1, 0, text=L"p_\text{b}", fontsize=14)
    if cap_case == "ranger" || cap_case == "ranger-pico"
        xlims!(ax, 0, 100)
        ylims!(ax, 0, 100)
    end
    axislegend(position = :lt, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8)) 
    f

end

# Figure


#save(join([figs_dir, "cap_model.pdf"]), f, pt_per_unit = 1)