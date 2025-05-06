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
using Interpolations
using CairoMakie
CairoMakie.activate!()

# Directories and options
dir_in = "C:/Users/lsalo/julia/output/wilmington/"
figs_dir = "C:/Users/lsalo/julia/graphics/figs/wilmington/paper1/" 

regime = "strike-slip"
a_phi = "1.75"
cap_case = "above-basement"  # entire_domain | above-basement | ranger | ranger-pico | ranger-pico-overPico | pico | overPico | bottomRanger
model_name = "mdpcm"   # mdpcm | cam-clay (isotropic, K=1)


# Read initial stress data and convert to cap stress space
centroids       = readdlm(join([dir_in, "element_centroids.dat"]), ',')
nodal_depth     = readdlm(join([dir_in, "nodal_depth.dat"]), ',')
if regime == "reverse"
    #peff_ini_nod    = readdlm(join([dir_in, "peff_ini_nod.inp"]), ',')
    #stress_ini      = readdlm(join([dir_in, "initial_stress_rotated_criticallyStressed_mu0.6.inp"]), ',')
    peff_ini_nod    = readdlm(join([dir_in, "peff_ini_nod_aphi2.5.inp"]), ',')
    stress_ini      = readdlm(join([dir_in, "initial_stress_rotated_criticallyStressed_mu0.6_aphi2.5.inp"]), ',')
elseif regime == "normal"
    peff_ini_nod    = readdlm(join([dir_in, "peff_ini_nod_normal.inp"]), ',')
    stress_ini      = readdlm(join([dir_in, "initial_stress_rotated_criticallyStressed_mu0.6_normal.inp"]), ',')
elseif regime == "strike-slip"
    #peff_ini_nod    = readdlm(join([dir_in, "peff_ini_nod_aphi1.5.inp"]), ',')
    #stress_ini      = readdlm(join([dir_in, "initial_stress_rotated_criticallyStressed_mu0.6_aphi1.5.inp"]), ',')
    peff_ini_nod    = readdlm(join([dir_in, "peff_ini_nod_aphi1.75_sed_mu0.3_bas_mu0.6.inp"]), ',')
    stress_ini      = readdlm(join([dir_in, "initial_stress_rotated_sed_mu0.3_bas_mu0.6_aphi1.75.inp"]), ',')
    #peff_ini_nod    = readdlm(join([dir_in, "peff_ini_nod_sed_mu0.3_aphi1.25_bas_mu0.6_aphi1.75.inp"]), ',')
    #stress_ini      = readdlm(join([dir_in, "initial_stress_rotated_sed_mu0.3_aphi1.25_bas_mu0.6_aphi1.75.inp"]), ',')
end

if cap_case == "ranger"
    id_cells = 1:67736     # elements of Ranger elementset, grid file WilmingtonMesh_AllFaults_Rotated_Coarse.inp
elseif cap_case == "entire_domain"
    id_cells = 1:length(centroids[:, 1])
elseif cap_case == "ranger-pico"
    id_cells = vec(readdlm(join([dir_in, "id_ranger_pico.dat"]), Int))              
elseif cap_case == "ranger-pico-overPico"
    id_cells = vec(readdlm(join([dir_in, "id_ranger_pico_overPico.dat"]), Int))    
elseif cap_case == "above-basement"
    id_cells = vec(readdlm(join([dir_in, "id_all_above_basement.dat"]), Int))       
elseif cap_case == "pico"
    id_cells = vec(readdlm(join([dir_in, "id_pico.dat"]), Int))  
elseif cap_case == "overPico"
    id_cells = vec(readdlm(join([dir_in, "id_overPico.dat"]), Int))  
elseif cap_case == "bottomRanger"
    id_cells = vec(readdlm(join([dir_in, "id_bottomRanger.dat"]), Int))  
end
stress_ini = stress_ini[id_cells, 2:end]/1e6   # MPa
depth_elems = centroids[id_cells, end]
p_ini = vec(-(1/3)*sum(stress_ini[:,1:3], dims=2))
q_ini = sqrt.( 0.5*((stress_ini[:,1] - stress_ini[:,2]).^2 .+ (stress_ini[:,2] - stress_ini[:,3]).^2 .+ (stress_ini[:,3] - stress_ini[:,1]).^2) ) # ok for shear stress = 0


# Define Cap model 
if model_name == "mdpcm"
    # Cap model parameters 
    d = 1e-3        # MPa
    b = 55.0        # 58.6 (\phi = 40) 52.41 deg (32), 55 (35)
    R = 2.0         # 0.4
    W = 0.27        # [-]
    D = 0.7143      # 1/kbar

    # Hardening curve
    # We define hydrostatic yield (p_b) as 118% of the initial p_eff
    if cap_case == "ranger"
        z_pt = -[930, 940, 1420, 1430, 2830, 2840]
        #z_pt = -[890, 900, 1420, 1430, 2890, 3000]
        #e_p = [0, 0.0025, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2]
        e_p = [0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]
        f_p = 1.7 #1.28 #1.12 #1.08 #1.01 (R=0.001) 1.02 (R=0.02) 1.08 (R=0.1) 1.15-1.18 (R=0.2) 1.41 (R=0.5) 1.35 (R=0.3, normal) 1.2 (R=0.3, strike-slip)
        f_K = 0.5
    elseif cap_case == "entire_domain"
        z_pt = -[0, 1, 990, 1000, 2890, 3000, 5890, 6000, 14995, 15005]
        e_p = [0, 0.0025, 0.005, 0.01, 0.025, 0.05, 0.1, 0.2]
        f_p = 1.18
        f_K = 0.25
    elseif cap_case == "ranger-pico"
        z_pt = -[400, 405, 890, 900, 1420, 1430, 2890, 3000]
        e_p = [0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]
        f_p = 1.18
        f_K = 0.25
    elseif cap_case == "ranger-pico-overPico"
        z_pt = -[0, 1, 400, 405, 890, 900, 1420, 1430, 2890, 3000]
        e_p = [0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]
        f_p = 1.15
        f_K = 0.25
    elseif cap_case == "above-basement"
        z_pt = -[5, 6, 890, 900, 1420, 1430, 2890, 3000, 4890, 5000, 8890, 9000]
        e_p = [0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]
        f_p = 1.7 #1.27 #1.17 #1.05 #1.12 #1.2 1.35 #1.02 #1.18 #1.41 #1.08 #1.15
        f_K = 0.5
    elseif cap_case == "pico"
        z_pt = -[475, 485, 950, 1000, 2620, 2650]
        #z_pt = -[450, 460, 950, 1000, 2620, 2650]
        e_p = [0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]
        f_p = 1.17
        f_K = 0.5
    elseif cap_case == "overPico"
        #z_pt = -[0, 1, 475, 485, 1360, 1370]
        z_pt = -[25, 30, 475, 485, 1355, 1360]
        e_p = [0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]
        f_p = 1.05 #1.41 #1.08 #1.15
        f_K = 0.5
    elseif cap_case == "bottomRanger"
        z_pt = -[1160, 1170, 1990, 2000, 4790, 5000]
        #z_pt = -[25, 30, 475, 485, 1360, 1370]
        e_p = [0, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]
        f_p = 1.05
        f_K = 0.5
    end
    fn_hardening = "hardening_curve_" * string(f_p) * "peff_" * string(R) * "R_" * string(f_K) * "Kosloff_" * cap_case * "_aphi" * a_phi * ".inp"
    nc = Int(length(z_pt)/2)
    p_b = Vector{}(undef, nc)
    peff_id = Vector{}(undef, nc)
    for i in 1:nc
        id = findfirst(x -> z_pt[2*i] <= x <= z_pt[2*i-1], nodal_depth[:, 2])
        if isnothing(id)
            #peff_id[i] = 0.5*minimum(peff_ini_nod[peff_ini_nod[:, end].>0, end])
            print("NOTE: no id found for depth ", z_pt[2*i], " to ", z_pt[2*i-1], ": interpolation used. \n")
            id_z = unique(i -> nodal_depth[i, 2], eachindex(nodal_depth[:, 2]))   # Get first indices
            z_unique = nodal_depth[id_z, 2] 
            peff_unique = peff_ini_nod[id_z, 2]    
            id_z_sorted = sortperm(z_unique)
            itp = interpolate((z_unique[id_z_sorted],), peff_unique[id_z_sorted], Gridded(Linear()))
            peff_id[i] = itp(z_pt[2*i-1])
        else
            peff_id[i] = peff_ini_nod[id, end]
        end
        p_b[i] = f_p*peff_id[i] .+ f_K*((1/3)*(1/D)*log.(e_p./W .+ 1)*1e8) # Pa
    end

    # write to file
    p_b_str = map(x -> @sprintf("%.2f", x), vcat(p_b...))
    peff_id_str = map(x -> @sprintf("%.2f", x), peff_id)
    to_save = hcat(p_b_str, repeat(e_p, nc, 1), zeros(nc*length(e_p)), repeat(peff_id_str, inner=length(e_p))) # NOTE: number formats like 1.315e7 won't work in Abaqus
    writedlm(join([dir_in, fn_hardening]), to_save, ",\t")

elseif model_name == "cam-clay"
    # *Clay plasticity parameters
    l = 0.2
    phi = 32 # friction angle
    M = round(6*sind(phi)/(3 - sind(phi)), digits=2)
    b = 0.25
    K = 1.0
    p_t = 0 # no cohesion, must be 0 for exponential hardening as defined here

    # hardening parameter (a0)
    if cap_case == "ranger"
        z_pt = -[890, 900, 1420, 1430, 2890, 3000]
        #z_pt = -[930, 935, 1420, 1430, 2835, 2840]
        fn_camclay = "clay_plasticity_pc1.19peff_b0.25.inp"
    end

    nc = Int(length(z_pt)/2)
    p_c = Vector{}(undef, nc)       # yield stress in hydrostatic compression
    a0 = Vector{}(undef, nc) 
    peff_id = Vector{}(undef, nc)
    for i in 1:nc
        id = findfirst(x -> z_pt[2*i] <= x <= z_pt[2*i-1], nodal_depth[:, 2])
        peff_id[i] = peff_ini_nod[id, end]
        p_c[i] = 1.19*peff_id[i]    # Pa
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
        #q_Fc = vcat(Vector{}(0:1:79), Vector{}(79.1:.1:79.5))   # reverse regime pb[end]
        #q_Fc = vcat(Vector{}(0:1:26), Vector{}(26.1:.1:26.1))  # reverse regime pb[1]
        #q_Fc = vcat(Vector{}(0:1:25), Vector{}(25.02:.02:26.66))   # normal
        #q_Fc = vcat(Vector{}(0:1:8), Vector{}(8.05:.05:9.15)) 
        #q_Fc = vcat(Vector{}(0:1:24), Vector{}(24:.01:25.3)) #strike-slip
        q_Fc = vcat(Vector{}(0:1:7), Vector{}(7.01:.01:8.7)) #strike-slip
    elseif cap_case == "entire_domain"
        p_Fs = Vector{}(0:10:400)
        q_Fc = vcat(Vector{}(0:10:400), Vector{}(400:1:407))
    elseif cap_case == "pico"
        p_Fs = Vector{}(0:10:100)
        #q_Fc = vcat(Vector{}(0:1:71), Vector{}(71.1:.1:71.2))
        #q_Fc = vcat(Vector{}(0:1:15), Vector{}(15.1:.1:15.5))
        q_Fc = vcat(Vector{}(0:1:28), Vector{}(28.01:.01:28.4))
        #q_Fc = vcat(Vector{}(0:.2:3), Vector{}(3.01:.01:4.4))
    elseif cap_case == "overPico"
        p_Fs = Vector{}(0:10:100)
        #q_Fc = vcat(Vector{}(0:1:37), Vector{}(37.1:.1:37.5))
        #q_Fc = vcat(Vector{}(0:.1:3), Vector{}(3.01:.01:3.1))
        #q_Fc = vcat(Vector{}(0:1:13), Vector{}(13.1:.1:13.7))
        q_Fc = Vector{}(0:.01:0.46)
    elseif cap_case == "bottomRanger"
        p_Fs = Vector{}(0:10:150)
        #q_Fc = vcat(Vector{}(0:1:80), Vector{}(80.1:.1:80.5))
        q_Fc = vcat(Vector{}(0:1:31), Vector{}(31.1:.1:31.3))
    end
    q_Fs = p_Fs.*tand(b) .+ d
    #p_a = (1e-6*p_b[end][1] - R*d) / (1 + R*tand(b))
    p_a = (1e-6*p_b[1][1] - R*d) / (1 + R*tand(b))
    p_Fc = p_a .+ sqrt.( abs.((R*(d + p_a*tand(b)))^2 .- (R.*q_Fc).^2) )

    f = Figure(size=(800, 400))
    #ax = Axis(f[1, 1], xlabel = L"p^\text{eff} = -\frac{\sigma_{kk}^\text{eff}}{3} \; \text{ [MPa]}", ylabel = L"q = \sqrt{3 J_2} \; \text{ [MPa]}", xlabelsize = 16, ylabelsize = 16, xticks = 0:25:125, yticks=0:25:125)
    ax = Axis(f[1, 1], xlabel = L"p^\text{eff} = -\frac{\sigma_{kk}^\text{eff}}{3} \; \text{ [MPa]}", ylabel = L"q = \sqrt{3 J_2} \; \text{ [MPa]}", xlabelsize = 16, ylabelsize = 16)
    ln1 = lines!(ax, p_Fs, q_Fs, color = RGBAf(0, 0, 0, 1), linewidth = 2, label=L"F_\text{s}")
    ln1 = lines!(ax, p_Fc, q_Fc, color = RGBAf(0.5, 0.5, 0.5, 1), linewidth = 1.5, linestyle = :solid, label=L"F_\text{c}")
    sc1 = scatter!(ax, p_ini, q_ini, markersize = 3, color=RGBAf(0,0,0,0), strokecolor=:blue, strokewidth=0.5, label=L"\sigma_0 \; \text{ (elements)}")
    #text!(ax, p_b0+1, 0, text=L"p_\text{b}", fontsize=14)
    if cap_case == "ranger" || cap_case == "ranger-pico"
        xlims!(ax, 0, 25)   # 30
        ylims!(ax, 0, 25)   # 30
    elseif cap_case == "entire_domain"
        xlims!(ax, 0, 400)
        ylims!(ax, 0, 450)
    elseif cap_case == "pico"
        xlims!(ax, 0, 40)
        ylims!(ax, 0, 40)
    elseif cap_case == "overPico"
        xlims!(ax, 0, 5)
        ylims!(ax, 0, 5)
    elseif cap_case == "bottomRanger"
        #xlims!(ax, 0, 80)
        #ylims!(ax, 0, 80)
        xlims!(ax, 0, 30)
        ylims!(ax, 0, 35)
    end
    axislegend(position = :lt, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8)) 

    ax = Axis(f[1, 2], xlabel = L"e^\text{in}_\text{v} \; \text{ [-]}", ylabel = L"p_\text{b} \; \text{ [MPa]}", xlabelsize = 16, ylabelsize = 16)
    ln1 = lines!(ax, e_p, p_b[1]*1e-6, color = RGBAf(0, 0, 0, 1), linewidth = 1.5)
    ln1 = lines!(ax, e_p, p_b[2]*1e-6, color = RGBAf(0, 0, 0, 1), linewidth = 1.5)
    ln1 = lines!(ax, e_p, p_b[3]*1e-6, color = RGBAf(0, 0, 0, 1), linewidth = 1.5)
    xlims!(ax, 0, 0.2)
    ylims!(ax, 0, 60)
    f

elseif model_name == "cam-clay"
    if cap_case == "ranger" || cap_case == "ranger-pico"
        p_Fs = Vector{}(0:10:125)
        #q_Fc = vcat(Vector{}(0:1:25), Vector{}(25.1:.1:25.2))
        #q_Fc = vcat(Vector{}(0:1:76), Vector{}(76.1:.1:76.9))
        #q_Fc = vcat(Vector{}(0:1:78), Vector{}(78.1:.1:78.6))
        
        #q_Fc = vcat(Vector{}(0:1:77), Vector{}(77.1:.1:77.9))
        #q_Fc = vcat(Vector{}(0:1:23), Vector{}(23.1:.1:23.8))

        # b=0.25, pc = 1.19peff:
        q_Fc = vcat(Vector{}(0:1:80), Vector{}(80.1:.1:80.2))
        #q_Fc = vcat(Vector{}(0:1:25), Vector{}(25.1:.1:25.5))
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