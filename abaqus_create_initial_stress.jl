#=
Script to initialize stresses at Wilmington using FlatMaxwell model from SCEC CSM.

Important notes:
    1. We removed character # from CSM stress .csv files, row 48, where the header is.
    2. In our convention, S1 = maximum compressive, S2 = intermediate and S3 = minimum compressive or most tensional. Note that this is different than the convention used in the CSM stress .csv files
    3. In Abaqus, Wilmington grid, s11 is // to x, s22 // y, and s33 // z.
=#

# Cd into current directory and activate the Julia environment in this directory
cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__)

# Startup
using Revise
using CSV
using CairoMakie
CairoMakie.activate!()
using LaTeXStrings
using GLM 
using Statistics
using DelimitedFiles
using Printf
#using Colors
#using ColorSchemes
using Distances

lib_dir = "C:/Users/lsalo/julia/jl-tools/"
includet(joinpath(lib_dir, "modules", "abaqus.jl"))
using .Abaqus

# Directories
dir_inp = "C:/Users/lsalo/Wilmington/AbaqusGrids/"
dir_data = "C:/Users/lsalo/Wilmington/2023_SCEC_CSM/SCEC_CSM_archive_2023/"
dir_out = "C:/Users/lsalo/julia/output/wilmington/"
figs_dir = "C:/Users/lsalo/julia/graphics/figs/wilmington/paper1/" 


# 1. STRESS MODELS DATA
filenames = ["FlatMaxwell", "Hardebeck_FM", "Luttrell-2017", "YH14-K", "YHSM-2013"]
stress_data_all = Vector{}(undef, length(filenames))
for (i, fn) in enumerate(filenames)
    if fn == "FlatMaxwell"
        stress_data_all[i] = CSV.File(join([dir_data, fn, ".csv"]), header=48, select=1:20, skipto=49)
    else
        stress_data_all[i] = CSV.File(join([dir_data, fn, ".csv"]), header=48, select=[1, 2, 3, 14], skipto=49)
    end
end

# Get subset: Wilmington area
min_lon, max_lon = -118.41947, -117.93033 
min_lat, max_lat = 33.58679, 33.93842
stress_data = Vector{}(undef, length(filenames))
for i in 1:length(filenames)
    row_ids = [min_lon <= row.LON <= max_lon && min_lat <= row.LAT <= max_lat && 0 <= row.DEP <= 16 for row in stress_data_all[i]]
    stress_data[i] = Matrix{Float64}(undef, sum(row_ids), length(stress_data_all[i].names))
    k = 1
    for (j, row) in enumerate(stress_data_all[i])
        if row_ids[j] == 1
            stress_data[i][k, :] = [getproperty(row, name) for name in stress_data_all[i].names]
            k += 1
        end
    end
end

# FlatMaxwell least squares linear fit through the origin
x = stress_data[1][:,3]
lm_s1 = glm(@formula(y ~ x + x2 + x3 + x4 + 0), (x=x, x2=x.^2, x3=x.^3, x4=x.^4, y=stress_data[1][:,20]), Normal(), IdentityLink())
lm_s2 = glm(@formula(y ~ x + x2 + x3 + x4 + 0), (x=x, x2=x.^2, x3=x.^3, x4=x.^4, y=stress_data[1][:,19]), Normal(), IdentityLink())
lm_s3 = glm(@formula(y ~ x + 0), (x=x, y=stress_data[1][:,18]), Normal(), IdentityLink())
#lm_sv = glm(@formula(y ~ x + 0), (x=x, y=stress_data[1][:,9]), Normal(), IdentityLink()) # basically equal to s3

# A_phi fit
x = vcat(stress_data[1][:, 3], stress_data[2][:, 3], stress_data[3][:, 3], stress_data[4][:, 3], stress_data[5][:, 3])
y = vcat(stress_data[1][:, 14], stress_data[2][:, 4], stress_data[3][:, 4], stress_data[4][:, 4], stress_data[5][:, 4])
lm_aphi = glm(@formula(y ~ x), (x=x, y=y), Normal(), IdentityLink())
x_fmx = stress_data[1][:, 3]
lm_aphi_fmx = glm(@formula(y ~ x + x2 + x3 + x4), (x=x_fmx, x2 = x_fmx.^2, x3 = x_fmx.^3, x4 = x_fmx.^4, y=stress_data[1][:, 14]), Normal(), IdentityLink())


# 2. PREDICT STRESS MAGNITUDES BASED ON CRITICALLY STRESSED CRUST
# Ref: Lund Snee Thesis, Chapter 7
# NOTE: pore pressures are negative for plots, but positive otherwise, hence the change of signs when adding pore pressure.
# NOTE: Lund Snee Thesis uses positive compression, we use positive tension
regime = "a_fit"
z_plot = Vector{}(LinRange(0, -15, 41))
p_plot = 1000*9.807.*z_plot*1e3./1e6                                    # [NEGATIVE] MPa
sv_plot = 2300*9.807.*z_plot*1e3./1e6
mu      = 0.6
ratio_s1s3  = ((mu^2 +1)^0.5 + mu)^2                                    # effective stress ratio (S1 - p) / (S3 - p)
if regime == "a_fit"
    a_phi       = predict(lm_aphi, (x=(-1*z_plot),))
    s3_plot     = sv_plot                                               # MPa
    s1_plot     = ratio_s1s3*(s3_plot .- p_plot) .+ p_plot              # Reverse: S_H
    s2_plot     = (a_phi .- 2.0).*(s1_plot .- s3_plot) .+ s3_plot       # Reverse: S_h
elseif regime == "normal"
    a_phi       = 0
    s1_plot     = sv_plot
    s3_plot     = (s1_plot .- p_plot)/ratio_s1s3 .+ p_plot
    s2_plot     = a_phi*(s1_plot .- s3_plot) .+ s3_plot    # a_phi = 0 --> s2 = s3
elseif regime == "strike-slip"
    a_phi       = 1.5
    s2_plot     = sv_plot
    #s1_plot     = (1/(2-a_phi)*s2_plot .- 1/(2-a_phi)*p_plot .+ p_plot .+ p_plot.*(1/(ratio_s1s3*(2-a_phi)) -1/ratio_s1s3)) ./ (1 + 1/(ratio_s1s3*(2-a_phi)) -1/ratio_s1s3)
    A           = a_phi + (-a_phi +1)/3.12 -2
    s1_plot     = ((A + 1).*p_plot - s2_plot) ./ A
    s3_plot     = (s1_plot .- p_plot) ./ ratio_s1s3 .+ p_plot 
elseif regime == "reverse"
    a_phi       = 2.5
    s3_plot     = sv_plot                                               # MPa
    s1_plot     = ratio_s1s3*(s3_plot .- p_plot) .+ p_plot              
    s2_plot     = (a_phi .- 2.0).*(s1_plot .- s3_plot) .+ s3_plot       
else
    error("Regime not recognized")    
end


# 3. ELEMENT CENTROIDS AND ASSIGN STRESS (HOMOGENEOUS STRESS STATE)
fn_inp  = "WilmingtonMesh_AllFaults_Rotated_Coarse.inp"
nodes_inp, elements_inp, elset_names_inp, nset_inp, nset_names_inp = get_connectivity_from_inp(dir_inp, fn_inp, "tetrahedral")
elements_all = vcat(elements_inp[1:end-1]...)
centroids = get_element_centroids_inp(elements_all, nodes_inp)

# Assign stress based on depth: tension postivie
# ABAQUS coordinate system: s11 // x, s22 // y, s33 // z
z = centroids[:, end]
p = -1000*9.807.*z                  # NOTE: Positive pore pressure here (Positive tension in Abaqus)  
sv = 2300*9.807.*z  
t = zeros(length(z))
z_nod = nodes_inp[:, end]
p_nod = -1000*9.807.*z_nod  
sv_nod = 2300*9.807.*z_nod  
if regime == "a_fit"
    # Elements (initial stresses)
    a_phi   = predict(lm_aphi, (x=(-1e-3*z),))  
    s3      = sv  
    s1      = ratio_s1s3*(s3 .+ p) .- p 
    s2      = (a_phi .- 2.0).*(s1 .- s3) .+ s3 

    # Nodes (initial mean effective stress)
    a_phi = predict(lm_aphi, (x=(-1e-3*z_nod),)) 
    s3_nod = sv_nod 
    s1_nod = ratio_s1s3*(s3_nod .+ p_nod) .- p_nod 
    s2_nod = (a_phi .- 2.0).*(s1_nod .- s3_nod) .+ s3_nod

elseif regime == "normal"
    a_phi   = 0
    s1      = sv
    s3      = (s1 .+ p)/ratio_s1s3 .- p
    s2      = a_phi*(s1 .- s3) .+ s3   
    
    s1_nod  = sv_nod
    s3_nod  = (s1_nod .+ p_nod)/ratio_s1s3 .- p_nod
    s2_nod  = a_phi*(s1_nod .- s3_nod) .+ s3_nod 
    
elseif regime == "strike-slip"
    a_phi   = 1.5
    s2      = sv
    #s1      = (1/(2-a_phi)*s2 .+ 1/(2-a_phi)*p .- p .- p.*(1/(ratio_s1s3*(2-a_phi)) -1/ratio_s1s3)) ./ (1 + 1/(ratio_s1s3*(2-a_phi)) -1/ratio_s1s3)
    s1      = ((A + 1).*(-p) - s2) ./ A
    s3      = (s1 .+ p)/ratio_s1s3 .- p
    
    s2_nod  = sv_nod
    #s1_nod  = (1/(2-a_phi)*s2_nod .+ 1/(2-a_phi)*p_nod .- p_nod .- p_nod.*(1/(ratio_s1s3*(2-a_phi)) -1/ratio_s1s3)) ./ (1 + 1/(ratio_s1s3*(2-a_phi)) -1/ratio_s1s3)
    s1_nod = ((A + 1).*(-p_nod) - s2_nod) ./ A
    s3_nod  = (s1_nod .+ p_nod)/ratio_s1s3 .- p_nod

elseif regime == "reverse"
    a_phi   = 2.5
    s3      = sv
    s1      = ratio_s1s3*(s3 .+ p) .- p 
    s2      = (a_phi .- 2.0).*(s1 .- s3) .+ s3 

    s3_nod  = sv_nod
    s1_nod  = ratio_s1s3*(s3_nod .+ p_nod) .- p_nod 
    s2_nod  = (a_phi .- 2.0).*(s1_nod .- s3_nod) .+ s3_nod

else
    error("Regime not recognized")
end

# 4. ELEMENT CENTROIDS AND ASSIGN STRESS (HETEROG. STRESS STATE, varying mobilized \mu)
z = centroids[:, end]
t = zeros(length(z))
names_sed = ["Ranger", "BottomRanger", "Pico", "AquiferFootwall", "OverburdenPico", "OverburdenCenterSouth", "OverburdenCenterNorth", "OverburdenHangingWall", "OverburdenRight", "OverburdenLeft"]
id_sed = findall(x -> x in names_sed, elset_names_inp)
elements_sed = vcat(elements_inp[id_sed]...)
id_elements_sed = elements_sed[:, 1]
centroids_sed = centroids[id_elements_sed, :]
id_bas = [3, 7, 8, 9, 10] # Basement elements
elements_bas = vcat(elements_inp[id_bas]...)
id_elements_bas = elements_bas[:, 1]
centroids_bas = centroids[elements_bas[:, 1], :]

nodes_sed = sort(unique(vcat(elements_sed[:, 2:end]...)))
id_nodes_sed = findall(x -> x in nodes_sed, nodes_inp[:, 1])
nodes_sed = nodes_inp[id_nodes_sed, :]
nodes_bas = sort(unique(vcat(elements_bas[:, 2:end]...)))   
id_nodes_bas = findall(x -> x in nodes_bas, nodes_inp[:, 1])
id_nodes_bas = findall(x -> !(x in id_nodes_sed), id_nodes_bas) # nodes at basement top are already in sediment, so we remove them
nodes_bas = nodes_inp[id_nodes_bas, :]

# Stress conditions
a_phi = 1.75
mu = [0.3, 0.6]                         # sedimentary section, basement
ratio_s1s3 = ((mu.^2 .+1).^0.5 .+ mu).^2  
A  = a_phi .+ (-a_phi .+1)./ratio_s1s3 .-2

# Element-based values
p = -1000*9.807.*z                      # NOTE: Positive pore pressure here (Positive tension in Abaqus)  
sv = 2300*9.807.*z 
s2 = sv
s1 = Vector{Float64}(undef, length(centroids[:, 1]))
s1[id_elements_sed] = ((A[1] + 1).*(-p[id_elements_sed]) - s2[id_elements_sed]) ./ A[1]
s1[id_elements_bas] = ((A[2] + 1).*(-p[id_elements_bas]) - s2[id_elements_bas]) ./ A[2]
s3 = Vector{Float64}(undef, length(centroids[:, 1]))
s3[id_elements_sed] = (s1[id_elements_sed] .+ p[id_elements_sed])/ratio_s1s3[1] .- p[id_elements_sed]
s3[id_elements_bas] = (s1[id_elements_bas] .+ p[id_elements_bas])/ratio_s1s3[2] .- p[id_elements_bas]

# Node-based values
z_nod = nodes_inp[:, end]
p_nod = -1000*9.807.*z_nod  
sv_nod = 2300*9.807.*z_nod  
s2_nod = sv_nod
s1_nod = Vector{Float64}(undef, length(nodes_inp[:, 1]))
s1_nod[id_nodes_sed] = ((A[1] + 1).*(-p_nod[id_nodes_sed]) - s2_nod[id_nodes_sed]) ./ A[1]
s1_nod[id_nodes_bas] = ((A[2] + 1).*(-p_nod[id_nodes_bas]) - s2_nod[id_nodes_bas]) ./ A[2]
s3_nod = Vector{Float64}(undef, length(nodes_inp[:, 1]))
s3_nod[id_nodes_sed]  = (s1_nod[id_nodes_sed] .+ p_nod[id_nodes_sed])/ratio_s1s3[1] .- p_nod[id_nodes_sed]
s3_nod[id_nodes_bas]  = (s1_nod[id_nodes_bas] .+ p_nod[id_nodes_bas])/ratio_s1s3[2] .- p_nod[id_nodes_bas]

# SAVE FILE FOR USE WITH ABAQUS
# Also get mean stress, so that the hardening curve is a function of depth 
peff_ini_nod = vec(-(1/3)*sum(hcat(s1_nod .+ p_nod, s2_nod .+ p_nod, s3_nod .+ p_nod), dims=2))

# NOTE: EFFECTIVE STRESSES!
if regime == "a_fit" || regime == "reverse"        # reverse, s11 = sxx, s22 = syy, s33 = szz
    stress_out = hcat(centroids[:,1], s1 .+ p, s2 .+ p, s3 .+ p, t, t, t)
    writedlm(join([dir_out, "initial_stress_rotated_criticallyStressed_mu0.6_aphi", a_phi, ".inp"]), stress_out, ',')
    fnp = join(["peff_ini_nod_aphi", a_phi, ".inp"])
elseif regime == "normal"   # normal, s11 = szz, s22 = sxx, s33 = syy
    stress_out = hcat(centroids[:,1], s2 .+ p, s3 .+ p, s1 .+ p, t, t, t)
    writedlm(join([dir_out, "initial_stress_rotated_criticallyStressed_mu0.6_normal.inp"]), stress_out, ',')
    fnp = "peff_ini_nod_normal.inp"
elseif regime == "strike-slip" # strike-slip, s11 = sxx, s22 = szz, s33 = syy
    stress_out = hcat(centroids[:,1], s1 .+ p, s3 .+ p, s2 .+ p, t, t, t)
    #writedlm(join([dir_out, "initial_stress_rotated_criticallyStressed_mu0.6_aphi", a_phi, ".inp"]), stress_out, ',')
    writedlm(join([dir_out, "initial_stress_rotated_sed_mu0.3_bas_mu0.6_aphi", a_phi, ".inp"]), stress_out, ',')
    #fnp = join(["peff_ini_nod_aphi", a_phi, ".inp"])
    fnp = join(["peff_ini_nod_aphi", a_phi, "_sed_mu0.3_bas_mu0.6.inp"])
end
writedlm(join([dir_out, fnp]), hcat(nodes_inp[:, 1], map(x -> @sprintf("%.2f", x), peff_ini_nod)), ",\t")

# Save centroids and nodal depths
writedlm(join([dir_out, "element_centroids.dat"]), centroids, ',')
writedlm(join([dir_out, "nodal_depth.dat"]), hcat(nodes_inp[:, 1], nodes_inp[:, end]), ',')

# Save Ranger-Pico element ids
id_rp = vcat(elements_inp[1][:, 1], elements_inp[4][:, 1])
writedlm(join([dir_out, "id_ranger_pico.dat"]), id_rp)
id_rpo = vcat(elements_inp[1][:, 1], elements_inp[4][:, 1], elements_inp[6][:, 1])
writedlm(join([dir_out, "id_ranger_pico_overPico.dat"]), id_rpo)
writedlm(join([dir_out, "id_pico.dat"]), elements_inp[4][:, 1])
writedlm(join([dir_out, "id_overPico.dat"]), elements_inp[6][:, 1])
writedlm(join([dir_out, "id_bottomRanger.dat"]), elements_inp[2][:, 1])

# 5. PLOTS
# Find center of the grid
xc, yc = mean(nodes_inp[:, 2]), mean(nodes_inp[:, 3])
target_xyz = Matrix{Float64}(undef, 76, 3)
target_xyz[:, 1] .= xc
target_xyz[:, 2] .= yc
target_xyz[:, 3] = LinRange(0, -15000, 76)
D = pairwise(Euclidean(), centroids[:, 2:end], target_xyz; dims=1)
id_closest_els = vec(mapslices(argmin, D; dims=1))

f = Figure(size=(500, 400))
ax = Axis(f[1, 1], xlabel = L"\text{Stress [MPa]}", ylabel = L"\text{Depth [km]}", xlabelsize = 16, ylabelsize = 16, xticks=[0, -100, -200, -300, -400, -500, -600], yticks=[0, 2.5, 5, 7.5, 10, 12.5, 15])
poly!(ax, Point2f.([0, -600, -600, 0, 0], [15, 15, 2.9, 2.9, 15]), color = RGBAf(0.8, 0.8, 0.8, 0.5))
lines!(ax, [-600, 0], [2.9, 2.9], color = RGBAf(0.5, 0.5, 0.5, 1), linewidth = 0.75, linestyle = :solid)
text!(ax, -300, 2.0, text=L"\text{Sediments, mobilized } \mu=0.3", fontsize=14, color = :black)
text!(ax, -300, 5.0, text=L"\text{Basement, mobilized } \mu=0.6", fontsize=14, color = :black)
scatter!(ax, s1[id_closest_els]./1e6, centroids[id_closest_els, end]./-1000, markersize = 8, color=RGBAf(0,0,0,0), strokecolor=:black, strokewidth=0.5, label=L"S_1")
scatter!(ax, s2[id_closest_els]./1e6, centroids[id_closest_els, end]./-1000, markersize = 8, color=RGBAf(0,0,0,0), strokecolor=:blue, strokewidth=0.5, label=L"S_2")
scatter!(ax, s3[id_closest_els]./1e6, centroids[id_closest_els, end]./-1000, markersize = 8, color=RGBAf(0,0,0,0), strokecolor=:red, strokewidth=0.5, label=L"S_3")
ylims!(ax, 15, 0)
xlims!(ax, 0, -550)
axislegend(position = :lb, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8)) 
f

x_plot = LinRange(0, 20, 41)
mycmap=cgrad(:nipy_spectral)
f = Figure(size=(900, 400))
ax = Axis(f[1, 1], xlabel = L"A_\phi \text{ [-]}", ylabel = L"\text{Depth [km]}", xlabelsize = 16, ylabelsize = 16)
band!(ax, [0, 1, 2, 3], [0, 0, 0, 0], [15, 15, 15, 15], color=[0, 1, 2, 3], colormap=reverse(cgrad(:roma)))
lines!(ax, [1, 1], [0, 15], color = RGBAf(0.5, 0.5, 0.5, 1), linewidth = 0.5)
lines!(ax, [2, 2], [0, 15], color = RGBAf(0.5, 0.5, 0.5, 1), linewidth = 0.5)
text!(ax, 0.3, 0.7, text=L"\text{Normal}", fontsize=14, color = :black)
text!(ax, 1.3, 0.7, text=L"\text{Strike slip}", fontsize=14, color = :black)
text!(ax, 2.3, 0.7, text=L"\text{Reverse}", fontsize=14, color = :black)
sca1 = scatter!(ax, stress_data[1][:, 14], stress_data[1][:, 3], markersize = 8, color=RGBAf(0,0,0,0), strokecolor=:black, strokewidth=0.5, label=L"\text{FlatMaxwell}")
sca2 = scatter!(ax, stress_data[2][:, 4], stress_data[2][:, 3], markersize = 8, color=RGBAf(0,0,0,0), strokecolor=:grey, strokewidth=0.5, label=L"\text{Hardebeck FM}")
sca3 = scatter!(ax, stress_data[3][:, 4], stress_data[3][:, 3], markersize = 8, color=RGBAf(0,0,0,0), strokecolor=:red, strokewidth=0.5, label=L"\text{Luttrell 2017}")
sca4 = scatter!(ax, stress_data[4][:, 4], stress_data[4][:, 3], markersize = 8, color=RGBAf(0,0,0,0), strokecolor=:green, strokewidth=0.5, label=L"\text{YH14 K}")
sca5 = scatter!(ax, stress_data[5][:, 4], stress_data[5][:, 3], markersize = 8, color=RGBAf(0,0,0,0), strokecolor=:blue, strokewidth=0.5, label=L"\text{YHSM 2013}")
#lin1 = lines!(ax, predict(lm_aphi, (x=x_plot,)), x_plot, color =RGBAf(0,0,0,1), linewidth = 1.5, label=L"P_1 \text{ model fit}")
ylims!(ax, 15, 0)
xlims!(ax, 0, 3)
axislegend(position = :lb, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8)) 

ax = Axis(f[1, 2], xlabel = L"\text{Stress [MPa]}", xlabelsize = 16)
#sca1 = scatter!(ax, stress_data[1][:, 20], stress_data[1][:, 3], markersize = 8, color=RGBAf(0,0,0,0), strokecolor=RGBAf(0,0,0,1), strokewidth=0.5, label=L"S_1 \; (\text{FlatMaxwell})")
#sca2 = scatter!(ax, stress_data[1][:, 19], stress_data[1][:, 3], markersize = 8, color=RGBAf(0,0,0,0), strokecolor=RGBAf(0.5,0.5,0.5,1), strokewidth=0.5, label=L"S_2")
#sca3 = scatter!(ax, stress_data[1][:, 18], stress_data[1][:, 3], markersize = 8, color=RGBAf(0,0,0,0), strokecolor=RGBAf(0.8,0.8,0.8,1), strokewidth=0.5, label=L"S_3")
#lin1 = lines!(ax, predict(lm_s1, (x=x_plot, x2=x_plot.^2, x3=x_plot.^3, x4=x_plot.^4,)), x_plot, color =RGBAf(0,0,0,1), linewidth = 1.5, label=L"P_4 \; (\text{model fit})")
#lin2 = lines!(ax, predict(lm_s2, (x=x_plot, x2=x_plot.^2, x3=x_plot.^3, x4=x_plot.^4,)), x_plot, color =RGBAf(0.5,0.5,0.5,1), linewidth = 1.5, label=L"P_4")
#lin3 = lines!(ax, predict(lm_s3, (x=x_plot,)), x_plot, color = RGBAf(0.8,0.8,0.8,1), linewidth = 2, label=L"P_1")
#lin4 = lines!(ax, predict(lm_sv, (x=x_plot,)), x_plot, color = RGBAf(0.8, 0.8, 0.8, 1.0), linewidth = 1) equal to s3
if regime == "afit" || regime == "reverse"
    s1_label = L"S_\text{H} \text{ (crit. stress)}"
    s2_label = L"S_\text{h} \text{ (crit. stress)}"
    s3_label = L"\rho g z \text{, } \rho = 2300 \text{ kg/m}^3"
elseif regime == "normal"
    s2_label = L"S_\text{H} \text{ (crit. stress)}"
    s3_label = L"S_\text{h} \text{ (crit. stress)}"
    s1_label = L"\rho g z \text{, } \rho = 2300 \text{ kg/m}^3"
elseif regime == "strike-slip"
    s1_label = L"S_\text{H} \text{ (crit. stress)}"
    s3_label = L"S_\text{h} \text{ (crit. stress)}"
    s2_label = L"\rho g z \text{, } \rho = 2300 \text{ kg/m}^3"
end
lin4 = lines!(ax, s3_plot, z_plot*-1, color = RGBAf(1.0, 0.0, 0.0, 1), linewidth = 1.5, linestyle = :dashdot, label=s3_label)
lin5 = lines!(ax, s2_plot, z_plot*-1, color = RGBAf(0.7, 0.0, 0.0, 1), linewidth = 1.5, linestyle = :dashdot, label=s2_label)
lin6 = lines!(ax, s1_plot, z_plot*-1, color = RGBAf(0.5, 0.0, 0.0, 1), linewidth = 1.5, linestyle = :dashdot, label=s1_label)
lin7 = lines!(ax, p_plot, z_plot*-1, color = :blue, linewidth = 1, linestyle = :solid, label=L"\rho_\text{w} g z")
ylims!(ax, 15, 0)
xlims!(ax, 0, -750)
axislegend(position = :rt, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8)) # rt = right-top
f

# Save
# save(join([figs_dir, "stress_state.pdf"]), f, pt_per_unit = 1)