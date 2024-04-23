#=
Script to initialize stresses at Wilmington using FlatMaxwell model from SCEC CSM.

Important notes:
    1. We removed character # from FlatMaxwell.csv row 48, where the header is.
    2. In our convention, S1 = maximum compressive, S2 = intermediate and S3 = minimum compressive or most tensional. Note that this is different than the convention used in FlatMaxwel.csv.

TBD: Add another option to compute based on frictional limit, follow eq. in chapter 7 of Lund Snee thesis. To do this, the idea is to get average A_phi from all CSM stress models, then select stress regime in Wilmington area (likely reverse), and compute Sh and SH from Sv. Sv can be computed from avg density used in Abaqus to facilitate geostatic equilibrium. Option b would be taking it from FlatMaxwell model.
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

#using DelimitedFiles
#using Statistics
#using Colors
#using ColorSchemes


#lib_dir = "C:/Users/lsalo/julia/jl-tools/"
#includet(joinpath(lib_dir, "modules", "faultOperations.jl"))
#using .FaultOperations

# Directories
dir_data = "C:/Users/lsalo/Wilmington/2023_SCEC_CSM/SCEC_CSM_archive_2023/"
dir_out = "C:/Users/lsalo/julia/output/wilmington/"
figs_dir = "C:/Users/lsalo/julia/graphics/figs/wilmington/test/" 

# Read dModel data 
stress_data_all = CSV.File(join([dir_data, "FlatMaxwell.csv"]), header=48, select=1:20, skipto=49)

# Get subset: Wilmington area
min_lon, max_lon = -118.41947, -117.93033 
min_lat, max_lat = 33.58679, 33.93842
row_ids = [min_lon <= row.LON <= max_lon && min_lat <= row.LAT <= max_lat && 0 <= row.DEP <= 16 for row in stress_data_all]
stress_data = Matrix{Float64}(undef, sum(row_ids), length(stress_data_all.names))
k = 1
for (i, row) in enumerate(stress_data_all)
    if row_ids[i] == 1
        stress_data[k, :] = [getproperty(row, name) for name in stress_data_all.names]
        k += 1
    end
end

# Least squares linear fit through the origin
x = stress_data[:,3]
lm_s1 = glm(@formula(y ~ x + x2 + x3 + x4 + 0), (x=x, x2=x.^2, x3=x.^3, x4=x.^4, y=stress_data[:,20]), Normal(), IdentityLink())
lm_s2 = glm(@formula(y ~ x + x2 + x3 + x4 + 0), (x=x, x2=x.^2, x3=x.^3, x4=x.^4, y=stress_data[:,19]), Normal(), IdentityLink())
lm_s3 = glm(@formula(y ~ x + 0), (x=x, y=stress_data[:,18]), Normal(), IdentityLink())
#lm_sv = glm(@formula(y ~ x + 0), (x=x, y=stress_data[:,9]), Normal(), IdentityLink()) # basically equal to s3

# Predict values based on sv and frictional constraints (Lund Snee Thesis, Chapter 7)


# Calculate element centroids for all elements in geomechanics mesh
# get inp file data
# find nodes of each element in reservoir elset
# nodal coordinates of each element and calculate centroids

# Assign stress based on depth

# Save file for use with Abaqus
# initial_stress_wilmington_rotated.inp

# Plots 
x_new = LinRange(0, 20, 41)
f = Figure()
ax = Axis(f[1, 1], xlabel = L"\text{Stress [MPa]}", ylabel = L"\text{Depth [km]}", xlabelsize = 16, ylabelsize = 16)
sca1 = scatter!(ax, stress_data[:, 20], stress_data[:, 3], markersize = 8, color=RGBAf(0,0,0,0), strokecolor=:black, strokewidth=0.5, label=L"\sigma_\text{11} \; (\text{FlatMaxwell})")
sca2 = scatter!(ax, stress_data[:, 19], stress_data[:, 3], markersize = 8, color=RGBAf(0,0,0,0), strokecolor=:blue, strokewidth=0.5, label=L"\sigma_\text{22}")
sca3 = scatter!(ax, stress_data[:, 18], stress_data[:, 3], markersize = 8, color=RGBAf(0,0,0,0), strokecolor=:grey, strokewidth=0.5, label=L"\sigma_\text{33}")
lin1 = lines!(ax, predict(lm_s1, (x=x_new, x2=x_new.^2, x3=x_new.^3, x4=x_new.^4,)), x_new, color = :black, linewidth = 1.5, label=L"P_4 \; (\text{model fit})")
lin2 = lines!(ax, predict(lm_s2, (x=x_new, x2=x_new.^2, x3=x_new.^3, x4=x_new.^4,)), x_new, color = :blue, linewidth = 1.5, label=L"P_4")
lin3 = lines!(ax, predict(lm_s3, (x=x_new,)), x_new, color = RGBAf(0.6, 0.6, 0.6, 1.0), linewidth = 2, label=L"P_1")
#lin4 = lines!(ax, predict(lm_sv, (x=x_new,)), x_new, color = RGBAf(0.8, 0.8, 0.8, 1.0), linewidth = 1) equal to s3
ylims!(ax, 15, 0)
xlims!(ax, 0, -600)
axislegend(position = :rt, framevisible=false) # rt = right-top
f