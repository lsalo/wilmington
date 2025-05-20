#=
Script to create well information for the Wilmington model, from raw data from CalGEMS
=#

# Cd into current directory and activate the Julia environment in this directory
cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__)

# Startup - modules by others
using Revise
using DelimitedFiles
using StatsBase
using Proj
using MAT

# Startup - modules by me
#lib_dir = "C:/Users/lsalo/julia/jl-tools/"
lib_dir = "/Users/lluis/Documents/julia/jl-tools/"
includet(joinpath(lib_dir, "modules", "arrayOps.jl"))
using .ArrayOps
includet(joinpath(lib_dir, "modules", "reservoirfluidspvt.jl"))
using .ReservoirFluidsPVT


# ─────────────────────────────────
# ▶ DIRECTORIES, OPTIONS AND DATA ◀
# ─────────────────────────────────
platform = Sys.ARCH
if platform == :x86_64 
    dir_grid    = "C:/Users/lsalo/matlab/mrst-dev/mrst-ls/ls-proj/wilmington/grid/"
elseif platform == :aarch64
    dir_data = "/Users/lluis/Dropbox (MIT)/home_copy/julia/data/wilmington/"
    dir_out = "/Users/lluis/Dropbox (MIT)/home_copy/julia/output/wilmington/"
    dir_figs = "/Users/lluis/Dropbox (MIT)/home_copy/julia/graphics/wilmington/test/"
    dir_grid = "/Users/lluis/Documents/MATLAB/mrst-dev/mrst-ls/ls-proj/wilmington/grid/"
end

# Options
tstep_freq = 3 # frequency of timesteps to save (1 = all (monthly), 2 = bi-monthly, 3 = quarterly, 6 = semi-annually)

# Read well info
coord_inj = readdlm(join([dir_data, "/raw_flowrates/Injection/coord.dat"]), '\t', Float64, '\n')
inj_monthly = readdlm(join([dir_data, "/raw_flowrates/Injection/waterInjection.dat"]), '\t', Float64, '\n') # tsteps x wells

coord_prod = readdlm(join([dir_data, "/raw_flowrates/Production/coord.dat"]), '\t', Float64, '\n')
prod_water_monthly = readdlm(join([dir_data, "/raw_flowrates/Production/waterProduction.dat"]), '\t', Float64, '\n') # tsteps x wells
prod_oil_monthly = readdlm(join([dir_data, "/raw_flowrates/Production/oilProduction.dat"]), '\t', Float64, '\n') # tsteps x wells
prod_gas_monthly = readdlm(join([dir_data, "/raw_flowrates/Production/gasProduction.dat"]), '\t', Float64, '\n') # tsteps x wells


# ─────────────────────────────
# ▶ CONVERT TO UTM COORDINATES ◀
# ─────────────────────────────
datum_name = "WGS84" # Appears to be the most consistent datum for the grid

# Define coordinate systems
wgs84 = Proj.CRS("EPSG:4326")     # expects (lat, lon)       
if datum_name == "WGS84"
    utm11_datum = Proj.CRS("EPSG:32611")    
elseif datum_name == "NAD27"
    utm11_datum = Proj.CRS("EPSG:26711")     
end
trans = Proj.Transformation(wgs84, utm11_datum)

# Transform
coord_inj_filled = fill_zero_rows(coord_inj; method="nearest_average") 
coord_inj_filled[:, [1, 2]] = coord_inj_filled[:, [2, 1]]               # swap lat and lon to agree with EPSG:4326
coord_inj_utm = trans.(eachrow(coord_inj_filled))                       # produces a vector of tuples
coord_inj_utm = hcat(first.(coord_inj_utm), last.(coord_inj_utm))       # convert to matrix

coord_prod[:, [1, 2]] = coord_prod[:, [2, 1]]   
coord_prod_utm = trans.(eachrow(coord_prod))
coord_prod_utm = hcat(first.(coord_prod_utm), last.(coord_prod_utm))

# ──────────────────────────────────────────────────
# ▶ CONVERT FLOWRATES TO M3/S (MRST uses SI units) ◀
# ──────────────────────────────────────────────────
inj_monthly_m3 = convert_volume_units(inj_monthly, "BBL", "M3")             # convert total monthly injection to m3
avg_monthly_days = 365.25 / 12
inj_monthly_m3_s = inj_monthly_m3 ./ (avg_monthly_days * 24 * 3600)         # monthly injection rate in m3/s

prod_water_monthly_m3 = convert_volume_units(prod_water_monthly, "BBL", "M3") 
prod_water_monthly_m3_s = prod_water_monthly_m3 ./ (avg_monthly_days * 24 * 3600) 
prod_oil_monthly_m3 = convert_volume_units(prod_oil_monthly, "BBL", "M3")  
prod_oil_monthly_m3_s = prod_oil_monthly_m3 ./ (avg_monthly_days * 24 * 3600)
prod_gas_monthly_m3 = convert_volume_units(prod_gas_monthly, "MSCF", "M3") 
prod_gas_monthly_m3_s = prod_gas_monthly_m3 ./ (avg_monthly_days * 24 * 3600) 

# ───────────
# ▶ UPSCALE ◀
# ───────────
# Import MRST grid centroids
fn_G = "wilmington_ranger_G_withGeometry.mat"
G = matread(join([dir_grid, fn_G]))["G"]
centroids = G["cells"]["centroids"]

# Space: consolidate injectors/producers in the same cell centroid
well_ids = find_ids_from_pairwise_distance(coord_inj_utm, centroids[:, 1:2])
inj_monthly_m3_s_wid = Matrix{Union{Int, Float64}}(undef, length(well_ids), size(inj_monthly)[1]+1)
inj_monthly_m3_s_wid[:, 1] = well_ids
inj_monthly_m3_s_wid[:, 2:end] = transpose(inj_monthly_m3_s) # NOTE: columns are well numbers after transpose
inj_monthly_m3_s_cid = reduce_matrix_with_repeated_ids(inj_monthly_m3_s_wid, "sum"; out_matrix_type="Int-Float") 

# Time: Reduce number of timesteps if desired, and save to file
if tstep_freq == 1

end


# ───────────────────────
# ▶ FAST VISUALIZATIONS ◀
# ───────────────────────
using CairoMakie
CairoMakie.activate!()

# Plot cell entroids and wells and coastline
coast_data = readdlm(join([dir_data, "cal_outline.dat"]))
cst = [coast_data[coast_data[:, 3] .== i, :] for i in 1:18] # third column is "category"
coastline = [cst[6], cst[17], cst[18]]

f = Figure(size=(600, 500))
ax = Axis(f[1, 1], xlabel = L"x \text{ [m], UTM 11N}", ylabel = L"y \text{ [m], UTM 11N}", xlabelsize = 16, ylabelsize = 16)
scatter!(ax, centroids[:, 1], centroids[:, 2], markersize = 0.5, color = :black, label=L"\text{Ranger cell centroids}")
scatter!(ax, coord_prod_utm[:, 1], coord_prod_utm[:, 2], markersize = 2, color = :red, label=L"\text{Producers}")
scatter!(ax, coord_inj_utm[:, 1], coord_inj_utm[:, 2], markersize = 2, color = :blue, label=L"\text{Injectors}")
count = 1
for line in coastline
    if count ==1
        lines!(ax, line[:, 1], line[:, 2], color=:black, linewidth=0.5, linestyle=:solid, label=L"\text{Coastline}")
        count += 1
    else
        lines!(ax, line[:, 1], line[:, 2], color=:black, linewidth=0.5, linestyle=:solid)
    end
end
xlims!(ax, 3.78e5, 4.0e5)
ylims!(ax, 3.728e6, 3.752e6)
axislegend(position = :lt, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8), labelsize=14) 
f
save(join([dir_figs, "cell_centroids_and_wellsWGS84.pdf"]), f, pt_per_unit = 1)