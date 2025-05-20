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

elseif platform == :aarch64
    dir_data = "/Users/lluis/Dropbox (MIT)/home_copy/julia/data/wilmington/"
    dir_out = "/Users/lluis/Dropbox (MIT)/home_copy/julia/output/wilmington/"
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
datum_name = "WGS84"

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
coord_inj_filled[:, [1, 2]] = coord_inj_filled[:, [2, 1]]   # swap lat and lon to agree with EPSG:4326
coord_inj_utm = trans.(eachrow(coord_inj_filled))           # produces a vector of tuples

coord_prod[:, [1, 2]] = coord_prod[:, [2, 1]]   
coord_prod_utm = trans.(eachrow(coord_prod))

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

# Space: consolidate injectors/producers in the same cell centroid

# Time: Reduce number of timesteps if desired 
