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

# ─────────────────────────────
# ▶ DIRECTORIES AND DATA ◀
# ─────────────────────────────
platform = Sys.ARCH
if platform == :x86_64 

elseif platform == :aarch64
    dir_data = "/Users/lluis/Dropbox (MIT)/home_copy/julia/data/wilmington/"
    dir_out = "/Users/lluis/Dropbox (MIT)/home_copy/julia/output/wilmington/"
end

# Read well coordinates
coord_inj = readdlm(join([dir_data, "/raw_flowrates/Injection/coord.dat"]), '\t', Float64, '\n')


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

# Transform
trans = Proj.Transformation(wgs84, utm11_datum)
coord_inj_filled = fill_zero_rows(coord_inj; method="nearest_average") 
coord_inj_filled[:, [1, 2]] = coord_inj_filled[:, [2, 1]]   # swap lat and lon to agree with EPSG:4326
coord_inj_utm = trans.(eachrow(coord_inj_filled))           # produces a vector of tuples

