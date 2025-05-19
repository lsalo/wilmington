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
wgs84 = Proj.CRS("EPSG:4326")                             
if datum_name == "WGS84"
    utm11_datum = Proj.CRS("EPSG:32611")    
elseif datum_name == "NAD27"
    utm11_datum = Proj.CRS("EPSG:26711")     
end

# Transform
trans = Proj.Transformation(wgs84, utm11_datum)
coord_inj_utm = trans.(eachrow(coord_inj))Y)
