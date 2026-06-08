# `include("config.jl")` at the top of every script; nothing else hardcodes a path.
const REPO_DIR  = @__DIR__
import Pkg
Pkg.activate(REPO_DIR)        

const TOOLS_DIR = joinpath(REPO_DIR, "jl-tools", "modules")   

# --- Dataset (separate Zenodo deposit) ---------------------------------------
# Resolve DATA_DIR: local override file → env var → ../data default (downloader gets the default).
const _LOCAL_CFG = joinpath(REPO_DIR, "config_local.jl")   # gitignored; not shipped
isfile(_LOCAL_CFG) && include(_LOCAL_CFG)

const DATA_DIR = @isdefined(LOCAL_DATA_DIR) ? LOCAL_DATA_DIR : get(ENV, "WILMINGTON_DATA", normpath(joinpath(REPO_DIR, "..", "data")))

# Assembled input roots — these consolidate the scattered absolute paths.
const ABAQUS_DIR = joinpath(DATA_DIR, "abaqus_grids")   
const MRST_GRID  = joinpath(DATA_DIR, "mrst", "grid")   
const MRST_PRESS = joinpath(DATA_DIR, "mrst", "pressure")
const SCEC_DIR   = joinpath(DATA_DIR, "scec_csm")       

datapath(p...) = joinpath(DATA_DIR, p...)

# --- Generated outputs (produced by running the code) ---------
const OUT_DIR  = get(ENV, "WILMINGTON_OUT",  joinpath(REPO_DIR, "output"))
const FIGS_DIR = get(ENV, "WILMINGTON_FIGS", joinpath(REPO_DIR, "figures"))
foreach(d -> isdir(d) || mkpath(d), (OUT_DIR, FIGS_DIR))

isdir(DATA_DIR) || @warn "DATA_DIR not found — set \$WILMINGTON_DATA or place data at ../data" DATA_DIR

# --- Tool modules ------------------------------------------------------------
const TOOLS = Dict(
    :arrayops  => "arrayOps.jl",          :fileio    => "fileIO.jl",
    :abaqus    => "abaqus.jl",            :mrst      => "mrst.jl",
    :cmg       => "cmg.jl",               :faultops  => "faultOperations.jl",
    :gmsh      => "gmsh.jl",              :pvt       => "reservoirFluidsPVT.jl",
    :makieplot => "makiePolyPlot.jl",
)
tool(k::Symbol) = joinpath(TOOLS_DIR, TOOLS[k])