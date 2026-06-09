# include("config.jl") at the top of every script; nothing else hardcodes a path.
const REPO_DIR  = @__DIR__
import Pkg
Pkg.activate(REPO_DIR)

const TOOLS_DIR = joinpath(REPO_DIR, "jl-tools", "modules")

# --- Dataset (separate Zenodo deposit) ---------------------------------------
# Resolve DATA_DIR: local override file → env var → ../data default (downloader gets the default).
const _LOCAL_CFG = joinpath(REPO_DIR, "config_local.jl")   # gitignored; not shipped
isfile(_LOCAL_CFG) && include(_LOCAL_CFG)

const DATA_DIR = @isdefined(LOCAL_DATA_DIR) ? LOCAL_DATA_DIR : get(ENV, "WILMINGTON_DATA", normpath(joinpath(REPO_DIR, "..", "data")))

datapath(p...) = joinpath(DATA_DIR, p...)

# --- Data roots (match the deposit tree; see data/README.md) -----------------
# Inputs
const GRIDS_ABAQUS = datapath("grids", "abaqus")            # mesh .inp files
const GRIDS_MRST   = datapath("grids", "mrst")              # G structs, node/face maps, StateFunction
const GEO_DIR      = datapath("geo")                        # cal_outline.dat
const FOR_MRST_FLUID = datapath("for_mrst_fluid_data")      # BENCH_SPE9_*.DATA deck

# Raw well records (CalGEM), monthly
const RAW_WELL_DATA = datapath("raw_well_data")
const RAW_INJ  = joinpath(RAW_WELL_DATA, "Injection")
const RAW_PROD = joinpath(RAW_WELL_DATA, "Production")

# Intermediates (produced by scripts, then read back)
const FOR_MRST_FLOWRATES  = datapath("for_mrst_flowrates")  # *_m3_s_mrst_trimester.dat
const FOR_ABAQUS_STRESS   = datapath("for_abaqus_stress", "initial")   # initial_stress_*, peff_*, centroids, ...
const FOR_ABAQUS_PRESSURE = datapath("for_abaqus_pressure")
const NODAL_DP_DIR = joinpath(FOR_ABAQUS_PRESSURE, "nodal_dp_1936_to_2020_tri")

# Simulation output (curated: MRST result_data, Abaqus extracted/)
const SIM_MRST = datapath("simulation_output", "mrst")
const MRST_RUN = "JS_LS_lowRs_props_owc1800_fTmultPermBot_BCpvMult_maxDBHPInf_minBHP0MPa_dtTrimester_t177"
const MRST_RESULT_DATA = joinpath(SIM_MRST, MRST_RUN, "result_data")   # states_p_Pa.mat, states_t_s.mat, ...

const SIM_ABAQUS = datapath("simulation_output", "abaqus")  # then variant + "extracted":
const ABAQUS_VARIANTS = ("plastic-130", "plastic-120")      # critically-stressed (CSM), strike-slip heterogeneous
abaqus_extracted(variant) = joinpath(SIM_ABAQUS, variant, "extracted")

# SCEC CSM is third-party and NOT in the deposit (cited at source). Only
# abaqus_create_initial_stress.jl needs it, and only when you regenerate the
# stress intermediates locally. 

# --- Generated outputs (produced by running the code) ------------------------
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