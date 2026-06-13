#=
Well information and well-location map for the Wilmington MRST model, from raw CalGEM data.

Pipeline: read raw monthly injection/production volumes + well lat/lon → drop inactive
wells → transform to UTM zone 11N (WGS84) → convert to SI rates (m^3/s) → assign each
well to an MRST grid cell by nearest centroid → consolidate co-located wells → resample
to the chosen timestep frequency.

Reproduces:
  - for_mrst_flowrates/{inj, prod_water, prod_oil, prod_gas, time}_m3_s_mrst_trimester.dat
    (the well-rate inputs to the reservoir simulation; verified bit-exact against the
     deposited files)
  - Fig. 7b: map view of the simulation domain (cell-centroid depths, well locations, coastline)

Inputs (deposited): raw_well_data/{Injection,Production}/, grids/mrst/wilmington_ranger_G_withGeometry.mat,
geo/cal_outline.dat. Set SAVE_TO_DEPOSIT=true only to regenerate deposited outputs.
=#
include(joinpath(@__DIR__, "..", "config.jl"))

using DelimitedFiles, StatsBase, Proj, MAT
include(tool(:arrayops)); using .ArrayOps                 # fill_zero_rows, find_ids_from_pairwise_distance, resample_matrix, reduce_matrix_with_repeated_ids
include(tool(:pvt));      using .ReservoirFluidsPVT        # convert_volume_units

## Options
tstep_freq = 3            # 1=monthly, 2=bimonthly, 3=quarterly(trimester), 6=semiannual
const SAVE_TO_DEPOSIT = false
figdir = SAVE_TO_DEPOSIT ? datapath("figures", "mrst") : FIGS_DIR
isdir(figdir) || mkpath(figdir)

## ── Read raw well data (CalGEM) ──────────────────────────────────────────────
coord_inj   = readdlm(joinpath(RAW_INJ, "coord.dat"), '\t', Float64, '\n')
inj_monthly = readdlm(joinpath(RAW_INJ, "waterInjection.dat"), '\t', Float64, '\n')   # tsteps × wells, BBL
inj_total   = vec(sum(inj_monthly, dims=1))

coord_prod         = readdlm(joinpath(RAW_PROD, "coord.dat"), '\t', Float64, '\n')
prod_water_monthly = readdlm(joinpath(RAW_PROD, "waterProduction.dat"), '\t', Float64, '\n')
prod_oil_monthly   = readdlm(joinpath(RAW_PROD, "oilProduction.dat"),   '\t', Float64, '\n')
prod_gas_monthly   = readdlm(joinpath(RAW_PROD, "gasProduction.dat"),   '\t', Float64, '\n')
prod_total         = vec(sum(prod_water_monthly + prod_oil_monthly + prod_gas_monthly, dims=1))

# Drop wells with no injection / production
inj_monthly = inj_monthly[:, inj_total .> 0];  coord_inj = coord_inj[inj_total .> 0, :]
prod_water_monthly = prod_water_monthly[:, prod_total .> 0]
prod_oil_monthly   = prod_oil_monthly[:,   prod_total .> 0]
prod_gas_monthly   = prod_gas_monthly[:,   prod_total .> 0]
coord_prod = coord_prod[prod_total .> 0, :]

time = readdlm(joinpath(RAW_INJ, "time.dat"), '\t', Float64, '\n')

## ── Lat/lon → UTM zone 11N ───────────────────────────────────────────────────
datum_name = "WGS84"
wgs84 = Proj.CRS("EPSG:4326")                                   # (lat, lon)
utm11_datum = datum_name == "WGS84" ? Proj.CRS("EPSG:32611") : Proj.CRS("EPSG:26711")  # WGS84 / NAD27
trans = Proj.Transformation(wgs84, utm11_datum)

coord_inj_filled = fill_zero_rows(coord_inj; method="nearest_average")
coord_inj_filled[:, [1, 2]] = coord_inj_filled[:, [2, 1]]       # swap to (lat, lon) for EPSG:4326
coord_inj_utm = trans.(eachrow(coord_inj_filled))
coord_inj_utm = hcat(first.(coord_inj_utm), last.(coord_inj_utm))

coord_prod[:, [1, 2]] = coord_prod[:, [2, 1]]
coord_prod_utm = trans.(eachrow(coord_prod))
coord_prod_utm = hcat(first.(coord_prod_utm), last.(coord_prod_utm))

## ── Flowrates → m³/s (MRST SI units) ─────────────────────────────────────────
avg_monthly_days = 365.25 / 12
to_rate(v) = v ./ (avg_monthly_days * 24 * 3600)
inj_monthly_m3_s        = to_rate(convert_volume_units(inj_monthly,        "BBL",  "M3"))
prod_water_monthly_m3_s = to_rate(convert_volume_units(prod_water_monthly, "BBL",  "M3"))
prod_oil_monthly_m3_s   = to_rate(convert_volume_units(prod_oil_monthly,   "BBL",  "M3"))
prod_gas_monthly_m3_s   = to_rate(convert_volume_units(prod_gas_monthly,   "MSCF", "M3"))

## ── Upscale: assign wells to MRST grid cells, consolidate co-located ─────────
G = matread(joinpath(GRIDS_MRST, "wilmington_ranger_G_withGeometry.mat"))["G"]
centroids = G["cells"]["centroids"]

function consolidate(coord_utm, rate_monthly_m3_s, nt)
    wid = find_ids_from_pairwise_distance(coord_utm, centroids[:, 1:2])
    M = Matrix{Union{Int,Float64}}(undef, length(wid), nt + 1)
    M[:, 1] = wid
    M[:, 2:end] = transpose(rate_monthly_m3_s)          # columns = wells after transpose
    reduce_matrix_with_repeated_ids(M, "sum"; out_matrix_type="Int-Float")
end

nt = size(inj_monthly, 1)
inj_cid        = consolidate(coord_inj_utm,  inj_monthly_m3_s,        nt)
prod_water_cid = consolidate(coord_prod_utm, prod_water_monthly_m3_s, nt)
prod_oil_cid   = consolidate(coord_prod_utm, prod_oil_monthly_m3_s,   nt)
prod_gas_cid   = consolidate(coord_prod_utm, prod_gas_monthly_m3_s,   nt)

wid_inj_mrst  = inj_cid[:, 1]
wid_prod_mrst = prod_water_cid[:, 1]

## ── Resample timesteps ───────────────────────────────────────────────────────
function resample_cid(cid, tstep_freq)
    tstep_freq == 1 && return cid
    nt_new = Int(ceil((size(cid, 2) - 1) / tstep_freq))
    out = Matrix{Union{Int,Float64}}(undef, size(cid, 1), nt_new + 1)
    out[:, 1:2] = cid[:, 1:2]                            # well ids + first column (1936.0)
    out[:, 3:end] = resample_matrix(cid[:, 3:end], tstep_freq)
    out
end
inj_fin        = resample_cid(inj_cid, tstep_freq)
prod_water_fin = resample_cid(prod_water_cid, tstep_freq)
prod_oil_fin   = resample_cid(prod_oil_cid, tstep_freq)
prod_gas_fin   = resample_cid(prod_gas_cid, tstep_freq)
time_fin       = tstep_freq == 1 ? time : time[1:tstep_freq:end]

## ── Save ─────────────────────────────────────────────────────────────────────
if SAVE_TO_DEPOSIT
    tstep_str = Dict(1=>"monthly", 2=>"bimonthly", 3=>"trimester", 6=>"semiannual")[tstep_freq]
    writedlm(joinpath(FOR_MRST_FLOWRATES, "inj_m3_s_mrst_$(tstep_str).dat"),        inj_fin,        '\t')
    writedlm(joinpath(FOR_MRST_FLOWRATES, "prod_water_m3_s_mrst_$(tstep_str).dat"), prod_water_fin, '\t')
    writedlm(joinpath(FOR_MRST_FLOWRATES, "prod_oil_m3_s_mrst_$(tstep_str).dat"),   prod_oil_fin,   '\t')
    writedlm(joinpath(FOR_MRST_FLOWRATES, "prod_gas_m3_s_mrst_$(tstep_str).dat"),   prod_gas_fin,   '\t')
    writedlm(joinpath(FOR_MRST_FLOWRATES, "time_mrst_$(tstep_str).dat"),            time_fin,       '\t')
end

# For internal reproducibility check
#for (name, mat) in (("inj", inj_fin), ("prod_water", prod_water_fin),
#                    ("prod_oil", prod_oil_fin), ("prod_gas", prod_gas_fin))
#    dep = readdlm(joinpath(FOR_MRST_FLOWRATES, "$(name=="inj" ? "inj" : name)_m3_s_mrst_trimester.dat"), '\t', Float64, '\n')
#    println(rpad(name, 12), " maxdiff = ", maximum(abs.(Float64.(mat) .- dep)), "  shape ", size(mat)==size(dep))
#end

## ── Fig7b map: cell centroids, wells, coastline ────────────────────────
using CairoMakie
CairoMakie.activate!()
coast = readdlm(joinpath(GEO_DIR, "cal_outline.dat"))
cst = [coast[coast[:, 3] .== i, :] for i in 1:18]
coastline = [cst[6], cst[17], cst[18]]
f = Figure(size=(600, 425))
ax = Axis(f[1, 1], xlabel=L"x \text{ [m], UTM 11N}", ylabel=L"y \text{ [m], UTM 11N}", xlabelsize=16, ylabelsize=16)
sc = scatter!(ax, centroids[:, 1], centroids[:, 2], markersize=1.5, color=centroids[:, 3], colormap=:tempo)
scatter!(ax, centroids[wid_prod_mrst, 1], centroids[wid_prod_mrst, 2], markersize=2, color=:red,  label=L"\text{Producers}")
scatter!(ax, centroids[wid_inj_mrst, 1],  centroids[wid_inj_mrst, 2],  markersize=2, color=:blue, label=L"\text{Injectors}")
for (i, line) in enumerate(coastline)
    lines!(ax, line[:, 1], line[:, 2], color=:black, linewidth=0.5, label=(i==1 ? L"\text{Coastline}" : nothing))
end
xlims!(ax, 3.78e5, 4.01e5); ylims!(ax, 3.728e6, 3.752e6)
axislegend(position=:rt, framevisible=true, framewidth=0.5, backgroundcolor=RGBAf(1,1,1,0.8), labelsize=14)
Colorbar(f[1, 2], sc, label="Depth [m]", height=Relative(0.5))
display(f)
SAVE_TO_DEPOSIT && save(joinpath(figdir, "cell_centroids_and_wells_WGS84.pdf"), f, pt_per_unit=1)