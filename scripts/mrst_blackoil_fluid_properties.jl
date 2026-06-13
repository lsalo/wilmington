#=
Black-oil fluid properties for the Wilmington reservoir simulation (MRST), ECLIPSE FIELD units.
  Water: 30,000 ppm brine (French & Burchfield, 1990)
  Gas:   95% CH4 / 5% CO2 (Sankur et al., 1986)
  Oil:   Ranger Fm, up to ~25° API (Mayuga, 1970; Nelson, 1976; California Division of Oil and Gas, 1992; Montgomery, 1998)

Scope:
  - This script REGENERATES the gas and brine surface densities, the brine table (PVTW),
    and the dry-gas table (PVDG) used in the deposited MRST .DATA deck.
  - The live-oil table (PVTO) in the .DATA was adapted from the SPE9 benchmark and
    hand-calibrated to better match Wilmington oil; it is NOT regenerated here. The
    oilPVTBW1992 call below produces a comparable PVTO for reference/validation only. Oil density and viscosity in the .DATA are validated against measured ranges (see Table 1 and Fig. 6 in the paper).
  - The oil surface density in the .DATA (56.3965 lb/ft^3) corresponds to API 25 (set below).
=#
include(joinpath(@__DIR__, "..", "config.jl"))

using DelimitedFiles, Interpolations, Statistics
include(tool(:pvt)); using .ReservoirFluidsPVT   # gasPVTBW1992, oilPVTBW1992

const SAVE_TO_DEPOSIT = false      # true → (re)write the PVT tables into the deposit

# Input brine PVTO table (from MATLAB pvtBrineWithCO2BlackOil)
b_PVTO = readdlm(joinpath(FOR_MRST_FLUID, "brine_PVTO_T323.15_P1P500_Sbrine30000ppm.txt"))
units  = "FIELD"                   # "FIELD" or "METRIC"

kgm3_to_lbft3 = 0.062428
api           = 25
temp_r_C      = 50
sal_ppm       = 30000
bar_to_psi    = 14.503773773

# 1. Surface densities (ECLIPSE standard conditions)
density_gas_kgm3   = 0.738                                   # AGA8 natural-gas density
density_brine_kgm3 = 1018.46                                 # Batzle & Wang (1992)
density_oil_kgm3   = 999.15*(141.5/(api + 131.5))            # API gravity definition
density = [density_oil_kgm3, density_brine_kgm3, density_gas_kgm3]
units == "FIELD" && (density = density .* kgm3_to_lbft3)

# 2. PVTW (brine), interpolated from the brine PVTO table; C and Cv by finite difference
p_out_bar = [1, 100, 200, 300, 400, 500]
PVTW = Matrix{Float64}(undef, length(p_out_bar), 5)
PVTW[:, 1] = p_out_bar
pv = Float64.(b_PVTO[2:end, 1])
for (id, p_val) in enumerate(p_out_bar)
    fb = linear_interpolation(pv, Float64.(b_PVTO[2:end, 2])); PVTW[id, 2] = fb(p_val)   # B [rm3/Sm3]
    fm = linear_interpolation(pv, Float64.(b_PVTO[2:end, 3])); PVTW[id, 4] = fm(p_val)   # mu [cP]
    index = searchsortedfirst(pv, p_val)
    index = index == 1 ? [1, 2] : [index-1, index+1]
    bw_calc = b_PVTO[1 .+ index, 2]; p_calc = pv[index]; bw_avg = mean(bw_calc)
    PVTW[id, 3] = (abs.(diff(bw_calc)/diff(p_calc)).*bw_avg)[1]                           # C [1/bar]
    mu_calc = b_PVTO[1 .+ index, 3]; mu_avg = mean(mu_calc)
    PVTW[id, 5] = (abs.(diff(mu_calc)/diff(p_calc)).*mu_avg)[1]                           # Cv [1/bar]
end
if units == "FIELD"
    PVTW[:, 1] = PVTW[:, 1] .* bar_to_psi
    PVTW[:, 3] = PVTW[:, 3] ./ bar_to_psi
    PVTW[:, 5] = PVTW[:, 5] ./ bar_to_psi
end

# 3. PVDG (dry gas)
pressure_bar = [1, 10, 20, 30, 50, 75, 100, 130, 160, 200, 250, 300, 350, 400, 450, 500]
gas_density, gas_FVF, gas_viscosity_cP = gasPVTBW1992(density_gas_kgm3, temp_r_C, pressure_bar, units)
PVDG = hcat(pressure_bar, gas_FVF, gas_viscosity_cP)
units == "FIELD" && (PVDG[:, 1] = PVDG[:, 1] .* bar_to_psi)

# 4. PVTO (live oil)
oil_dens_kgm3, oil_Rs, p, oil_FVF, oil_visc_cP = oilPVTBW1992(density_gas_kgm3, api, temp_r_C, pressure_bar, units)
PVTO = hcat(oil_Rs, p, oil_FVF, oil_visc_cP)

# Save PVT tables into the deposit (fluid-model inputs alongside the brine table + .DATA)
if SAVE_TO_DEPOSIT
    writedlm(joinpath(FOR_MRST_FLUID, "PVTW.dat"), PVTW, '\t')
    writedlm(joinpath(FOR_MRST_FLUID, "PVDG.dat"), PVDG, '\t')
    writedlm(joinpath(FOR_MRST_FLUID, "PVTO.dat"), PVTO, '\t')
end