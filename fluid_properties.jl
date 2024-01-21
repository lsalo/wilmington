#
# Generate fluid properties (PVTO, PVDG, PVTW) for Wilmington field
# Reservoir simulation using a black oil model in MRST
# We use ECLIPSE "FIELD" units for fluid propeties
# Water is a brine with 30,000ppm salinity (French & Burchfield, 1990)
# Gas is a mixture of 95% methane and 5% CO2 (Sankur et al., 1986)
# Oil averages 16 deg. API (Nelson, 1976)
#

# Startup
using DelimitedFiles
using Interpolations
using Statistics

save_data = false

# Paths and variables
in_dir = "C:\\Users\\lsalo\\matlab\\mrst-dev\\mrst-ls\\ls-proj\\wilmington\\input\\"
PVTO = readdlm(join([in_dir, "brine_PVTO_T323.15_P1P500_Sbrine30000ppm.txt"]))
units = "FIELD"                                             # "FIELD" or "METRIC"

kgm3_to_lbft3   = 0.062428
temp_std_C      = 15.56
p_std_bar       = 1.01325
api             = 18;
temp_r_C        = 50
sal_ppm         = 30000
bar_to_psi      = 14.503773773

# 1. Surface densities (ECLIPSE standard conditions)
density_gas_kgm3 = 0.738                                    # AGA8 equation for natural gas dens. (https://www.unitrove.com/engineering/tools/gas/natural-gas-density)
density_brine_kgm3 = 1018.46                                # Batzle & Wang (1992)
density_oil_kgm3 = 999.15*(141.5/(api + 131.5))             # API gravity number definition. pure H20 density = 999.15 kg/m3
density = [density_oil_kgm3, density_brine_kgm3, density_gas_kgm3]
if units == "FIELD"
    density = density.*kgm3_to_lbft3
end


# 2. PVTW (Brine PVT fcns)
# We generate brine props using pvtBrineWithCO2BlackOil (MATLAB), and then convert from PVTO format to PVTW and units
# Note: PVTO is assumed to be in ECLIPSE "METRIC" units
# Note: PVTO is assumed to have at least one value above the maximum requested PVTW p value

# Calculate
p_out_bar = [1, 100, 200, 300, 400, 500]
PVTW = Matrix{Float64}(undef, length(p_out_bar), 5)
PVTW[:,1] = p_out_bar
pv = Float64.(PVTO[2:end, 1])
for (id, p_val) in enumerate(p_out_bar)
    # interpolate for B, mu
    fb = linear_interpolation(pv, Float64.(PVTO[2:end, 2]))
    PVTW[id, 2] = fb(p_val)                                              # B [rm3/Sm3]
    fm = linear_interpolation(pv, Float64.(PVTO[2:end, 3]))
    PVTW[id, 4] = fm(p_val)                                              # mu [cP]
    
    # Calculate C (compressibility) and Cv (viscosibility)
    index = searchsortedfirst(pv, p_val)
    if index == 1
        index = [1,2] 
    else
        index = [index-1, index+1]
    end 
    bw_calc     = PVTO[1 .+ index, 2]
    p_calc      = pv[index]
    bw_avg      = mean(bw_calc)
    PVTW[id, 3] = (abs.(diff(bw_calc)/diff(p_calc)).*bw_avg)[1]         # C [1/bar]
    mu_calc     = PVTO[1 .+ index, 3]
    mu_avg      = mean(mu_calc)
    PVTW[id, 5] = (abs.(diff(mu_calc)/diff(p_calc)).*mu_avg)[1]         # Cv [1/bar]
end

# Convert units if not "METRIC"
if units == "FIELD"
    PVTW[:,1] = PVTW[:,1].*bar_to_psi  
    PVTW[:,3] = PVTW[:,3]./bar_to_psi
    PVTW[:,5] = PVTW[:,5]./bar_to_psi
end


# 3. PVDG (dry gas)

# 4. PVTO (live oil)

# Save output to text file
if save_data
    DelimitedFiles.writedlm(join([in_dir, "PVTW.dat"]), PVTW, '\t') 
    # PVDG
    # PVTO
end