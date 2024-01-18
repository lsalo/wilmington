#
# Generate fluid properties (PVTO, PVDG, PVTW) for Wilmington field
# Reservoir simulation using a black oil model in MRST
# We use ECLIPSE "FIELD" units for fluid propeties
# Water is a brine with 30,000ppm salinity (French & Burchfield, 1990)
# Gas is a mixture of 95% methane and 5% CO2 (Sankur et al., 1986)
# Oil averages 16 deg. API (Nelson, 1976)
#

# Startup
save_data = false

# Surface densities (ECLIPSE standard conditions)
kgm3_to_lbft3 = 0.062428
temp_std_C = 15.56
p_std_bar = 1.01325
api = 18;
density_gas_kgm3 = 0.738                                    # AGA8 equation for natural gas dens. (https://www.unitrove.com/engineering/tools/gas/natural-gas-density)
density_brine_kgm3 = 1018.46                                # Batzle & Wang (1992)
density_oil_kgm3 = 999.15*(141.5/(api + 131.5))             # API gravity number definition. pure H20 density = 999.15 kg/m3
density = [density_oil_kgm3, density_brine_kgm3, density_gas_kgm3]*kgm3_to_lbft3

# PVTW


# PVDG (dry gas)

# PVTO (live oil)

# Save output to text file
if save_data
    # TBD
end