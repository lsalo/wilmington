module ReservoirFluidsPVT

export convert_volume_units, gasPVTBW1992, oilPVTBW1992

function convert_volume_units(values, units_in, units_out)
    """
    Convert units of a value according to the specified units.
    INPUTS:
    value:  scalar, vector, or matrix set of values
    units:  string:
            "MSCF" for gas volume in thousand standard cubic feet
            "BBL" for oil/water volume in barrels
            "M3" for volume in cubic meters
                
    OUTPUTS:
    converted_values: values in the specified units
    """
    if units_in == "BBL" && units_out == "M3"
        converted_values = values .* 0.1589873 
    elseif units_in == "M3" && units_out == "BBL"
        converted_values = values .* 6.2898107704 
    elseif units_in == "MSCF" && units_out == "M3"
        converted_values = values .* 28.3168466
    elseif units_in == "M3" && units_out == "MSCF"
        converted_values = values .* 0.0353146667215 
    end

    return converted_values
end

function gasPVTBW1992(gas_density_std_kgm3, temperature_C, pressure_bar, units)
    """
    Output gas density, FVF, and viscosity for a given temperature and pressure,
    according to Batzle & Wang, Geophysics (1992), Eq. 10 and 13.
    p_pr, t_pr are computed according to Thomas et al., JPT (1970)

    INPUTS:
    gas_density_std_kgm3:   scalar value [kg/m3]
    temperature_C:          scalar [ºC]
    pressure_bar:           scalar or vector with p values [bar]
    units:                  string, "FIELD" or "METRIC" (according to ECLIPSE)

    OUTPUTS:
    gas_dens_kgm3:          vector with density for each pressure [kg/m3]
    gas_FVF:                vector with formation volume factors for each pressure [units]
    gas_visc_cP:            vector with dynamic viscosity for each pressure [cP]
    """
    r_ct    = 83.1446261815324             # cm3⋅bar⋅K^(−1)⋅mol^(−1)
    t_abs   = 273.15 + temperature_C
    p_MPa   = pressure_bar*0.1
    air_density_std_kgm3 = 1.225
    spec_grav   = gas_density_std_kgm3/air_density_std_kgm3
    p_pr        = p_MPa ./ (4.892 - 0.4048*spec_grav)
    t_pr        = t_abs / (94.72 + 170.75*spec_grav)
    
    # Density and FVF
    e_fac_exp   = (0.45 + 8*(0.56 - 1/t_pr)^2)*p_pr.^1.2/t_pr
    e_fac       = 0.109*(3.85 - t_pr)^2 * exp.(-e_fac_exp)
    z_fac       = (0.03 + 0.00527*(3.5 -t_pr)^3)*p_pr .+ (0.642*t_pr - 0.007*t_pr^4 - 0.52) .+ e_fac    # gas compressibility factor
    gas_dens_kgm3 = 1000*28.8*spec_grav*pressure_bar ./ (z_fac*r_ct*t_abs)       
    if units == "METRIC"    # rm3/sm3
        gas_FVF         = gas_density_std_kgm3 ./ gas_dens_kgm3
    elseif units == "FIELD" # rb/Mscf
        kgm3_to_kgrb    = 1/6.2898107704
        kgm3_to_kgMscf  = 1/0.0353147
        gas_FVF         = (gas_density_std_kgm3*kgm3_to_kgMscf) ./ (gas_dens_kgm3*kgm3_to_kgrb)
    end

    # Viscosity
    gas_visc_atm_cP = 1e-4*(t_pr*(28 + 48*spec_grav - 5*spec_grav^2) - 6.47/spec_grav^2 + 35/spec_grav + 1.14*spec_grav - 15.55)
    f1 = (1057 - 8.08*t_pr) ./ p_pr
    f2 = (796*p_pr.^0.5 .- 704) ./ ((t_pr - 1)^0.7 * (p_pr .+ 1))
    gas_visc_cP = gas_visc_atm_cP*1e-3*p_pr.*(f1 .+ f2 .- 3.24*t_pr .- 38)

    return gas_dens_kgm3, gas_FVF, gas_visc_cP

end

function oilPVTBW1992(gas_density_std_kgm3, api, temp_C, pressure_bar, units)
    """
    Output oil Rs, density, and FVF, for a given temperature and pressure,
    according to Batzle & Wang, Geophysics (1992), Eq. 24 and 26. In this paper:
    - Rs, FVF, density are computed according to Standing (1962)
    - Dead oil viscosity is computed from relationships in Beggs & Robinson, JPT (1975) and Beal, Petr. Trans. (1946)
    We then calculate live oil viscosity following Beggs & Robinson, JPT (1975), which uses the equation of 
    Chew & Connally, Trans AIME (1959) with updated parameters.

    INPUTS:
    gas_density_std_kgm3:   scalar value [kg/m3]
    api:                    scalar value [º api]
    temperature_C:          scalar [ºC]
    pressure_bar:           scalar or vector with p values [bar]
    units:                  string, "FIELD" or "METRIC" (according to ECLIPSE)

    OUTPUTS:
    oil_Rs:                 vector with gas-oil solution ratio for each pressure [units]
    oil_B:                  vector with formation volume factors for each pressure [units]
    oil_dens_kgm3:          vector with density for each pressure [kg/m3]
    oil_visc_cP:            vector with dynamic viscosity for each pressure [cP]
    """
    p_MPa                   = pressure_bar*0.1
    air_density_std_kgm3    = 1.225
    spec_grav               = gas_density_std_kgm3/air_density_std_kgm3
    oil_dens_std_gcm3       = 0.99915*(141.5/(api + 131.5)) 
    scf_in_stb              = 5.6145833333333    
    t_F                     = temp_C*1.8 + 32  
    bar_to_psi              = 14.503773773

    # Rs (Eq. 21)
    oil_Rs = 2.03*spec_grav*(p_MPa*exp(0.02878*api - 0.00377*temp_C)).^1.205                            # L/L
    oil_Rs_scfSTB = oil_Rs*scf_in_stb

    # B0 (Eq. 23)
    oil_B0 = 0.972 .+ 3.8e-4*(2.4*oil_Rs*(spec_grav/oil_dens_std_gcm3)^0.5 .+ temp_C .+ 17.8).^1.175     # V / V (any same units)

    # density (Eq. 24) and FVF
    oil_dens_gcm3 = (oil_dens_std_gcm3 .+ 1.2e-3*spec_grav*oil_Rs) ./ oil_B0
    
    # Viscosity (Beggs & Robinson, JPT, 1975; Beal, 1949)
    z       = 3.0324 - 0.02023*api
    x       = 10^z*t_F^-1.163
    mu_t    = 10^x - 1                      # cP, dead-oil, p ignored
    i       = 10^(18.6*(0.1*log10(mu_t) + (log10(mu_t) + 2)^-0.1 - 0.985))
    mu_od   = mu_t .+ 0.145*p_MPa*i         # cP, dead-oil, with p
    a       = 10.715*(oil_Rs_scfSTB .+ 100).^-0.515
    b       = 5.44*(oil_Rs_scfSTB .+ 150).^-0.338
    oil_visc_cP = a.*mu_od.^b               # cP, live-oil, with p

    # Add undersaturated values for highest Rs 
    # For density, we use Eq. 18 to determine the density increase from p_max to p_add, and add it to live oil density
    p_max   = p_MPa[end]
    p_add   = p_max .+ [0, 5, 10]
    rho_max = oil_dens_gcm3[end]
    rho_p   = oil_dens_std_gcm3 .+ (2.77e-3*p_add .- 1.71e-7*p_add.^3)*(oil_dens_std_gcm3 - 1.15)^2 .+ 3.49e-4*p_add
    rho_p_add = rho_p[2:3] .- rho_p[1]
    oil_dens_gcm3 = vcat(oil_dens_gcm3, oil_dens_gcm3[end] .+ rho_p_add)
    # For viscosity, we use Beal's correlation to compute mu_od, then we multiply live oil visc by the same factor
    mu_od_p = mu_t .+ 0.145*p_add*i 
    mu_f = mu_od_p[2:3] ./ mu_od_p[1]
    oil_visc_cP = vcat(oil_visc_cP, oil_visc_cP[end].*mu_f)

    # FVF and output
    oil_FVF = oil_dens_std_gcm3 ./ oil_dens_gcm3
    oil_dens_kgm3 = oil_dens_gcm3*1e3
    if units == "METRIC"    
        oil_Rs  = vcat(oil_Rs, oil_Rs[end], oil_Rs[end])
        p_bar   = vcat(pressure_bar, p_add[2:3]*10)
        return oil_dens_kgm3, oil_Rs, p_bar, oil_FVF, oil_visc_cP
    elseif units == "FIELD"
        oil_Rs_MscfSTB  = vcat(oil_Rs_scfSTB, oil_Rs_scfSTB[end], oil_Rs_scfSTB[end])/1000
        p_psi           = vcat(pressure_bar, p_add[2:3]*10)*bar_to_psi
        return oil_dens_kgm3, oil_Rs_MscfSTB, p_psi, oil_FVF, oil_visc_cP
    end

end


end