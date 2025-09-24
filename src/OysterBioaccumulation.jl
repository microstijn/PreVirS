
module OysterBioaccumulation

using DataFrames
using ..oyster
using ..VirusFateModel 

export simulate_oyster_concentration

function calculate_pseudofeces_fraction(env_conds::EnvironmentalConditions, oyster_params::OysterParameters)
    # Extract TSS and convert from kg/m³ to mg/L
    tss_mg_per_l = env_conds.tss * 1000.0

    if tss_mg_per_l <= oyster_params.tss_rejection_threshold
        return 0.0 # No rejection at low TSS
    elseif tss_mg_per_l >= oyster_params.tss_clogging_threshold
        return 1.0 # Maximum rejection at high TSS
    else
        # Linear ramp between the thresholds
        return (tss_mg_per_l - oyster_params.tss_rejection_threshold) /
               (oyster_params.tss_clogging_threshold - oyster_params.tss_rejection_threshold)
    end
end

function calculate_depuration_rate(env_conds::EnvironmentalConditions, oyster_params::OysterParameters)
    # Extract temperature from the struct
    temp_c = env_conds.temperature
    return oyster_params.k_dep_20 * oyster_params.theta_dep^(temp_c - 20.0)
end

function simulate_oyster_concentration(
        oyster_params::OysterParameters, 
        env_df::DataFrame, 
        water_results::WaterSimulationResult
    )
    c_oyster = 0.0
    oyster_results = [c_oyster]

    for i in 1:(nrow(env_df) - 1)
        env_row = env_df[i, :]
        dt_days = (env_df[i+1, :Hour] - env_row.Hour) / 24.0
        
        c_dissolved = water_results.dissolved_conc[i]
        c_sorbed = water_results.sorbed_conc[i]
        
        current_env = EnvironmentalConditions(
            env_row.Temperature_C, env_row.Salinity_ppt, 25.0,
            env_row.TSS_kg_m3, env_row.Surface_UVB_W_m2, 1.5, 5.0
        )
        
        fr_l_hr = oyster.calculate_filtration_rate(oyster_params, current_env.temperature, current_env.salinity, current_env.tss * 1000.0)
        fr_m3_day = fr_l_hr * 24.0 / 1000.0

        f_pseudo = calculate_pseudofeces_fraction(current_env, oyster_params)
        depuration_rate = calculate_depuration_rate(current_env, oyster_params)
        
        uptake_rate_free = fr_m3_day * c_dissolved
        uptake_rate_sorbed = fr_m3_day * c_sorbed * (1.0 - f_pseudo)
        total_uptake_rate = uptake_rate_free + uptake_rate_sorbed
        
        dC_oyster_dt = (total_uptake_rate / oyster_params.dry_weight_g) - (depuration_rate * c_oyster)
        c_oyster += dC_oyster_dt * dt_days
        c_oyster = max(0.0, c_oyster)
        
        push!(oyster_results, c_oyster)
    end
    
    return oyster_results
end

end # module