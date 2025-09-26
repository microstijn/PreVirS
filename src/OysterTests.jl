# Filename: PreVirS/test/test_simulations.jl

module OysterTests

using Accessors # Use the modern library for modifying structs
using Test

using ..oyster
using ..oysterData
using ..VirusFateModel
using ..OysterBioaccumulation

export run_simulation_tests

"""
    run_simulation_tests()

Runs an extensive suite of sensitivity tests for the PreVirS simulation logic.
"""
function run_simulation_tests()
    @testset "Simulation Logic Sensitivity Tests" begin
        # (Setup is unchanged)
        duration_days_short = 1; env_df_short = generate_synthetic_data(duration_days_short)
        influx_df_short = generate_virus_influx(duration_days_short, [1e4], [0], [24], 0.0)
        p_virus_base = VirusParameters(0.2, 1.07, 0.01, 0.05, 1.0, 0.2, 0.05, 0.5)
        p_oyster_base = OysterParameters(1.0, 0.1, 1.05, 100.0, 200.0, 0.1, 0.8)
        initial_cond = [0.0, 0.0, 0.0]
        res_base_water = simulate_water_dynamics(p_virus_base, env_df_short, influx_df_short, initial_cond)
        res_base_oyster = simulate_oyster_concentration(p_oyster_base, env_df_short, res_base_water)
        final_total_base = last(res_base_water.dissolved_conc) + last(res_base_water.sorbed_conc)
        final_sorbed_base = last(res_base_water.sorbed_conc)
        final_settled_base = last(res_base_water.settled_conc)
        final_oyster_base = last(res_base_oyster)

        @testset "VirusParameter Sensitivities" begin
            # (Tests for k20, k_I, adsorption, desorption, settling, protection are unchanged)
            p_mod_k20 = @set p_virus_base.k20 = 0.8; res_mod_k20 = simulate_water_dynamics(p_mod_k20, env_df_short, influx_df_short, initial_cond)
            @test last(res_mod_k20.dissolved_conc) + last(res_mod_k20.sorbed_conc) < final_total_base

            # Test theta: Higher temp sensitivity + COLD day (<20C) -> MORE virus (slower decay)
            p_mod = @set p_virus_base.theta = 1.1
            res_mod = simulate_water_dynamics(p_mod, env_df_short, influx_df_short, initial_cond)
            @test last(res_mod.dissolved_conc) + last(res_mod.sorbed_conc) > final_total_base
            
            # Test alpha: Higher salinity sensitivity -> LESS virus
            p_mod = @set p_virus_base.alpha = 0.05
            res_mod = simulate_water_dynamics(p_mod, env_df_short, influx_df_short, initial_cond)
            @test last(res_mod.dissolved_conc) + last(res_mod.sorbed_conc) < final_total_base
            
             p_mod_kI = @set p_virus_base.k_I = 0.2; res_mod_kI = simulate_water_dynamics(p_mod_kI, env_df_short, influx_df_short, initial_cond); @test last(res_mod_kI.dissolved_conc) + last(res_mod_kI.sorbed_conc) < final_total_base
             p_mod_ads = @set p_virus_base.adsorption_rate = 5.0; res_mod_ads = simulate_water_dynamics(p_mod_ads, env_df_short, influx_df_short, initial_cond); @test last(res_mod_ads.sorbed_conc) > final_sorbed_base
             p_mod_des = @set p_virus_base.desorption_rate = 0.8; res_mod_des = simulate_water_dynamics(p_mod_des, env_df_short, influx_df_short, initial_cond); @test last(res_mod_des.sorbed_conc) < final_sorbed_base
             p_mod_set = @set p_virus_base.settling_velocity = 0.2; res_mod_set = simulate_water_dynamics(p_mod_set, env_df_short, influx_df_short, initial_cond); @test last(res_mod_set.settled_conc) > final_settled_base
             p_mod_pro = @set p_virus_base.sorbed_protection_factor = 0.9; res_mod_pro = simulate_water_dynamics(p_mod_pro, env_df_short, influx_df_short, initial_cond); @test last(res_mod_pro.sorbed_conc) > final_sorbed_base
        end

        @testset "OysterParameter Sensitivities" begin
            # (Tests for dry_weight, k_dep_20, efficiencies are unchanged)
             p_mod_dw = @set p_oyster_base.dry_weight_g = 5.0; res_mod_dw = simulate_oyster_concentration(p_mod_dw, env_df_short, res_base_water); @test last(res_mod_dw) < final_oyster_base
             p_mod_kdep = @set p_oyster_base.k_dep_20 = 0.5; res_mod_kdep = simulate_oyster_concentration(p_mod_kdep, env_df_short, res_base_water); @test last(res_mod_kdep) < final_oyster_base

            # --- CORRECTED TEST LOGIC for theta_dep ---
            # Test theta_dep: Higher temp sensitivity + COLD day (<20C) -> MORE virus in oyster (slower depuration)
            p_mod = @set p_oyster_base.theta_dep = 1.1
            res_mod = simulate_oyster_concentration(p_mod, env_df_short, res_base_water)
            @test last(res_mod) > final_oyster_base
            
            # Test tss_rejection_threshold: Lower threshold -> LESS virus in oyster (rejects earlier)
            p_mod = @set p_oyster_base.tss_rejection_threshold = 5.0 # Lower than baseline TSS
            res_mod = simulate_oyster_concentration(p_mod, env_df_short, res_base_water)
            @test last(res_mod) < final_oyster_base
            -
            # Test tss_clogging_threshold: Lower threshold -> LESS virus in oyster (clogs earlier)
            # Must also lower rejection_threshold to activate the mechanism
            p_mod = @set p_oyster_base.tss_rejection_threshold = 5.0
            p_mod = @set p_mod.tss_clogging_threshold = 10.0
            res_mod = simulate_oyster_concentration(p_mod, env_df_short, res_base_water)
            @test last(res_mod) < final_oyster_base
            
             p_mod_eff = @set p_oyster_base.efficiency_free = 0.5; res_mod_eff = simulate_oyster_concentration(p_mod_eff, env_df_short, res_base_water); @test last(res_mod_eff) > final_oyster_base
             p_mod_efs = @set p_oyster_base.efficiency_sorbed = 0.95; res_mod_efs = simulate_oyster_concentration(p_mod_efs, env_df_short, res_base_water); @test last(res_mod_efs) > final_oyster_base
        end
        # (Sewer Overflow test is unchanged)
        @testset "Sewer Overflow Event Impact" begin
            duration_days_long = 7
            env_df_long = generate_synthetic_data(duration_days_long)
            influx_df_long = generate_virus_influx(duration_days_long, [1e6], [30], [6], 1e2)
            water_results = simulate_water_dynamics(p_virus_base, env_df_long, influx_df_long, initial_cond)
            oyster_results = simulate_oyster_concentration(p_oyster_base, env_df_long, water_results)
            pre_cso_oyster_conc = oyster_results[31]
            peak_oyster_conc = maximum(oyster_results)
            peak_time_index = findmax(oyster_results)[2]
            peak_hour = water_results.time[peak_time_index]
            @test peak_oyster_conc > 5 * pre_cso_oyster_conc
            @test peak_hour >= 30
        end
    end
end

end