using Pkg
#Pkg.develop(path=joinpath(@__DIR__, ".."))
Pkg.activate(joinpath(@__DIR__, ".."))

using Revise
using PreVirS
using UnicodePlots

# set length of simulation data
simulation_days = 50

# genereate sim data
env_df = generate_synthetic_data(simulation_days)

baseline_influx = 1e2  # Low, constant background level (vg/m³)
cso_influx = [1e8, 1e8, 1e8]       # High concentration during overflow (vg/m³)
cso_start_hour = [30, 200, 1000]    # 6 AM on Day 2
cso_duration_hours = [12, 12, 24]

influx_df = generate_virus_influx(simulation_days, cso_influx, cso_start_hour, cso_duration_hours, baseline_influx)

# set virus and oyster params
virus_params = VirusParameters(
    0.23,       # k20: Reference decay rate at 20°C
    1.076,      # theta: Temperature adjustment coefficient for decay
    0.0,        # alpha: Salinity effect coefficient
    0.05,       # k_I: UV light decay coefficient
    1.0,        # adsorption_rate: Rate of attachment to particles
    0.2,        # desorption_rate: Rate of detachment from particles
    0.05,       # settling_velocity: Sinking speed of sorbed virus [m/day] 
    0.5         # sorbed_protection_factor: (50%) reduction in decay when sorbed 
)

oyster_params = OysterParameters(
    1.0,        # dry_weight_g: Dry weight of the oyster in grams 
    0.107,      # k_dep_20: Depuration (elimination) rate at 20°C 
    1.055,      # theta_dep: Temperature coefficient for depuration 
    100.0,      # tss_rejection_threshold: TSS level to start rejecting particles [mg/L]
    200.0,      # tss_clogging_threshold: TSS level for maximum particle rejection [mg/L]
    0.5,        # efficiency_free: (50%) uptake efficiency for free viruses
    1           # efficiency_sorbed: (100%) uptake efficiency for sorbed viruses 
)

# run sim on sim data
water_sim = simulate_water_dynamics(virus_params, env_df, influx_df, [100, 0.0, 0.0]);
oyster_sim = simulate_oyster_concentration(oyster_params, env_df, water_sim)

# plotting 
begin
    p_water = lineplot(
        water_sim.time,
        water_sim.dissolved_conc,
        name="Virus (dissolved)",
        title="Virus in Water Column",
        xlabel="Time (hours)",
        ylabel="Conc (vg/m³)",
        yscale=:log10,
        width=300,
        height=30,
        ylim = (10^0, 10^8),
        canvas = BrailleCanvas
    )

    lineplot!(p_water, water_sim.time, water_sim.sorbed_conc; name = "Virus (sorbed)")
    #lineplot!(p_water, water_sim.time, water_sim.dissolved_conc; name = "Virus (free)")
    lineplot!(p_water, water_sim.time, water_sim.settled_conc; name = "Virus (settled)")

    lineplot!(
        p_water,
        water_sim.time,
        oyster_sim;
        name ="Virus in oyster"
    )

    lineplot!(
        p_water,
        water_sim.time,
        repeat([1000], size(water_sim.time, 1)),
        name ="Virus limit (1000?)"
    )

end

    lineplot(
        water_sim.time,
        oyster_sim;
        name ="Virus in oyster",
        ylabel="Conc (vg/g)",
        width=300,
        height=10
    )


run_simulation_tests();