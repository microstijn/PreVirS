using Pkg
#Pkg.develop(path=joinpath(@__DIR__, ".."))
Pkg.activate(joinpath(@__DIR__, ".."))

using Revise
using PreVirS
using UnicodePlots

simulation_days = 7

env_df = generate_synthetic_data(simulation_days)
influx_df = generate_virus_influx(simulation_days)

virus_params = VirusParameters(0.23, 1.076, 0.0, 0.05, 1.0, 0.2, 0.01, 0.5)
oyster_params = OysterParameters(1.0, 0.107, 1.055, 100.0, 200.0)

water_sim = simulate_water_dynamics(virus_params, env_df, influx_df, [100, 0.0, 0.0]);

oyster_sim = simulate_oyster_concentration(oyster_params, env_df, water_sim)

p_water = lineplot(
    water_sim.time,
    water_sim.dissolved_conc .+ water_sim.sorbed_conc,
    name="Virus (dissolved + particulate)",
    title="Virus in Water Column",
    xlabel="Time (hours)",
    ylabel="Conc (vg/m³)",
    yscale=:log10,
    width=70,
    height=15
)

lineplot!(p_water, water_sim.time, water_sim.sorbed_conc; name = "Virus (sorbed)")
lineplot!(p_water, water_sim.time, water_sim.dissolved_conc; name = "Virus (free)")
#lineplot!(p_water, water_sim.time, water_sim.settled_conc; name = "Virus (settled)")

lineplot!(
    p_water,
    water_sim.time,
    oyster_sim;
    name ="Virus in oyster"
)




