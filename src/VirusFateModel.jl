# VirusFateModel.jl

module VirusFateModel

# Dependencies 
using DifferentialEquations
using DataFrames
using ..oyster

# Exports
export VirusParameters
export EnvironmentalConditions
export WaterSimulationResult
export run_simulation
export simulate_water_dynamics

# -structs

"""
    VirusParameters

Holds all kinetic and physical parameters for a specific virus type. 
This struct defines the intrinsic properties of the virus and its interactions.

# Fields
- `k20::Float64`: Reference decay rate at 20°C. Unit: \$day^{-1}\$.
- `theta::Float64`: Dimensionless temperature adjustment coefficient for the Arrhenius model. Unit: dimensionless.
- `alpha::Float64`: Dimensionless calibration coefficient for the effect of salinity on decay. Unit: dimensionless.
- `k_I::Float64`: Light-dependent decay coefficient for UV inactivation. Unit: \$(W/m^2)^{-1} \\cdot day^{-1}\$.
- `adsorption_rate::Float64`: Second-order rate constant for virus adsorption to suspended solids. Unit: \$m^3/(kg \\cdot day)\$.
- `desorption_rate::Float64`: First-order rate constant for virus desorption from suspended solids. Unit: \$day^{-1}\$.
- `settling_velocity::Float64`: Gravitational settling speed of virus particles (especially relevant for vesicles). Unit: \$m/day\$.
- `sorbed_protection_factor::Float64`: Dimensionless factor (0-1) that reduces the decay rate for viruses sorbed to particles. Unit: dimensionless.
"""
struct VirusParameters
    k20::Float64
    theta::Float64
    alpha::Float64
    k_I::Float64
    adsorption_rate::Float64
    desorption_rate::Float64
    settling_velocity::Float64
    sorbed_protection_factor::Float64
end

"""
    EnvironmentalConditions

Holds all time-variable environmental parameters for the aquatic system.
This struct defines the state of the environment at a specific point in time.

# Fields
- `temperature::Float64`: Ambient water temperature. Unit: °C.
- `salinity::Float64`: Ambient water salinity. Unit: ppt (parts per thousand).
- `salinity_ref::Float64`: Reference salinity at which decay rates were measured. Unit: ppt.
- `tss::Float64`: Concentration of total suspended solids. Unit: \$kg/m^3\$.
- `surface_uvb::Float64`: Intensity of UVB radiation at the water surface. Unit: \$W/m^2\$.
- `light_extinction_coeff::Float64`: Coefficient for light attenuation in the water column. Unit: \$m^{-1}\$.
- `water_depth::Float64`: Total height of the water column. Unit: \$m\$.
"""
struct EnvironmentalConditions
    temperature::Float64
    salinity::Float64
    salinity_ref::Float64
    tss::Float64
    surface_uvb::Float64
    light_extinction_coeff::Float64
    water_depth::Float64
end

"""
    WaterSimulationResult

Holds the time-series results from a virus in water simulation.

# Fields
- `time::Vector{Float64}`: Time points of the simulation (hours).
- `dissolved_conc::Vector{Float64}`: Concentration of dissolved virus (vg/m³).
- `sorbed_conc::Vector{Float64}`: Concentration of virus sorbed to suspended solids (vg/m³).
"""
struct WaterSimulationResult
    time::Vector{Float64}
    dissolved_conc::Vector{Float64}
    sorbed_conc::Vector{Float64}
    settled_conc::Vector{Float64}
end

# internals

function calculate_temp_dependent_decay(virus::VirusParameters, env::EnvironmentalConditions)::Float64
    return virus.k20 * virus.theta^(env.temperature - 20.0)
end

function calculate_salinity_modifier(virus::VirusParameters, env::EnvironmentalConditions)::Float64
    return 1.0 - virus.alpha * (env.salinity - env.salinity_ref)
end

function calculate_k_uv_avg(virus::VirusParameters, env::EnvironmentalConditions)::Float64
    I_0 = env.surface_uvb
    K_ext = env.light_extinction_coeff
    H = env.water_depth

    if K_ext <= 0.0 || H <= 0.0
        return virus.k_I * I_0
    end
    
    avg_uvb_intensity = I_0 * (1.0 - exp(-K_ext * H)) / (K_ext * H)
    return virus.k_I * avg_uvb_intensity
end

# internal ode

function virus_dynamics_water!(du, u, p, t)
    # The state vector `u` now has three components
    # u = [C_dissolved(vg/m³), C_sorbed_TSS(vg/m³), C_settled_cumulative(vg/m²)]
    C_dissolved, C_sorbed_tss, C_settled_cumulative = u
    virus_params, env_conds = p

    # Calculate total decay rate in water
    k_background = virus_params.k20 * virus_params.theta^(env_conds.temperature - 20.0)
    k_total_water = k_background + virus_params.k_I * env_conds.surface_uvb

    # Fluxes
    adsorption_flux = virus_params.adsorption_rate * C_dissolved * env_conds.tss
    desorption_flux = virus_params.desorption_rate * C_sorbed_tss
    settling_rate_loss = virus_params.settling_velocity / env_conds.water_depth
    decay_sorbed = k_total_water * virus_params.sorbed_protection_factor

    # Calculate the flux of viruses settling to the bottom
    # This is the amount of virus (vg) settling per square meter per day.
    settling_flux_to_bed = virus_params.settling_velocity * C_sorbed_tss 

    # Equations of change
    # 1. Dissolved concentration
    du[1] = -k_total_water * C_dissolved - adsorption_flux + desorption_flux
    # 2. Sorbed concentration
    du[2] = -decay_sorbed * C_sorbed_tss + adsorption_flux - desorption_flux - settling_rate_loss * C_sorbed_tss
    # 3. The accumulator for settled viruses
    du[3] = settling_flux_to_bed
end

# public 

"""
    run_water_simulation(virus_params, env_conditions, u0, t_span)

Runs a virus fate and transport simulation.

# Arguments
- `virus_params::VirusParameters`: The parameters for the virus type.
- `env_conditions::EnvironmentalConditions`: The environmental conditions.
- `u0::Vector{Float64}`: The initial concentrations `[dissolved, sorbed]`.
- `t_span::Tuple{Float64, Float64}`: The simulation time span `(start, end)`.

# Returns
- An `ODESolution` object containing the simulation results.
"""
function run_water_simulation(virus_params, env_conditions, u0, t_span)
    params = (virus_params, env_conditions)
    prob = ODEProblem(virus_dynamics_water!, u0, t_span, params)
    sol = solve(prob, Tsit5(), save_everystep=false, save_start=false, save_end=true)
    return sol
end

"""
    simulate_water_dynamics(virus_params, env_df)

Runs a 24-hour simulation, now returning time, total concentration, and sorbed concentration.
"""
function simulate_water_dynamics(
        virus_params::VirusParameters, 
        env_df::DataFrame,
        influx_df::DataFrame,
        initial_concentrations::Vector{Float64}
    )
    current_concentrations = copy(initial_concentrations)

    # Arrays to store results for all three pools
    time_results = [0.0]
    dissolved_results = [current_concentrations[1]] # Dissolved + Sorbed TSS
    sorbed_results = [current_concentrations[2]] # Sorbed
    settled_results = [current_concentrations[3]] # settled

    for i in 1:(nrow(env_df) - 1)
        current_concentrations[1] += influx_df.Influx_vg_m3[i]
        row = env_df[i, :]
        t_start, t_end = Float64(row.Hour), Float64(env_df[i+1, :Hour])
        
        current_env = EnvironmentalConditions(
            row.Temperature_C, 
            row.Salinity_ppt,
            25.0,
            row.TSS_kg_m3,
            row.Surface_UVB_W_m2,
            1.5,
            5.0
        )
        
        solution_step = run_water_simulation(virus_params, current_env, current_concentrations, (t_start, t_end))
        current_concentrations = solution_step.u[1]
        
        push!(time_results, t_end)
        push!(dissolved_results, current_concentrations[1])
        push!(sorbed_results, current_concentrations[2])
        push!(settled_results, current_concentrations[3]) 
    end
    
    return WaterSimulationResult(time_results, dissolved_results, sorbed_results, settled_results)

end

end # module VirusFateModel
