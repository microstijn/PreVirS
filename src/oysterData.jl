module oysterData

using UnicodePlots
using Dates
using ..oyster

export simulate_yearly_environment
export plot_simulation_timeseries


"""
    simulate_yearly_environment(year::Int)

Generates a year of hourly environmental data with seasonal and daily cycles.

# Arguments
- `year::Int`: The year to simulate.

# Returns
- A tuple containing:
    - `timestamps::Vector{DateTime}`: Hourly timestamps for the year.
    - `temps_c::Vector{Float64}`: Simulated hourly temperature (°C).
    - `salinities_psu::Vector{Float64}`: Simulated hourly salinity (psu).
    - `tss_mg_per_l::Vector{Float64}`: Simulated hourly TSS (mg/L).
"""
function simulate_yearly_environment(year::Int)
    @info "Simulating Hourly Environmental Data for $year"
    timestamps = DateTime(year, 1, 1, 0):Hour(1):DateTime(year, 12, 31, 23)
    num_hours = length(timestamps)

    # Initialize vectors to store the data
    temps_c = zeros(num_hours)
    salinities_psu = zeros(num_hours)
    tss_mg_per_l = zeros(num_hours)

    # --- Model Parameters ---
    # Temperature (°C)
    mean_annual_temp = 15.0
    seasonal_temp_amplitude = 10.0 # Fluctuation from the mean
    daily_temp_amplitude = 1.5   # Day/night fluctuation

    # Salinity (psu)
    mean_annual_salinity = 18.0
    seasonal_salinity_amplitude = 4.0 # Lower in spring, higher in summer
    salinity_phase_shift_days = 60 # Peak salinity lags peak temp

    # TSS (mg/L)
    baseline_tss = 5.0
    storm_event_probability = 0.01 # 1% chance of a storm starting each hour
    storm_duration_hours = 6

    # --- Simulation Loop ---
    in_storm = 0
    for (i, ts) in enumerate(timestamps)
        day_of_year = Dates.dayofyear(ts)
        hour_of_day = Dates.hour(ts)

        # 1. Temperature Model (Seasonal + Daily Sine Waves)
        seasonal_temp = mean_annual_temp - seasonal_temp_amplitude * cos((2π * day_of_year) / 365.25)
        daily_temp = -daily_temp_amplitude * cos((2π * hour_of_day) / 24)
        temps_c[i] = seasonal_temp + daily_temp

        # 2. Salinity Model (Seasonal Sine Wave)
        salinities_psu[i] = mean_annual_salinity - seasonal_salinity_amplitude * cos((2π * (day_of_year - salinity_phase_shift_days)) / 365.25)

        # 3. TSS Model (Baseline + Random Storm Events)
        if in_storm > 0
            # Continue ongoing storm
            tss_mg_per_l[i] = baseline_tss + rand(20:80)
            in_storm -= 1
        elseif rand() < storm_event_probability
            # Start a new storm
            in_storm = storm_duration_hours
            tss_mg_per_l[i] = baseline_tss + rand(20:80)
        else
            # Normal conditions
            tss_mg_per_l[i] = baseline_tss + rand() * 5.0
        end
    end

    @info "Simulation complete."
    return timestamps, temps_c, salinities_psu, tss_mg_per_l
end

"""
    plot_simulation_timeseries(timestamps, temps, salinities, tsss, filtration_rates)

Plots the simulated environmental data and the resulting filtration rate over time.
"""
function plot_simulation_timeseries(timestamps, temps, salinities, tsss, filtration_rates)
    println("\n--- Plotting Yearly Simulation Results ---")

    # To make plotting faster, we can plot every 24th hour (daily average)
    indices = 1:24:length(timestamps)

    # Helper function to create plots
    function create_plot(y_data, title, y_label)
        # The `xlabels` keyword is not supported; it has been removed.
        # The x-axis will be labeled with the start and end of the index range.
        lineplot(
            indices, y_data[indices],
            title=title,
            xlabel="Hour of the Year (sampled daily)",
            ylabel=y_label,
            width=80,
            height=15
        )
    end

    p_temp = create_plot(temps, "Simulated Temperature", "Temp (°C)")
    println(p_temp)

    p_sal = create_plot(salinities, "Simulated Salinity", "Salinity (psu)")
    println(p_sal)

    p_tss = create_plot(tsss, "Simulated TSS", "TSS (mg/L)")
    println(p_tss)

    p_fr = create_plot(filtration_rates, "Calculated Filtration Rate", "FR (L/hr)")
    println(p_fr)
end

end