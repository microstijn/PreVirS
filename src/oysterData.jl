module oysterData

using UnicodePlots
using Dates
using DataFrames
using ..oyster

export simulate_yearly_environment
export generate_synthetic_data
export generate_virus_influx


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
    generate_synthetic_data(duration_days::Int)

Creates a synthetic `DataFrame` for environmental conditions over a specified duration.
Daily cycles for temperature and UVB are repeated each day.
"""
function generate_synthetic_data(duration_days::Int)
    hours = 0.0:1.0:(duration_days * 24)

    # (Model parameters are the same as before)
    salinity_base = 18.0; salinity_amplitude = 7.0
    temp_base = 15.0; temp_amplitude = 1.5; temp_peak_hour = 16.0
    tss_base = 8.0; tss_amplitude = 7.0
    uvb_peak = 25.0; sunrise = 6.0; sunset = 20.0
    daylight_hours = sunset - sunrise

    # Generate Time-Series Data, using modulo (%) for daily cycles
    salinity = [salinity_base + salinity_amplitude * cos(2 * π * h / 12.4) for h in hours]
    # Use h % 24 to repeat the daily temperature cycle
    temperature = [temp_base - temp_amplitude * cos(2 * π * ((h % 24) - (temp_peak_hour - 12)) / 24) for h in hours]
    tss_mg_L = [tss_base + tss_amplitude * (sin(2 * π * h / 12.4))^2 for h in hours]
    tss_kg_m3 = tss_mg_L / 1000.0
    # Use h % 24 to repeat the daily solar cycle
    uvb = [(h % 24) > sunrise && (h % 24) < sunset ? uvb_peak * sin(π * ((h % 24) - sunrise) / daylight_hours) : 0.0 for h in hours]

    env_data = DataFrame(
        Hour = hours,
        Temperature_C = temperature,
        Salinity_ppt = salinity,
        TSS_kg_m3 = tss_kg_m3,
        Surface_UVB_W_m2 = uvb
    )
    return env_data
end

"""
    generate_virus_influx(duration_days::Int)

Generates a synthetic time-series of virus influx.
Includes a low baseline and a high-concentration sewer overflow (CSO) event
on the second day.
"""
function generate_virus_influx(duration_days::Int)
    hours = 0.0:1.0:(duration_days * 24)

    # --- Influx Parameters ---
    baseline_influx = 1e2  # Low, constant background level (vg/m³)
    cso_influx = 1e6       # High concentration during overflow (vg/m³)
    cso_start_hour = 30    # 6 AM on Day 2
    cso_duration_hours = 6

    # Initialize with baseline
    influx = fill(baseline_influx, length(hours))

    # Add the CSO spike
    cso_end_hour = cso_start_hour + cso_duration_hours
    for i in eachindex(hours)
        if hours[i] >= cso_start_hour && hours[i] < cso_end_hour
            influx[i] = cso_influx
        end
    end
    
    return DataFrame(Hour = hours, Influx_vg_m3 = influx)
end


end