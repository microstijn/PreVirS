# Inital version of a few functions that can calculate excretion rate of oysters. 
# see https://github.com/microstijn/PreVirS/ for a description.
module oyster

export calculate_filtration_rate

"""
    _f_temp(temperature_c)

Calculate the temperature limitation factor `f(T)`.
This is a Gaussian-like function centered at 27°C.
"""
function _f_temp(temperature_c::Real)
    return exp(-0.006 * (temperature_c - 27.0)^2)
end

"""
    _f_salinity(salinity_psu)

Calculate the salinity limitation factor `f(S)`.
This is a piecewise function representing the oyster's tolerance to salinity.
"""
function _f_salinity(salinity_psu::Real)
    if salinity_psu < 5.0
        return 0.0
    elseif 5.0 <= salinity_psu <= 12.0
        return 0.0926 * (salinity_psu - 0.0139)
    else # salinity_psu > 12.0
        return 1.0
    end
end

"""
    _f_tss(tss_mg_per_l)

Calculate the Total Suspended Solids limitation factor `f(TSS)`.
This piecewise function models the effect of particle concentration on filtration.
"""
function _f_tss(tss_mg_per_l::Real)
    if tss_mg_per_l < 4.0
        return 0.1
    elseif 4.0 <= tss_mg_per_l <= 25.0
        return 1.0
    else # tss_mg_per_l > 25.0
        return 10.364 * log(tss_mg_per_l)^(-2.0477)
    end
end

# the main function below. 

"""
    calculate_filtration_rate(dry_weight_g, temperature_c, salinity_psu, tss_mg_per_l)

Calculates the filtration rate (FR) of an eastern oyster (*Crassostrea virginica*).

The model is based on the general formula:
`FR = 0.17 * W_dw^0.75 * f(T) * f(S) * f(TSS)`

# Arguments
- `dry_weight_g::Real`: Dry weight of the oyster tissue in grams (g).
- `temperature_c::Real`: Water temperature in Celsius (°C).
- `salinity_psu::Real`: Salinity in practical salinity units (psu), equivalent to ppt.
- `tss_mg_per_l::Real`: Total suspended solids in milligrams per liter (mg/L).

# Returns
- `Real`: The calculated filtration rate in liters per hour per oyster (L·hr⁻¹).

# Examples
```jldoctest
julia> calculate_filtration_rate(0.5, 25.0, 15.0, 10.0)
0.10123187284893093
```
"""
function calculate_filtration_rate(
    dry_weight_g::Real,
    temperature_c::Real,
    salinity_psu::Real,
    tss_mg_per_l::Real
)
    # Allometric scaling factor for oyster size
    allometric_scaling = dry_weight_g^0.75

    # Calculate the final filtration rate by combining all factors
    fr = 0.17 *
         allometric_scaling *
         _f_temp(temperature_c) *
         _f_salinity(salinity_psu) *
         _f_tss(tss_mg_per_l)

    return fr
end

"""
    generate_gridded_dataset(steps_per_variable::Int)

Generates a gridded dataset by creating a mesh over the realistic ranges
of all four input variables. This is ideal for systematic model exploration.

# Arguments
- `steps_per_variable::Int`: The number of steps for each variable's range. The total
  number of points will be `steps_per_variable ^ 4`.

# Returns
- `Vector{Vector{Float64}}`: A vector of vectors, where each inner vector is a data point.
- `Vector{String}`: The header for the data columns.
"""
function generate_gridded_dataset(steps_per_variable::Int)
    total_points = steps_per_variable^4
    println("Generating gridded dataset with $steps_per_variable steps per variable ($total_points total points)...")

    # Define realistic ranges for the input variables
    weight_range = (0.1, 5.0)   # g
    temp_range = (5.0, 32.0)    # °C
    salinity_range = (3.0, 30.0) # psu
    tss_range = (2.0, 100.0)   # mg/L

    # Create linearly spaced points for each variable's range
    weights = LinRange(weight_range[1], weight_range[2], steps_per_variable)
    temps   = LinRange(temp_range[1], temp_range[2], steps_per_variable)
    sals    = LinRange(salinity_range[1], salinity_range[2], steps_per_variable)
    tsss    = LinRange(tss_range[1], tss_range[2], steps_per_variable)

    # Define the header
    header = ["dry_weight_g", "temperature_c", "salinity_psu", "tss_mg_per_l", "filtration_rate_l_hr"]

    # Initialize the dataset
    dataset = Vector{Vector{Float64}}()
    sizehint!(dataset, total_points) # Pre-allocate memory for efficiency

    # Iterate through every combination of the variables to create the grid
    for w in weights
        for t in temps
            for s in sals
                for ts in tsss
                    # Calculate the filtration rate for this point in the grid
                    fr = calculate_filtration_rate(w, t, s, ts)

                    data_row = [round(w, digits=4),
                                round(t, digits=2),
                                round(s, digits=2),
                                round(ts, digits=2),
                                round(fr, digits=6)]

                    push!(dataset, data_row)
                end
            end
        end
    end

    println("Successfully generated gridded dataset.")
    return dataset, header
end

