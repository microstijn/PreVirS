# Inital version of a few functions that can calculate excretion rate of oysters. 
# see https://github.com/microstijn/PreVirS/ for a description.
module oyster

export OysterParameters
export f_temp
export f_salinity
export f_tss
export calculate_filtration_rate


"""
    OysterParameters

Represents the core physical and biological parameters for an oyster.
"""
struct OysterParameters
    # Physical Properties
    "Dry weight of the oyster tissue in grams (g)."
    dry_weight_g::Float64
    
    # Bioaccumulation Properties
    "Reference depuration (elimination) rate at 20°C. Unit: \$day^{-1}\$."
    k_dep_20::Float64
    "Temperature adjustment coefficient for depuration. Unit: dimensionless."
    theta_dep::Float64
    "TSS concentration at which pseudofeces rejection begins. Unit: mg/L."
    tss_rejection_threshold::Float64
    "TSS concentration at which rejection is maximal (gills clog). Unit: mg/L."
    tss_clogging_threshold::Float64
end


"""
    f_temp(temperature_c)

Calculate the temperature limitation factor `f(T)`.
This is a Gaussian-like function centered at 27°C.
"""
function f_temp(temperature_c::Real)
    return exp(-0.006 * (temperature_c - 27.0)^2)
end

"""
    f_salinity(salinity_psu)

Calculate the salinity limitation factor `f(S)`.
This is a piecewise function representing the oyster's tolerance to salinity.
"""
function f_salinity(salinity_psu::Real)
    if salinity_psu < 5.0
        return 0.0
    elseif 5.0 <= salinity_psu <= 12.0
        return 0.0926 * (salinity_psu - 0.0139)
    else # salinity_psu > 12.0
        return 1.0
    end
end

"""
    f_tss(tss_mg_per_l)

Calculate the Total Suspended Solids limitation factor `f(TSS)`.
This piecewise function models the effect of particle concentration on filtration.
"""
function f_tss(tss_mg_per_l::Real)
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
    oyster_params::OysterParameters,
    temperature_c::Real,
    salinity_psu::Real,
    tss_mg_per_l::Real
    )
    allometric_scaling = oyster_params.dry_weight_g^0.75

    fr = 0.17 *
         allometric_scaling *
         f_temp(temperature_c) *
         f_salinity(salinity_psu) *
         f_tss(tss_mg_per_l)

    return fr
end

end