using Pkg
#Pkg.develop(path=joinpath(@__DIR__, ".."))
Pkg.activate(joinpath(@__DIR__, ".."))

using Revise
using PreVirS

dataset, header = generate_oyster_dataset(1000)

simulation_year = 2024
fixed_oyster_weight_g = 0.5 

# Generate the environmental data
timestamps, temps, salinities, tsss = simulate_yearly_environment(simulation_year)
# Calculate the filtration rate for every hour
# The '.' broadcasts the function over the input vectors

filtration_rates = calculate_filtration_rate.(fixed_oyster_weight_g, temps, salinities, tsss)

plot_simulation_timeseries(timestamps, temps, salinities, tsss, filtration_rates)