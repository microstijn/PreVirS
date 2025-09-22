module PreVirS
    include("oyster.jl")
    include("oysterData.jl")

    using .oyster
    export f_temp
    export f_salinity
    export f_tss
    export calculate_filtration_rate

    using .oysterData
    export simulate_yearly_environment
    export plot_simulation_timeseries

end