module PreVirS
    include("oyster.jl")
    include("oysterData.jl")
    include("VirusFateModel.jl")
    include("OysterBioaccumulation.jl")

    using .oyster
    export OysterParameters
    export f_temp
    export f_salinity
    export f_tss
    export calculate_filtration_rate

    using .oysterData
    export simulate_yearly_environment
    export generate_synthetic_data
    export generate_virus_influx

    using .VirusFateModel
    export VirusParameters # struct
    export EnvironmentalConditions # struct
    export WaterSimulationResult # struct
    export run_water_simulation # main sim virus decay
    export simulate_water_dynamics #24  hour sim

    using .OysterBioaccumulation
    export simulate_oyster_concentration # main sim oyster accumulation

end