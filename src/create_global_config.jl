function create_global_config(
    window_yml_fname::String,
    gas_yml_fname::String,
    solar_model_fname::String,
    N_RT_level::Integer,
    N_met_level::Integer
)

    spectral_windows = create_windows(window_yml_fname)

    # Short check .. for now let's demand that ALL spectral windows have the same
    # spectral axis (either wavenumbers or wavelengths)
    have_wavelengths = (x -> isa(x.ww_unit, Unitful.LengthUnits)).(spectral_windows) |> all
    have_wavenumbers = (x -> isa(x.ww_unit, Unitful.WavenumberUnits)).(spectral_windows) |> all

    if (have_wavelengths)
        @info "All spectral windows are wavelenghts: ✓"
        spectral_unit = :Wavelength
    elseif all(have_wavenumbers)
        @info "All spectral windows are wavenumbers: ✓"
        spectral_unit = :Wavenumber
    else

        @error "Your supplied spectral windows have a mix between \
            wavelengths and wavenumbers: " * (
                ["\n$(swin) → $(swin.ww_unit)" for swin in spectral_windows] |> join
            )

        exit(1)
    end

    @info "Creating surfaces"
    surfaces = [
        (:Lambert, 1) for _ in 1:length(spectral_windows)
    ]

    @info "Creating spectroscopy and gases"
    gases = create_gases(
        gas_yml_fname,
        spectral_unit,
        N_RT_level,
    )

    @info "Creating atmosphere configuration"
    atmosphere_config = RESimulatorCore.SimulatorAtmosphereConfig(
        N_met_level=N_met_level,
        N_RT_level=N_RT_level,
        elements=[RE.RayleighScattering(), gases...]
        # leave units with default settings for now
    )

    @info "Creating solar configuration"
    solar_config = RESimulatorCore.SimulatorSolarConfig(
        solar_model_type="TSIS",
        solar_model_path=solar_model_fname,
    )

    @info "Creating RT configuration"
    #=
    RT_config = RESimulatorCore.SimulatorRTConfig(
        models=[:XRTM for _ in 1:length(spectral_windows)],
        model_options=[
            [
                Dict(
                    "solvers" => ["single"],
                    "add" => true, # zero-out RT container
                    "sun_normalized" => false, # use the actual solar spectrum
                    "streams" => 8,
                    "options" => [
                        "output_at_levels",
                        "source_solar",
                        "vector",
                        "psa", # pseudo-spherical approximation for incoming ray
                        "sfi",
                    ]
                ),
                #=
                Dict(
                    "solvers" => ["two_stream"],
                    "add" => true, # zero-out RT container
                    "sun_normalized" => false, # use the actual solar spectrum
                    "streams" => 2,
                    "options" => [
                        "output_at_levels",
                        "source_solar",
                        #"vector",
                        "psa", # pseudo-spherical approximation for incoming ray
                        "sfi",
                    ]
                ),
                =#
            ] for _ in 1:length(spectral_windows)
        ]
    )
    =#

    RT_config = RESimulatorCore.SimulatorRTConfig(
        models=[:BeerLambert for _ in 1:length(spectral_windows)],
        model_options=[[Dict()] for _ in 1:length(spectral_windows)]
    )


    # Create a global config, this is from RESimulatorCore
    @info "Creating global configuration"
    global_config = RESimulatorCore.SimulatorGlobalConfig(
        spectral_unit=spectral_unit,
        spectral_windows=spectral_windows,
        surfaces=surfaces,
        atmosphere=atmosphere_config,
        RT=RT_config,
        solar=solar_config
    )

    return global_config

end
