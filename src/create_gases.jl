function create_gases(
        fname::String,
        spectral_unit::Symbol,
        levels::Integer;
        distributed::Bool=true
    )

    yml = YAML.load_file(fname)
    gases = RE.GasAbsorber[]

    for gas_name in keys(yml)

        @assert haskey(yml[gas_name], "type") "Gas $(gas_name) does not have required key `type`!"
        @assert haskey(yml[gas_name], "spectroscopy") "Gas $(gas_name) does not have required key `spectroscopy`!"
        @assert haskey(yml[gas_name], "unit") "Gas $(gas_name) does not have required key `unit`!"

        # What type of spectroscopy do we use for this gas?
        spec_type = yml[gas_name]["type"]

        # Load spectroscopy table
        if spec_type == "ABSCO"

            spec = RE.load_ABSCO_spectroscopy(
                yml[gas_name]["spectroscopy"];
                spectral_unit=spectral_unit,
                distributed=distributed
            )

        else
            @error "Sorry, spectroscopy type $(type) is not yet supported!"
            exit(1)
        end


        vmr_levels = zeros(levels)
        # Try and parse the supplied unit
        if isnothing(yml[gas_name]["unit"])
            vmr_unit = Unitful.NoUnits
        else
            vmr_unit = uparse(yml[gas_name]["unit"])
        end

        gas = RE.GasAbsorber(
            gas_name,
            spec,
            vmr_levels,
            vmr_unit
        )

        push!(gases, gas)
        @info " .. created gas: $(gas)"
    end

    return gases

end