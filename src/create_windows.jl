function create_windows(fname::String)

    yml = YAML.load_file(fname)

    swins = RE.SpectralWindow[]

    for wname in keys(yml)

        # Try wavelength
        if haskey(yml[wname], "wavelength_unit")
            ww_key = "wavelength"
        elseif haskey(yml[wname], "wavenumber_unit")
            ww_key ="wavenumber"
        else
            error("Window specification must have either `wavenumber_unit` or \
                `wavelength_unit`!")
        end

        ww_unit = uparse(yml[wname]["$(ww_key)_unit"])
        ww_min = yml[wname]["$(ww_key)_min"]
        ww_max = yml[wname]["$(ww_key)_max"]
        ww_ref = yml[wname]["$(ww_key)_ref"]
        ww_buf = yml[wname]["$(ww_key)_buffer"]
        ww_spacing = yml[wname]["$(ww_key)_spacing"]

        # Create the grid here
        ww_grid = collect(ww_min-ww_buf:ww_spacing:ww_max+ww_buf+ww_spacing)

        this_swin = RE.SpectralWindow(
            wname,
            ww_min,
            ww_max,
            ww_grid,
            ww_unit,
            ww_ref
        )

        push!(swins, this_swin)

    end

    return swins

end