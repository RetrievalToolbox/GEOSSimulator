"""
Turns a string like `2020-03-01_18:00:00` into a `DateTime` object.
"""
function parse_time_string(s::String)

    year = parse(Int, s[1:4])
    # 5: dash
    month = parse(Int, s[6:7])
    # 8: dash
    day = parse(Int, s[9:10])
    #11: underscore
    hour = parse(Int, s[12:13])
    #14: colon
    minute = parse(Int, s[15:16])
    #17: colon
    second = parse(Int, s[18:19])

    return DateTime(year, month, day, hour, minute, second)

end

function build_tree(
    WRF_longitudes,
    WRF_latitudes
)

    # Build 1D array of lon/lat coordinates
    lonlat = hcat(vec(WRF_longitudes[:,:,1]), vec(WRF_latitudes[:,:,1]))'

    # Produce tree and return
    return BallTree(lonlat, Haversine())

end

"""
Calculates temperature from WRF-produced potential temperature perturbation
"""
function T_from_Δθ(
    ΔΘ::Unitful.Temperature,
    p0::Unitful.Pressure,
    p::Unitful.Pressure
    )

    R_air = 287.052874247u"J/kg/K"
    cp_air = 1004.0u"J/kg/K"

    Θ = ΔΘ + 300.0u"K"

    return Θ / ((p0 / p)^(R_air / cp_air))
end

function grab_scalar_from_2d(arr, idx_vector, w)

    @assert length(size(arr)) == 2

    result = zero(eltype(arr))
    for (i, idx) in enumerate(idx_vector)
        result += arr[idx] * w[i]
    end
    return result

end

function grab_vector_from_3d(arr, idx_vector, w)

    @assert length(size(arr)) == 3

    result = zeros(eltype(arr), size(arr, 3))
    for (i, idx) in enumerate(idx_vector)
        @views result[:] .+= arr[idx, :] .* w[i]
    end
    return result

end


function get_unit(nc_var)

    # Grab the unit attribute as a string
    unit_str = nc_var.attrib["units"]

    # Regex replace notation to work with Unitful:
    # e.g.:
    # kg kg-1 => kg kg^-1
    # kg m s-2 => kg m s^-2
    # etc..

    # Also replace spaces with multiplication, so that
    # kg kg^-1 => kg * kg^-1

    unit_str = replace(unit_str,
        r"(-?\d)" => s"^\1",
        " " => " * "
    )

    # Other known replacements
    unit_str = replace(unit_str,
        "ppmv" => "ppm",
        "ppbv" => "ppb",
        "pptv" => "ppt"
    )


    return uparse(unit_str)

end


function generate_scenes_from_WRF(
    global_config::RESimulatorCore.SimulatorGlobalConfig,
    WRF_fname::String,
    buffer::EarthAtmosphereBuffer,
    lonlat_array;
    NN::Integer=1
    )


    # Create an empty vector of scene configurations
    WRF_scenes = RESimulatorCore.SimulatorSceneConfig[]


    nc = NCDataset(WRF_fname)

    # How many time steps in this file?
    N_time = size(nc["Times"], 2)
    N_x = nc.dim["west_east"]
    N_y = nc.dim["south_north"]

    WRF_dims = (N_x, N_y)
    WRF_times = nc["Times"][:,:] # char array, WRF_times
    WRF_longitudes = nc["XLONG"][:,:,1] # x, y, time
    WRF_latitudes = nc["XLAT"][:,:,1] # x, y, time

    # Generate coordinate tree to use for spatial interpolation
    coord_tree = build_tree(WRF_longitudes, WRF_latitudes)

    # x, y, time
    WRF_surf_altitudes = nc["HGT"][:,:,1] * get_unit(nc["HGT"])
    # x, y, time
    WRF_surf_pressure = nc["PSFC"][:,:,1] * get_unit(nc["PSFC"])

    # WRF pressure levels are a sum of base pressure (PB) and perturbation pressure (P)
    WRF_pressure = (
        nc["P"][:,:,:,1] .+ # x, y, level, time
        nc["PB"][:,:,:,1] # x, y, level, time
    ) * get_unit(nc["P"])
    # must reverse order
    reverse!(WRF_pressure, dims=3)


    # WRF water vapor mixing ratio (mass of water / mass of dry air)
    # x, y, level, time
    WRF_H2O_VMR = nc["QVAPOR"][:,:,:,1] * get_unit(nc["QVAPOR"])
    # ===> convert to specific humidity Q (mass of water / mass of moist air)
    WRF_Q = @. WRF_H2O_VMR / (1 + WRF_H2O_VMR)
    # must reverse order
    reverse!(WRF_Q, dims=3)

    # Let's add up all CO2 contributions
    # (TODO: no idea how these work .. adding them up is clearly wrong)
    WRF_CO2 = (
        nc["CO2_BCK"][:,:,:,1] # background
        #nc["CO2_BIO"][:,:,:,1] .+ # biogenic
        #nc["CO2_OCE"][:,:,:,1] .+ # ocean
        #nc["CO2_ANT"][:,:,:,1] .+ # anthropogenic
        #nc["CO2_BBU"][:,:,:,1] # biomass burning
    ) * get_unit(nc["CO2_BCK"])
    # must reverse order
    reverse!(WRF_CO2, dims=3)



    # WRF potential temperature perturbation
    # (this must be constructed)
    WRF_pot_ΔT = nc["T"][:,:,:,1] * get_unit(nc["T"])
    # must reverse order
    reverse!(WRF_pot_ΔT, dims=3)

    # Surface albedo
    WRF_albedo = nc["ALBEDO"][:,:,1]

    # Obtain all gas objects into a gas-name, gas-object lookup
    gas_list = filter(x -> x isa RE.GasAbsorber, buffer.scene.atmosphere.atm_elements)
    gas_dict = Dict(g.gas_name => g for g in gas_list)

    @info "Generating scenes .."

    # Pre-allocate some vectors, we just re-use them inside the loop
    weights = zeros(NN)

    for idx_scene in axes(lonlat_array, 2)

        # ==============================
        # Reading fields from WRF arrays
        # ==============================
        scene_datetime = WRF_times[:,1] |> String |> parse_time_string
        scene_lon = lonlat_array[1, idx_scene]
        scene_lat = lonlat_array[2, idx_scene]

        # Calculate array indices from lon/lat
        idx_flat, distances = knn(coord_tree, [scene_lon, scene_lat], NN)
        @views @. weights[:] = 1 / (distances ^ 2)
        @views weights[:] ./= sum(weights)

        # Flatten index
        idx_WRF = CartesianIndices(WRF_dims)[idx_flat]

        # Sample from array
        scene_altitude = grab_scalar_from_2d(WRF_surf_altitudes, idx_WRF, weights)

        # =======================
        # Generate for this scene
        # =======================

        # Grab the meteorological pressure level
        met_pressure_levels = grab_vector_from_3d(WRF_pressure, idx_WRF, weights)

        psurf = grab_scalar_from_2d(WRF_surf_pressure, idx_WRF, weights)
        #pressure_levels = RE.create_ACOS_pressure_grid(psurf)
        pressure_levels = copy(met_pressure_levels)

        # Grab the specific humidity
        specific_humidity_levels = grab_vector_from_3d(WRF_Q, idx_WRF, weights)

        # .. and Δθ
        pot_ΔT_levels = grab_vector_from_3d(WRF_pot_ΔT, idx_WRF, weights)
        # which we now turn into T
        temperature_levels = T_from_Δθ.(pot_ΔT_levels, Ref(psurf), met_pressure_levels)

        # Create surface parameters
        # (just use same albedo for all bands for now..)
        #surface_parameters = create_surface_parameters(global_config.spectral_windows)
        albedo = grab_scalar_from_2d(WRF_albedo, idx_WRF, weights)
        surface_parameters = [(albedo,) for swin in global_config.spectral_windows]

        # Gas profile dictionary
        vmr_levels = Dict{String, Vector{Float64}}()

        # We set O2 always as 0.2095 parts
        if "O2" in keys(gas_dict)
            vmr_levels["O2"] = fill(0.2095, global_config.atmosphere.N_RT_level)
        end

        if "CO2" in keys(gas_dict)
            vmr_levels["CO2"] = grab_vector_from_3d(WRF_CO2, idx_WRF, weights)
        end

        # This is where we create the scene configuration ..
        this_scene = RESimulatorCore.SimulatorSceneConfig(
            date=scene_datetime,
            ######################
            solar_zenith_angle=0.0,
            solar_azimuth_angle=0.0,
            #######################
            loc_longitude=scene_lon,
            loc_latitude=scene_lat,
            loc_altitude=scene_altitude,
            #####################################
            surface_parameters=surface_parameters,
            ###############################
            pressure_levels=pressure_levels,
            #####################
            vmr_levels=vmr_levels,
            #######################################
            met_pressure_levels=met_pressure_levels,
            specific_humidity_levels=specific_humidity_levels,
            temperature_levels=temperature_levels
        )

        # .. and move it into the list
        push!(WRF_scenes, this_scene)

    end

    close(nc)

    return WRF_scenes

end