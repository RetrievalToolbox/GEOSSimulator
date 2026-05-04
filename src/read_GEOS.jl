"""
Calculates the coordinate tree for a GEOS file, needed for spatial interpolation at the
requrested sample locations.
"""
function build_tree(nc)

    # Find out if this is a cube-sphere file or not
    cube_sphere = "nf" in keys(nc.dim)

    if cube_sphere

        # Create a ball tree for every face
        face_trees = BallTree[]
        for idx_face in 1:6
            lonlat = hcat(vec(nc["lons"][:,:,idx_face]), vec(nc["lats"][:,:,idx_face]))'
            tree = BallTree(lonlat, Haversine())
            push!(face_trees, tree)
        end

        return (true, face_trees)

    else

        lon = nc["lon"][:]
        Nlon = length(lon)
        lat = nc["lat"][:]
        Nlat = length(lat)

        lonlat = zeros(2, Nlon * Nlat)
        idx_2d = CartesianIndices((Nlon, Nlat))[1:Nlon*Nlat]

        lonlat[1, :] = lon[(x -> x.I[1]).(idx_2d)]
        lonlat[2, :] = lat[(x -> x.I[2]).(idx_2d)]

        tree = BallTree(lonlat, Haversine())

        return (false, tree)

    end

end

"""
Pre-calculates the NN spatial indices needed to sample the model locations, along with the
distances between the sample locations and the NN closest model grid points.
"""
function calculate_spatial_indices(
    longitude,
    latitude,
    tree::Vector{BallTree},
    GEOS_dims::Tuple,
    NN::Integer
    )

    # This function is for calculating the spatial indices for a cube-sphere situation
    N_face = length(tree)

    indices = zeros(Int, NN, N_face)
    distances = zeros(NN, N_face)

    for i_face in 1:N_face

        # Calculate nearest neighbors
        nn = NearestNeighbors.knn(tree[i_face], [longitude, latitude], NN)

        indices[:, i_face] .= nn[1]
        distances[:, i_face] .= nn[2]

    end

    # Sort, find NN smallest distances, obtain linear index
    linear_indices = partialsortperm(vec(distances), 1:NN)

    # What's the index within the "local" distance matrix (neighbor, face)
    small_idx = CartesianIndices((NN, N_face))[linear_indices]
    face_idx = map(x -> x.I[2], small_idx)

    # What's the index in the GEOS model file PER FACE (xdim, ydim)
    GEOS_idx_per_face = CartesianIndices(GEOS_dims)[indices[small_idx]]

    # And finally construct full index for (xdim, ydim, face)
    GEOS_idx = [
        CartesianIndex(GEOS_idx_per_face[i].I..., face_idx[i])
        for i in 1:NN
    ]

    return GEOS_idx, distances[small_idx]

end

function calculate_spatial_indices(
    longitude,
    latitude,
    tree::BallTree,
    GEOS_dims::Tuple,
    NN::Integer
    )

    indices = zeros(Int, NN)
    distances = zeros(NN)

    # Calculate nearest neighbors
    indices, distances = NearestNeighbors.knn(tree, [longitude, latitude], NN)

    # What's the index in the GEOS model file for the regular lon/lat grid
    GEOS_idx = CartesianIndices(GEOS_dims)[indices]

    return GEOS_idx, distances

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


"""
Helper function to sum over the first dimension and return a collapsed vector
"""
wsum(x::Matrix) = dropdims(sum(x, dims=1), dims=1)

"""
AK/BK coefficients from BradW. These are not (yet) included in the GEOS-Carb
product files.
"""
AK = [1, 2.00000023841858, 3.27000045776367, 4.75850105285645,
    6.60000133514404, 8.93450164794922, 11.9703016281128, 15.9495029449463,
    21.1349029541016, 27.8526058197021, 36.5041084289551, 47.5806083679199,
    61.6779098510742, 79.5134124755859, 101.944023132324, 130.051025390625,
    165.079025268555, 208.497039794922, 262.021057128906, 327.64306640625,
    407.657104492188, 504.680114746094, 621.680114746094, 761.984191894531,
    929.294189453125, 1127.69018554688, 1364.34020996094, 1645.71032714844,
    1979.16040039062, 2373.04052734375, 2836.78051757812, 3381.00073242188,
    4017.541015625, 4764.39111328125, 5638.791015625, 6660.34130859375,
    7851.2314453125, 9236.572265625, 10866.3017578125, 12783.703125,
    15039.302734375, 17693.00390625, 20119.201171875, 21686.501953125,
    22436.30078125, 22389.80078125, 21877.59765625, 21214.998046875,
    20325.8984375, 19309.6953125, 18161.896484375, 16960.896484375,
    15625.99609375, 14290.9951171875, 12869.59375, 11895.8623046875,
    10918.1708984375, 9936.521484375, 8909.9921875, 7883.421875,
    7062.1982421875, 6436.263671875, 5805.3212890625, 5169.61083984375,
    4533.90087890625, 3898.20092773438, 3257.08081054688, 2609.20068359375,
    1961.310546875, 1313.48034667969, 659.375244140625, 4.80482578277588, 0]u"Pa"

BK = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8.17541323527848e-09,
    0.00696002459153533, 0.0280100405216217, 0.0637200623750687,
    0.113602079451084, 0.156224086880684, 0.200350105762482,
    0.246741116046906, 0.294403105974197, 0.343381136655807,
    0.392891138792038, 0.44374018907547, 0.494590193033218,
    0.546304166316986, 0.581041514873505, 0.615818440914154,
    0.650634944438934, 0.685899913311005, 0.721165955066681,
    0.749378204345703, 0.770637512207031, 0.791946947574615,
    0.81330394744873, 0.834660947322845, 0.856018006801605,
    0.877429008483887, 0.898908019065857, 0.920387029647827,
    0.941865026950836, 0.963406026363373, 0.984951972961426, 1
]

function generate_scenes_from_GEOS(
    global_config::RESimulatorCore.SimulatorGlobalConfig,
    GEOSIT_fdict::Dict,
    GEOSCARB_fdict::Dict,
    buffer::EarthAtmosphereBuffer,
    lonlat_array,
    dtime::DateTime;
    NN::Integer=1
    )

    #=
        GEOS IT section
        ===============
    =#

    nc_asm_lvl = NCDataset(GEOSIT_fdict["asm_lvl"])

    # Create the coordinate trees for interpolation
    GEOSIT_cube_sphere, GEOSIT_tree = build_tree(nc_asm_lvl)

    # How many spatial dimensions do we have?
    GEOSIT_ndim_spatial = GEOSIT_cube_sphere ? 3 : 2

    GEOSIT_slicer_lvl = [Colon() for _ in 1:GEOSIT_ndim_spatial]

    # Need the GEOS spatial dimensions to construct proper indexing
    GEOSIT_dims = size(nc_asm_lvl["lons"])[1:2]

    # Load the needed variable arrays into memory..
    # (consumes more memory than reading element by element, but faster)
    GEOSIT_PS = nc_asm_lvl["PS"][GEOSIT_slicer_lvl...,1] * get_unit(nc_asm_lvl["PS"])
    GEOSIT_PL = nc_asm_lvl["PL"][GEOSIT_slicer_lvl...,:,1] * get_unit(nc_asm_lvl["PL"])
    GEOSIT_DELP = nc_asm_lvl["DELP"][GEOSIT_slicer_lvl...,:,1] * get_unit(nc_asm_lvl["DELP"])
    GEOSIT_T = nc_asm_lvl["T"][GEOSIT_slicer_lvl...,:,1] * get_unit(nc_asm_lvl["T"])
    GEOSIT_QV = nc_asm_lvl["QV"][GEOSIT_slicer_lvl...,:,1] * get_unit(nc_asm_lvl["QV"])

    #=
        GEOS CARB section
        ===============
    =#

    nc_carb_lvl = NCDataset(GEOSCARB_fdict["carb_lvl"])
    # Create the coordinate trees for interpolation
    GEOSCARB_cube_sphere, GEOSCARB_tree = build_tree(nc_carb_lvl)
    GEOSCARB_ndim_spatial = GEOSCARB_cube_sphere ? 3 : 2
    GEOSCARB_slicer_lvl = [Colon() for _ in 1:GEOSCARB_ndim_spatial]
    GEOSCARB_dims = (nc_carb_lvl.dim["lon"], nc_carb_lvl.dim["lat"])

    GEOSCARB_PS = nc_carb_lvl["PS"][GEOSCARB_slicer_lvl...,1] * get_unit(nc_carb_lvl["PS"])
    GEOSCARB_CO2 = nc_carb_lvl["CO2"][GEOSCARB_slicer_lvl...,:,1] * get_unit(nc_carb_lvl["CO2"])

    # Obtain all gas objects into a gas-name, gas-object lookup
    gas_list = filter(x -> x isa RE.GasAbsorber, buffer.scene.atmosphere.atm_elements)
    gas_dict = Dict(g.gas_name => g for g in gas_list)

    # Interpolation weights - inverse distance squared!
    # (pre-allocate)
    weights = zeros(NN)

    # Create an empty vector of scene configurations
    GEOS_scenes = RESimulatorCore.SimulatorSceneConfig[]

    @info "Generating scene information from GEOS output"
    @showprogress for idx_scene in axes(lonlat_array, 2)

        scene_lon = lonlat_array[1, idx_scene]
        scene_lat = lonlat_array[2, idx_scene]

        # Calculate the spatial indices and distances between model grid and
        # the sample locations

        GEOSIT_idx_space, GEOSIT_distances = calculate_spatial_indices(
            scene_lon, scene_lat,
            GEOSIT_tree, GEOSIT_dims, NN
        )


        @views @. weights[:] = 1 / (GEOSIT_distances ^ 2)
        weights[:] ./= sum(weights)

        #=
            Sample GEOS IT data, and use weights to spatially interpolate at
            the actual sampling location.
        =#

        # Mid-level pressure for MET profiles
        met_pressure = GEOSIT_PL[GEOSIT_idx_space,:] .* weights |> wsum

        # Temperature
        temperature = GEOSIT_T[GEOSIT_idx_space,:] .* weights |> wsum

        # QV (specific humidity)
        specific_humidity = GEOSIT_QV[GEOSIT_idx_space,:] .* weights |> wsum


        # Construct RT pressure grid via surface pressure and Δp
        # Surface pressure
        PS = GEOSIT_PS[GEOSIT_idx_space] .* weights |> sum
        # Δp
        DELP = GEOSIT_DELP[GEOSIT_idx_space,:] .* weights |> wsum

        pressure_levels = zeros(global_config.atmosphere.N_RT_level) *
            global_config.atmosphere.pressure_unit

        pressure_levels[end] = PS

        for i in global_config.atmosphere.N_RT_level-1:-1:1
            pressure_levels[i] = max(pressure_levels[i+1] - DELP[i], 0.001u"Pa")
        end

        GEOSCARB_idx_space, GEOSCARB_distances = calculate_spatial_indices(
            scene_lon, scene_lat,
            GEOSCARB_tree, GEOSCARB_dims, NN
        )

        @views @. weights[:] = 1 / (GEOSCARB_distances ^ 2)
        weights[:] ./= sum(weights)

        # GEOS CARB model levels, created from PS in Pa and AK, BK
        CO2_PS = GEOSCARB_PS[GEOSCARB_idx_space,:] .* weights |> wsum
        CO2_PL_EDGES = CO2_PS .* BK .+ AK # this is unit-aware!
        CO2_PL = @. 0.5 * (CO2_PL_EDGES[1:end-1] + CO2_PL_EDGES[2:end])

        # CO2 is on MODEL LEVELS
        CO2 = GEOSCARB_CO2[GEOSCARB_idx_space,:] .* weights |> wsum

        # ... which we must inter/extrapolate onto the RETRIEVAL GRID
        CO2_itp = linear_interpolation(
            CO2_PL, CO2,
            extrapolation_bc=Flat()
        )
        CO2_on_rgrid = CO2_itp.(pressure_levels)

        # Add fake surface for now..
        surface_parameters = [(0.25,) for swin in global_config.spectral_windows]

        # Gas profile dictionary
        vmr_levels = Dict{String, Vector{Float64}}()

        # We set O2 always as 0.2095 parts
        if "O2" in keys(gas_dict)
            gas = gas_dict["O2"]
            vmr_levels["O2"] = ustrip.(
                Ref(gas.vmr_unit),
                fill(0.2095, global_config.atmosphere.N_RT_level)
            )
        end

        if "CO2" in keys(gas_dict)
            gas = gas_dict["CO2"]
            vmr_levels["CO2"] = ustrip.(
                Ref(gas.vmr_unit),
                CO2_on_rgrid
            )
        end

        # This is where we create the scene configuration ..

        this_scene = RESimulatorCore.SimulatorSceneConfig(
            date=dtime,
            ######################
            solar_zenith_angle=0.0,
            solar_azimuth_angle=0.0,
            #######################
            loc_longitude=scene_lon,
            loc_latitude=scene_lat,
            loc_elevation=0.0u"m",
            #####################################
            surface_parameters=surface_parameters,
            ###############################
            pressure_levels=pressure_levels,
            #####################
            vmr_levels=vmr_levels,
            #######################################
            met_pressure=met_pressure,
            specific_humidity=specific_humidity,
            temperature=temperature
        )

        # .. and move it into the list
        push!(GEOS_scenes, this_scene)

    end

    #=
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
    =#

    return GEOS_scenes

end
