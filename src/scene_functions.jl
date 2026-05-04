function solar_zenith_from_scene(scene::RESimulatorCore.SimulatorSceneConfig)
    sza, saa = RE.calculate_solar_angles(scene.loc_longitude, scene.loc_latitude, scene.date)
    return sza
end

function valid_solar_zenith_from_scene(scene::RESimulatorCore.SimulatorSceneConfig)
    sza = solar_zenith_from_scene(scene)

    return (sza >= 0) & (sza < 90)
end

"""
Filters a vector of scenes by whether their location/time combo produces a valid solar
zenith angle. Scenes with valid solar zenith angles are returned.
"""
function filter_by_solar_zenith(scenes::Vector{RESimulatorCore.SimulatorSceneConfig})

    return filter(x -> valid_solar_zenith_from_scene(x), scenes)
   
end
