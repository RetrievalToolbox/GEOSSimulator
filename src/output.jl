function unit_str(un::Unitful.Units)

    # First, grab the units as a non-fancy string (no UTF exponents)
    u1 = sprint(show, un, context=:fancy_exponent => false)
    # Now remove all ^
    u2 = replace(u1, "^" => "")

    return u2

end


function write_out(
    fname::String,
    global_config::RESimulatorCore.SimulatorGlobalConfig,
    scene_configs::Vector{RESimulatorCore.SimulatorSceneConfig},
    results::Vector,
    buffer::RE.EarthAtmosphereBuffer
    )

    @info "Writing out to NetCDF file at: $(fname)"

    # Create the scene dimension
    scene_dim = length(results)
    pol_dim = 3 # hard-code for now?

    NCDataset(fname, "c") do ds


        # Drop in the per-scene data. This is handled by RESimulatorCore
        RESimulatorCore.write_scenes_into_nc(ds, global_config, scene_configs)

        # The other results (radiances mostly) we have to do ourselves

        # File-wide dimensions
        defDim(ds, "scene", scene_dim)
        defDim(ds, "polarization", pol_dim)


        for (i_swin, swin) in enumerate(global_config.spectral_windows)
            @info "Writing out $(swin)"
            spec_dim = swin.N_hires

            # Create a group for every spectral window
            grp = defGroup(ds, swin.window_name)

            # Create the per-window spectral coordinate
            defDim(grp, "spectral", spec_dim)

            # Write out the spectral grid
            defVar(grp, "spectral", swin.ww_grid, ("spectral", ),
                attrib = OrderedDict(
                    "units" => unit_str(swin.ww_unit),
                    "polarization_order" => "I, Q, U",
                ))

            # Create the radiance output variable
            rad = defVar(grp, "radiance", Float64, ("spectral", "polarization", "scene"),
                attrib = OrderedDict(
                    "units" => unit_str(buffer.rt[swin].radiance_unit),
                )
            )

            for i_scene in 1:scene_dim
                rad[:,:,i_scene] .= results[i_scene][i_swin]
            end


        end

    end
end