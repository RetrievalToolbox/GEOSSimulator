module GEOSSimulator

    using ArgParse
    using Dates
    using Distances
    using Distributed
    #using FastInterpolations
    using Interpolations
    using NCDatasets
    using NearestNeighbors
    using OrderedCollections
    using Printf
    using ProgressMeter
    using Unitful
    using YAML

    using RetrievalToolbox
    const RE = RetrievalToolbox

    using RESimulatorCore


    include("file_functions.jl")
    include("create_windows.jl")
    include("create_gases.jl")
    include("surface_parameters.jl")

    include("create_global_config.jl")
    include("read_GEOS.jl")
    include("output.jl")


    function parse_commandline()

        s = ArgParseSettings(description = "WRF RT Simulator")

        @add_arg_table! s begin
            "--GEOSIT"
                help = "Path to the top-level GEOS IT directory (e.g. /css/gmao/geos-it/products)"
                arg_type = String
                required = true
            "--GEOSCARB"
            help = "Path to the top-level GEOS CARB directory (e.g. /css/gmao/geos_carb/pub/carbon/vNRT/netcdf/ghg_inst_3hr_glo_L2880x1441_v72)"
                arg_type = String
                required = true
            "-Y"
                help = "Year"
                arg_type = Int
                required = true
            "-M"
                help = "Month"
                arg_type = Int
                required = true
            "-D"
                help = "Day"
                arg_type = Int
                required = true
            "--hour"
                help = "Hour"
                arg_type = Int
                required = true
            "--minute"
                help = "Minute"
                arg_type = Int
                required = true
            "--windows", "-w"
                help = "Path to windows YML file"
                required = true
            "--gases", "-g"
                help = "Path to gases YML file"
                required = true
            "--TSIS"
                help = "Path to TSIS solar model file"
                required = true
            "--coords"
                help = "Path to CSV file containing the lon/lat pairs"
                required = true
            "--output", "-o"
                help = "Path to save the simulated radiances"
                required = true
            "--uplooking"
                help = "Uplooking mode for (e.g. TCCON observations)"
                arg_type = Bool
                default = false
            "--procs", "-p"
                help = "Number of *additional* processes to use (default 0)"
                arg_type = Int
                default = 0
            "--neighbors", "-N"
                help = "Number of nearest neighbors for spatial interpolation (default 1)"
                arg_type = Int
                default = 1
        end

        return parse_args(ARGS, s)
    end

    function main()

        # Check for XRTM
        if !haskey(ENV, "XRTM_PATH")
            error("XRTM_PATH environmental variable not set! \
                Please set it to the location of the XRTM library.")
            exit(1)
        end


        # Parse arguments
        args = parse_commandline()

        # NN must be >= 0
        if args["neighbors"] < 1
            error("Number of nearest neighbors for spatial interpolation must be >= 1.")
        end

        if args["procs"] > 0
            @info "Adding $(args["procs"]) extra processes."
            addprocs(args["procs"])
        end

        @everywhere @eval using RetrievalToolbox
        @everywhere @eval const RE = RetrievalToolbox
        @everywhere @eval using RESimulatorCore
        @everywhere @eval using GEOSSimulator

        #=
            FIND GEOS IT FILES
            ==================
        =#

        # Geos filenames for released produces are:
        # (e.g.)
        # GEOS.it.asm.asm_inst_1hr_glo_L576x361_slv.GEOS5294.2025-04-03T1900.V01.nc4
        # so the timestamp needs to be in form YYYY-MM-DDThhmm

        GEOSIT_tstamp = @sprintf(
            "%04d-%02d-%02dT%02d%02d",
            args["Y"],
            args["M"],
            args["D"],
            args["hour"],
            args["minute"],
        )

        GEOSIT_dir = joinpath(
            args["GEOSIT"],
            @sprintf("Y%04d", args["Y"]),
            @sprintf("M%02d", args["M"]),
            @sprintf("D%02d", args["D"])
        )

        GEOSIT_pattern = Dict(
            # needed for T, PS, QV profiles
            "asm_lvl" => Regex("asm\\.asm_inst_.+_C\\d+x\\d+x6_v\\d+.+$(GEOSIT_tstamp).+"),
            # add aerosols etc. later..
        )

        # Find the GEOSIT files needed to be read in for this date:
        GEOSIT_fdict = create_GEOS_filedict(
            GEOSIT_dir,
            GEOSIT_tstamp,
            GEOSIT_pattern
        )


        #=
            FIND GEOS CARB FILES
            ====================
        =#

        # GEOS CARB timestamps are like 20260331_2100z

        GEOSCARB_tstamp = @sprintf(
            "%04d%02d%02d_%02d%02dz",
            args["Y"],
            args["M"],
            args["D"],
            args["hour"],
            args["minute"],
        )

        GEOSCARB_dir = joinpath(
            args["GEOSCARB"],
            @sprintf("%04d", args["Y"]),
            @sprintf("%02d", args["M"])
        )

        GEOSCARB_pattern = Dict(
            # needed for CO2, CO
            "carb_lvl" => Regex(".+ghg_inst_3hr_glo_L\\d+x\\d+_v\\d+.+$(GEOSCARB_tstamp).+"),
        )

        GEOSCARB_fdict = create_GEOS_filedict(
            GEOSCARB_dir,
            GEOSCARB_tstamp,
            GEOSCARB_pattern
        )

        # Create the global config on the main worker, and then copy to all others. This
        # ensures that the spectroscopy objects will be shared across all workers, while
        # only one copy has to exist in memory.
        global_config = GEOSSimulator.create_global_config(
            args["windows"],
            args["gases"],
            args["TSIS"],
            73,
            72
        )

        # Let all workers have their own local copy
        @everywhere global_config = $global_config

        # Create the RetrievalToolbox buffer. THIS MUST BE DONE SEPARATELY ON ALL WORKERS
        # - OTHERWISE THE CRUCIAL OBJECT HASHING MECHANISM DOES NOT WORK!
        buffer = RESimulatorCore.create_buffer(global_config)
        @everywhere buffer = $buffer

        # By default, a SatelliteObserver is instantiated, but we can turn that into
        # an UplookingObserver, if needed.
        if args["uplooking"]
            @info "Running in uplooking mode!"
            @everywhere buffer.scene.observer = RE.UplookingGroundObserver()
        end
        #=
            Ingest coordinates to be sampled
        =#

        _raw_coord_txt = readlines(args["coords"])
        lonlat_array = zeros(2, length(_raw_coord_txt))
        for i in 1:length(_raw_coord_txt)
            lonlat_array[:,i] = parse.(Ref(Float64), split(_raw_coord_txt[i], ","))
        end

        # Create DateTime from args
        # (for now we use a single time for all obs for this run)

        dtime = DateTime(args["Y"], args["M"], args["D"], args["hour"], args["minute"])
        
        all_scenes = GEOSSimulator.generate_scenes_from_GEOS(
            global_config,
            GEOSIT_fdict,
            GEOSCARB_fdict,
            buffer,
            lonlat_array,
            dtime;
            NN=args["neighbors"]
        )

        @sync @everywhere @info "(synchronizing)"

        @info "Processing $(length(all_scenes)) scenes!"
        if nworkers() > 1
            @info "(parallel processing: $(nprocs()))"
            @time results = @showprogress showspeed=true @distributed (vcat) for scene in all_scenes
                [RESimulatorCore.process_scene!(buffer, scene)]
            end
        else
            results = []
            @showprogress showspeed=true for scene in all_scenes
                push!(results, RESimulatorCore.process_scene!(buffer, scene))
            end
        end

        GEOSSimulator.write_out(args["output"], global_config, all_scenes, results, buffer)

    end

end
