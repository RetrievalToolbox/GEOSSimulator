"""
"""
function create_GEOS_filedict(
    tdir::String,
    tstamp::String,
    pat_dict::Dict{String, Regex},
    )

    # Check if the directory exists:
    if !isdir(tdir)
        error("Directory $(tdir) does not exist!")
    end

    # Grab all files in the y/m/d directory
    flist_date = readdir(tdir)
    @info "N = $(length(flist_date)) files found in $(tdir)"

    # Make sure it is not empty!
    if length(flist_date) == 0
        error("No files found in: $(geos_dir). Check paths!")
    end

    @info "Reading from following files:"

    # Dictionary to hold the file names
    file_dict = Dict{String, String}()
    for (key, pat) in pat_dict

        # Pick out subset of files from `flist_date` where we see a regex match
        filt_result = filter(x -> !isnothing(match(pat, x)), flist_date)

        if length(filt_result) != 1
            error("Sorry, need exactly one result for pattern for [$(key) -> $(pat)]: \n" *
                "$(filt_result)"
                )
        end

        file_dict[key] = joinpath(tdir, filt_result[1])

        @info "[$(key)] => $(file_dict[key])"

    end

    # Return the dictionary
    return file_dict

end
