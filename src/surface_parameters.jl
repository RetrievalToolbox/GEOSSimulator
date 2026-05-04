function create_surface_parameters(spectral_windows::Vector{RE.SpectralWindow})
    return [(0.25,) for swin in spectral_windows]
end