echo "### Deploying GEOSSimulator"


juliaup -V # Check if JuliaUp is available at all
if [ $? -ne 0 ]; then
    echo "### JuliaUp not found - installing"
    curl -fsSL https://install.julialang.org | sh
    echo "### Julia was installed. Please re-start or re-source your shell!"
    exit
fi

# Install Julia packages
echo "### Installing required Julia packages."
julia --project="./" -e 'using Pkg; Pkg.instantiate();'
