
# GRAPE Light Capsule (Julia installed)

Assumes you have Julia already installed and on path.

To reproduce the example:

```bash
julia --project=. -e "using Pkg; Pkg.instantiate()"

julia --project=. examples/example_simplified_no_tetrad.jl

julia --project=. examples/example_ParkerSolarProbe_NoTetrad.jl

