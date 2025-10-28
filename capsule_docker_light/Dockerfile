# GRAPE capsule — Docker Light (Windows/Linux compatible)
FROM julia:1.11-bullseye

WORKDIR /grape
COPY . .

RUN julia --project=. -e "using Pkg; Pkg.instantiate()"

CMD ["julia", "--project=.", "examples/example_ParkerSolarProbe.jl"]
