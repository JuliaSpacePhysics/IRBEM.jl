# IRBEM.jl

[![DOI](https://zenodo.org/badge/978838161.svg)](https://doi.org/10.5281/zenodo.17577161)

The [IRBEM (International Radiation Belt Environment Modeling) Fortran library](https://prbem.github.io/IRBEM/) is dedicated to radiation belt modeling. It facilitates the calculation of magnetic coordinates and drift shells using various external magnetic field models. This Julia package wraps it.

## Features

- Find points of interest on the field line
- Compute magnetic field, derivatives and gradients
- Field tracing
- Coordinates transformations
- Single points or arrays of points

## Usage

```julia
using Pkg; Pkg.add("IRBEM")
using IRBEM, Dates

kext = T89
t = DateTime("2015-02-02T06:12:43")
𝐫 = [651, 63, 15.9] # GDZ: alt [km], lat, lon [deg]

# Define magnetic field model inputs
maginput = (; Kp = 40.0)

# Compute L* and related parameters
make_lstar(t, 𝐫, maginput; kext)

trace_field_line(t, 𝐫, maginput; kext)

find_magequator(t, 𝐫, maginput; kext) # -> a namedtuple with fields Bmin and XGEO

get_mlt(t, [2.2, 2.8, 0.35]) # GEO [Re]

# Transform raw components between IRBEM coordinate systems
coord_trans("GDZ", "GSM", t, 𝐫)
```

### Coordinates

Positions are raw components (vector, 3×n matrix, vector of vectors or tuples) in the IRBEM coordinate system named by `coord`
(`"GDZ"` default, `"GEO"`, `"GSM"`, `"GSE"`, `"SM"`, `"GEI"`, `"MAG"`, `"SPH"`, `"RLL"`, `"HEE"`, `"HAE"`, `"HEEQ"`, `"TOD"`, `"J2000"`, `"TEME"`),
in [IRBEM's component order](https://prbem.github.io/IRBEM/api/general_information.html#coordinate-systems).

[GeoCotrans.jl](https://github.com/JuliaSpacePhysics/GeoCotrans.jl) `CoordinateVector`s are also accepted; their type gives the system, so omit `coord`.

## Acknowledgments

The IRBEM library is developed and maintained by the IRBEM-LIB development team.
