# IRBEM.jl

[![DOI](https://zenodo.org/badge/978838161.svg)](https://doi.org/10.5281/zenodo.17577161)
[![Coverage](https://codecov.io/gh/JuliaSpacePhysics/IRBEM.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSpacePhysics/IRBEM.jl)

A Julia wrapper for the [IRBEM (International Radiation Belt Environment Modeling) Fortran library](https://prbem.github.io/IRBEM/).

**Documentation**: [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSpacePhysics.github.io/IRBEM.jl/dev/)

## Overview

The IRBEM library is a set of source codes dedicated to radiation belt modeling. It facilitates the calculation of magnetic coordinates and drift shells using various external magnetic field models. This Julia package provides a wrapper around the IRBEM library, allowing Julia users to access its functionality.

## Features

- Computing magnetic field coordinates
- Find points of interest on the field line
- Compute magnetic field, derivatives and gradients
- Field tracing
- Coordinates transformations
- Flexible interface (Julia-style and Python-style `Dict` inputs, single points or arrays)
- Thread-safe: concurrent calls are serialized around the (stateful) Fortran library

## Usage

```julia
using Pkg; Pkg.add("IRBEM")
using IRBEM
using IRBEM: transform, get_mlt, find_magequator # not exported: GeoCotrans exports its own
using Dates

# Initialize the magnetic field model
kext = T89
t = DateTime("2015-02-02T06:12:43")
𝐫 = GDZ(63, 15.9, 651) # (lat [deg], lon [deg], alt [km])

# Define magnetic field model inputs
maginput = (; Kp = 40.0)

# Compute L* and related parameters
make_lstar(t, 𝐫, maginput; kext)

# Trace a field line
trace_field_line(t, 𝐫, maginput; kext)

find_magequator(t, 𝐫, maginput; kext) # -> a namedtuple with fields Bmin and XGEO

# Calculate MLT
get_mlt(t, 𝐫)

# Transform coordinates
transform(GDZ => GSM, 𝐫, t)
transform(GSM, 𝐫, t)                         # -> GSM CoordinateVector
transform("GEO" => "J2000", rand(3, 10), t)  # 3×n (or a vector of 3-vectors) with one time or one per sample; `dim = 1` for n×3
# transform(GEO => GSM, A) for a time-keyed container A (e.g. a Ti×3 DimArray) keeps its dims and metadata
```

### Coordinates

Positions are [GeoCotrans.jl](https://github.com/JuliaSpacePhysics/GeoCotrans.jl) `CoordinateVector`s, and IRBEM reuses its frames,
so `using GeoCotrans, IRBEM` does not clash:

- Frames: `GEO`, `GSM`, `GSE`, `SM`, `GEI`, `MAG`, `J2000` (GeoCotrans) and `HEE`, `HAE`, `HEEQ`, `TEME` (IRBEM only); `GSM(x, y, z)` etc. build a Cartesian position.
- `GDZ(lat, lon, alt)`: geodetic GEO; `CoordinateVector{GEO, Spherical}(r, θ, φ)`: spherical, colatitude and longitude in radians.
- `CoordinateVector`s carry their system: `transform(GDZ, x, t)` or `transform((GEO(), Spherical()), x, t)` return them in these conventions.
- Raw components (vectors, matrices) are described by a frame (Cartesian) or an IRBEM `sysaxes` name (`"GDZ"`, `:sph`, case-insensitive) or code, in IRBEM's component order:
  `GDZ` is `(alt, lat, lon)`, `SPH`/`RLL` are `(r, lat, lon)` in degrees. So are the `coord` arguments of `make_lstar` & co. and `MagneticField(; sysaxes)`;
  `RLL` and `TOD` exist only as names. Typed representations (`GDZ`, `(frame, representation)`) of raw components throw, as they would contradict the names' order.

## Acknowledgments

- The IRBEM library is developed and maintained by the IRBEM-LIB development team.
- This Julia wrapper is inspired by the existing Python and MATLAB wrappers for IRBEM.
