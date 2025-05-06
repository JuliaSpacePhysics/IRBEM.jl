# IRBEM.jl

A Julia wrapper for the IRBEM (International Radiation Belt Environment Modeling) Fortran library.

## Overview

The IRBEM library is a set of source codes dedicated to radiation belt modeling. It facilitates the calculation of magnetic coordinates and drift shells using various external magnetic field models. This Julia package provides a wrapper around the IRBEM library, allowing Julia users to access its functionality.

## Installation

```julia
using Pkg
Pkg.add("IRBEM")
```

## Features

- Compute magnetic coordinates at spacecraft positions
- Trace drift shells for particles
- Find mirror points along field lines
- Find footprints of field lines
- Trace complete field lines
- Find magnetic equator coordinates
- Compute magnetic field vectors
- Calculate Magnetic Local Time (MLT)

## Usage

```julia
using IRBEM
using Dates

# Initialize the magnetic field model
model = MagneticField(verbose=false, kext="T89")

# Define position and time
X = Dict(
    "dateTime" => DateTime("2015-02-02T06:12:43"),
    "x1" => 600.0,  # km
    "x2" => 60.0,   # lat
    "x3" => 50.0    # lon
)

# Define magnetic field model inputs
maginput = Dict("Kp" => 40.0)

# Compute L* and related parameters
result = make_lstar(model, X, maginput)
println("L* = ", result["Lstar"])
println("Lm = ", result["Lm"])
println("Blocal = ", result["blocal"], " nT")
println("Bmin = ", result["bmin"], " nT")

# Trace a field line
field_line = trace_field_line(model, X, maginput)
println("Number of points in field line: ", field_line["Nposit"])

# Find the magnetic equator
equator = find_magequator(model, X, maginput)
println("Bmin at equator = ", equator["bmin"], " nT")
println("Equator position (GEO) = ", equator["XGEO"])

# Calculate MLT
mlt = get_mlt(model, X)
println("MLT = ", mlt, " hours")
```

## Building from Source

This package uses BinaryBuilder to provide pre-built binaries for the IRBEM library. If you want to build the library from source, you can use the provided `build_tarballs.jl` script:

```bash
julia --project=. build_tarballs.jl
```

## Coordinate Systems

The following coordinate systems are supported:

- `GDZ`: Geodetic (altitude, latitude, east longitude - km, deg, deg)
- `GEO`: Cartesian GEO (Earth radii)
- `GSM`: Cartesian GSM (Earth radii)
- `GSE`: Cartesian GSE (Earth radii)
- `SM`: Cartesian SM (Earth radii)
- `GEI`: Cartesian GEI (Earth radii)
- `MAG`: Cartesian MAG (Earth radii)
- `SPH`: Spherical GEO (radial distance, latitude, east longitude - Earth radii, deg, deg)
- `RLL`: Spherical GEO (radial distance, latitude, east longitude - Earth radii, deg, deg)

## External Magnetic Field Models

The following external magnetic field models are supported:

- `None`: No external field
- `MF75`: Mead & Fairfield 1975
- `TS87`: Tsyganenko 1987 short
- `TL87`: Tsyganenko 1987 long
- `T89`: Tsyganenko 1989
- `OPQ77`: Olson & Pfitzer quiet 1977
- `OPD88`: Olson & Pfitzer dynamic 1988
- `T96`: Tsyganenko 1996
- `OM97`: Ostapenko & Maltsev 1997
- `T01`: Tsyganenko 2001
- `T01S`: Tsyganenko 2001 storm
- `T04`: Tsyganenko 2004
- `A00`: Alexeev 2000
- `T07`: Tsyganenko 2007
- `MT`: Mead & Tsyganenko

## License

This package is distributed under the GNU Lesser General Public License (LGPL) v3.0, the same as the original IRBEM library.

## Acknowledgments

- The IRBEM library is developed and maintained by the IRBEM-LIB development team.
- This Julia wrapper is inspired by the existing Python and MATLAB wrappers for IRBEM.
