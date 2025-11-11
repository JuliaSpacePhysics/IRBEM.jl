# IRBEM.jl

[![Build Status](https://github.com/JuliaSpacePhysics/IRBEM.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JuliaSpacePhysics/IRBEM.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JuliaSpacePhysics/IRBEM.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaSpacePhysics/IRBEM.jl)

A Julia wrapper for the [IRBEM (International Radiation Belt Environment Modeling) Fortran library](https://prbem.github.io/IRBEM/).

**Installation**: at the Julia REPL, run `using Pkg; Pkg.add("IRBEM")`

**Documentation**: [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaSpacePhysics.github.io/IRBEM.jl/dev/)

## Overview

The IRBEM library is a set of source codes dedicated to radiation belt modeling. It facilitates the calculation of magnetic coordinates and drift shells using various external magnetic field models. This Julia package provides a wrapper around the IRBEM library, allowing Julia users to access its functionality.

## Features

- Computing magnetic field coordinates
- Find points of interest on the field line
- Compute magnetic field, derivatives and gradients
- Field tracing
- Coordinates transformations
- Flexible interface (Multiple dispatch with Fortran-style, Python-style, and Julia-style)

## Usage

```julia
using IRBEM
using Dates

# Initialize the magnetic field model
kext = T89
t = DateTime("2015-02-02T06:12:43")
𝐫 = GDZ(651, 63, 15.9)

# Define magnetic field model inputs
maginput = (; Kp = 40.0)

# Compute L* and related parameters
make_lstar(t, 𝐫, maginput; kext)

# Trace a field line
trace_field_line(t, 𝐫, maginput; kext)

# Find the magnetic equator
find_magequator(t, 𝐫, maginput; kext)

# Calculate MLT
get_mlt(t, 𝐫)
```

## Acknowledgments

- The IRBEM library is developed and maintained by the IRBEM-LIB development team.
- This Julia wrapper is inspired by the existing Python and MATLAB wrappers for IRBEM.
