"""
    MagneticField

A wrapper for IRBEM's magnetic field functions.

# Functions
- `make_lstar`: Compute magnetic coordinates at a spacecraft position
- `drift_shell`: Trace a full drift shell for particles with mirror point at input location
- `find_mirror_point`: Find magnitude and location of mirror point along field line
- `find_foot_point`: Find footprint of field line in a given hemisphere
- `trace_field_line`: Trace a full field line crossing the input position
- `find_magequator`: Find coordinates of magnetic equator from field line tracing
- `get_field_multi`: Compute GEO vector of magnetic field at input location
- `get_mlt`: Get Magnetic Local Time from GEO position and date
"""
module IRBEM

using Dates
using Printf
using LinearAlgebra
using IRBEM_jll

export MagneticField,
    make_lstar,
    drift_shell,
    drift_bounce_orbit,
    find_mirror_point,
    find_foot_point,
    trace_field_line,
    find_magequator,
    get_field_multi,
    get_mlt

# Physical constants
const Re = 6371.0  # Earth radius in km
const c = 3.0e8    # Speed of light in m/s

# External magnetic field model lookup table
const EXT_MODELS = ["None", "MF75", "TS87", "TL87", "T89", "OPQ77", "OPD88", "T96",
    "OM97", "T01", "T01S", "T04", "A00", "T07", "MT"]

include("magnetic_field.jl")
include("utils.jl")

end # module
