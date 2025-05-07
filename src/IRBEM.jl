"""
A wrapper for International Radiation Belt Environment Modeling (IRBEM) library

# Functions

## Computing magnetic field coordinates
- `make_lstar`: Compute magnetic coordinates at a spacecraft position

## Points of interest on the field line
- `find_mirror_point`: Find magnitude and location of mirror point along field line
- `find_foot_point`: Find footprint of field line in a given hemisphere
- `trace_field_line`: Trace a full field line crossing the input position
- `find_magequator`: Find coordinates of magnetic equator from field line tracing

## Magnetic field computation
- `get_field_multi`: Compute GEO vector of magnetic field at input location
- `get_mlt`: Get Magnetic Local Time from GEO position and date

## Field tracing
- `drift_shell`: Trace a full drift shell for particles with mirror point at input location
- `drift_bounce_orbit`: Trace a full bounce orbit for particles with mirror point at input location

## Coordinates transformations
- `transform`: Transform coordinates from one system to another

# References
- [IRBEM Documentation](https://prbem.github.io/IRBEM/)
- [IRBEM GitHub](https://github.com/PRBEM/IRBEM)
"""
module IRBEM
using Dates
using IRBEM_jll

export MagneticField
export make_lstar, get_field_multi, get_mlt
export find_mirror_point, find_magequator, find_foot_point
export trace_field_line, drift_shell, drift_bounce_orbit
export transform

const NTIME_MAX = Ref{Int32}()

# External magnetic field model lookup table
const EXT_MODELS = ["None", "MF75", "TS87", "TL87", "T89", "OPQ77", "OPD88", "T96",
    "OM97", "T01", "T01S", "T04", "A00", "T07", "MT"]

function __init__()
    @ccall libirbem.get_irbem_ntime_max1_(NTIME_MAX::Ref{Int32})::Cvoid
end

include("lib.jl")
include("utils.jl")
include("types.jl")
include("magnetic_field.jl")
include("find_points.jl")
include("tracing.jl")
include("coordinates.jl")

end # module
