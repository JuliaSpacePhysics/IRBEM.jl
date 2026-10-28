"""
A wrapper for International Radiation Belt Environment Modeling (IRBEM) library

See the [Documentation](https://juliaspacephysics.github.io/IRBEM.jl/dev/) for more information.

# Functions

## Computing magnetic field coordinates
- [`make_lstar`](@ref): Compute magnetic coordinates at a spacecraft position
- [`landi2lstar`](@ref): Fast empirical L* for locally mirroring particles
- [`get_mlt`](@ref): Get Magnetic Local Time from GEO position and date
- [`get_hemi`](@ref): Get magnetic hemisphere of input location

## Points of interest on the field line
- [`find_mirror_point`](@ref): Find magnitude and location of mirror point along field line
- [`find_foot_point`](@ref): Find footprint of field line in a given hemisphere
- [`find_magequator`](@ref): Find coordinates of magnetic equator from field line tracing

## Magnetic field computation
- [`get_field_multi`](@ref): Compute GEO vector of magnetic field at input location
- [`get_bderivs`](@ref): Compute the magnetic field and its 1st-order derivatives at each input location

## Field tracing
- [`trace_field_line`](@ref): Trace a full field line crossing the input position
- [`drift_shell`](@ref): Trace a full drift shell for particles with mirror point at input location
- [`drift_bounce_orbit`](@ref): Trace a full bounce orbit for particles with mirror point at input location

## Coordinates transformations
- [`coord_trans`](@ref): Transform coordinates from one system to another

## Library information
- [`get_igrf_version`](@ref): Returns the version number of the IGRF model
- [`irbem_fortran_version`](@ref): Provides the repository version number of the fortran source code
- [`irbem_fortran_release`](@ref): Provides the repository release tag of the fortran source code

Outputs that IRBEM could not compute (its `-1e31` fill value) are `NaN`.
Only one IRBEM routine runs at a time: calls from multiple threads are serialized.

# References
- [IRBEM Documentation](https://prbem.github.io/IRBEM/)
- [IRBEM GitHub](https://github.com/PRBEM/IRBEM)
- [spacepy.irbempy](https://spacepy.github.io/autosummary/spacepy.irbempy.html)
"""
module IRBEM
using Dates
using IRBEM_jll
using StaticArrays
using TimeseriesUtilities: mapsamples, times

export MagInput
export make_lstar, landi2lstar, get_field_multi, get_mlt, get_hemi
export get_bderivs
export find_mirror_point, find_magequator, find_foot_point
export trace_field_line, drift_shell, drift_bounce_orbit
export coord_trans
export get_igrf_version, irbem_fortran_version, irbem_fortran_release
export MF75, TS87, TL87, T89, OPQ77, OPD88, T96, OM97, T01, T01S, T04, A00, T07, MT

include("types.jl")
include("lib.jl")
include("utils.jl")

const SIG = """time, x, maginput=(; ); coord="GDZ", kext=OPQ77, options=(0, 0, 0, 0, 0)"""
const SIG_DOC = """
- `time`: Date and time (DateTime, Vector{DateTime}, or String)
- `x`: Position as a 3-vector, a 3×n array, or a vector of 3-vectors (one per time)
- `coord` (optional): IRBEM name of the coordinate system of raw `x` (default: "GDZ"); not allowed for a `CoordinateVector`
- `kext` (optional): External field model (default: `OPQ77`)
- `options` (optional): IRBEM model options, 5 integers (default: all 0)
- `maginput`: Named tuple, dictionary or `MagInput` with magnetic field model inputs (optional).
  Values are scalars (shared by all times) or vectors (one per time); unset inputs are missing (`-9999`).
"""

include("magnetic_field.jl")
include("find_points.jl")
include("tracing.jl")
include("coordinates.jl")
include("info.jl")
include("workload.jl")

end
