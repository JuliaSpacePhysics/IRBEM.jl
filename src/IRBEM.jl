"""
A wrapper for International Radiation Belt Environment Modeling (IRBEM) library

See the [Documentation](https://juliaspacephysics.github.io/IRBEM.jl/dev/) for more information.

# Functions

## Computing magnetic field coordinates
- [`make_lstar`](@ref): Compute magnetic coordinates at a spacecraft position
- [`get_mlt`](@ref): Get Magnetic Local Time from GEO position and date

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
- [`transform`](@ref): Transform coordinates from one system to another

## Library information
- [`get_igrf_version`](@ref): Returns the version number of the IGRF model
- [`irbem_fortran_version`](@ref): Provides the repository version number of the fortran source code
- [`irbem_fortran_release`](@ref): Provides the repository release tag of the fortran source code

# References
- [IRBEM Documentation](https://prbem.github.io/IRBEM/)
- [IRBEM GitHub](https://github.com/PRBEM/IRBEM)
- [spacepy.irbempy](https://spacepy.github.io/autosummary/spacepy.irbempy.html)
"""
module IRBEM
using Dates
using IRBEM_jll
using StaticArrays

export MagInput, MagneticField
export make_lstar, get_field_multi, get_mlt
export get_bderivs
export find_mirror_point, find_magequator, find_foot_point
export trace_field_line, drift_shell, drift_bounce_orbit
export transform
export get_igrf_version, irbem_fortran_version, irbem_fortran_release
export MF75, TS87, TL87, T89, OPQ77, OPD88, T96, OM97, T01, T01S, T04, A00, T07, MT

include("types.jl")
include("lib.jl")
include("utils.jl")
include("const.jl")

const KEXT = Ref{Int32}(OPQ77)
const OPTIONS = Ref{Vector{Int32}}([0, 0, 0, 0, 0])

const SIG1 = """time, x, [coord="GDZ",] maginput=(; ); kext=KEXT[], options=OPTIONS[]"""
const SIG2 = """model::MagneticField, X, maginput=(; )"""
const SIG_DOC = """
## Signature 1 (preferred):
- `time`: Date and time (DateTime, Vector{DateTime}, or String)
- `x`:  Position coordinates as a 3×n array or a tuple/vector of vectors.
      If the element type of `x` is `CoordinateVector`, `coord` is not needed.
- `coord` (optional): String specifying the coordinate system (default: "GDZ")
- `kext` (optional): External field model selection (default: KEXT[])
- `options` (optional): Model options (default: OPTIONS[])

## Signature 2:
- `model::MagneticField`: The magnetic field model
- `X`: Dictionary with keys:
  - `dateTime` or `Time`: Date and time (DateTime or String)
  - `x1`, `x2`, `x3`: Position coordinates in the system specified by `sysaxes`

## Common arguments:
- `maginput`: Named tuple or dictionary or `MagInput` with magnetic field model inputs (optional)
"""

include("magnetic_field.jl")
include("find_points.jl")
include("tracing.jl")
include("coordinates.jl")
include("info.jl")
include("python.jl")
include("workload.jl")

end
