"""
External Magnetic Field Model

See [IRBEM Documentation](https://prbem.github.io/IRBEM/api/general_information.html#kext) for more details.
"""
@enum ExternalFieldModel begin
    None = 0
    MF75 = 1
    TS87 = 2
    TL87 = 3
    T89 = 4
    OPQ77 = 5
    OPD88 = 6
    T96 = 7
    OM97 = 8
    T01 = 9
    T01S = 10
    T04 = 11
    A00 = 12
    T07 = 13
    MT = 14 # Mead-Tsyganenko
end

@doc "1: Mead & Fairfield [1975], uses 0 ≤ Kp ≤ 9, valid for rGEO ≤17 Re" MF75
@doc "2: Tsyganenko short [1987], uses 0 ≤ Kp ≤ 9, valid for rGEO ≤30 Re" TS87
@doc "3: Tsyganenko long [1987], uses 0 ≤ Kp ≤ 9, valid for rGEO ≤70 Re" TL87
@doc "4: Tsyganenko [1989], uses 0 ≤ Kp ≤ 9, valid for rGEO ≤70 Re" T89
@doc "5: Olson & Pfitzer quiet [1977], valid for rGEO ≤15 Re" OPQ77
@doc "6: Olson & Pfitzer dynamic [1988]" OPD88
@doc "7: Tsyganenko [1996]" T96
@doc "8: Ostapenko & Maltsev [1997]" OM97
@doc "9: Tsyganenko [2001], uses -50 ≤ Dst ≤ 20, 0.5 ≤ Pdyn ≤ 5, |By| ≤ 5, |Bz| ≤ 5, 0 ≤ G1 ≤ 10, 0 ≤ G2 ≤ 10, valid for xGSM ≥-15 Re" T01
@doc "10: Tsyganenko [2001], storm time, uses Dst, Pdyn, By, Bz, G2, G3, valid for xGSM ≥-15 Re" T01S
@doc "11: Tsyganenko [2004], uses Dst, Pdyn, By, Bz, W1, W2, W3, W4, W5, W6, valid for xGSM ≥-15 Re" T04
@doc "12: Alexeev [2000], also known as Paraboloid model, uses Dsw, Vsw, Dst, Bz, AL" A00
@doc "13: Tsyganenko [2007]" T07
@doc "14: Mead-Tsyganenko, uses Kp, onera model where the Tsyganenko 89 model is best fitted by a Mead model" MT

Base.convert(::Type{Int32}, x::ExternalFieldModel) = Int32(x)

mutable struct MagneticField
    kext::Int32
    options::Vector{Int32}
    sysaxes::Int32
end

function MagneticField(; options = [0, 0, 0, 0, 0], kext = OPQ77, sysaxes = "GDZ")
    return MagneticField(parse_kext(kext), options, coord_sys(sysaxes))
end

"""
Magnetic field inputs. Unset fields are `-9999`, which IRBEM treats as missing.

See [IRBEM Documentation](https://prbem.github.io/IRBEM/api/general_information.html#magnetic-field-inputs)
"""
@kwdef mutable struct MagInput
    Kp::Float64 = -9999.0 # Kp × 10 (e.g. 40 for Kp = 4)
    Dst::Float64 = -9999.0
    Dsw::Float64 = -9999.0 # solar wind density (cm-3)
    Vsw::Float64 = -9999.0 # solar wind velocity (km/s)
    Pdyn::Float64 = -9999.0 # solar wind dynamic pressure (nPa)
    ByIMF::Float64 = -9999.0 # GSM y component of interplanetary magnetic field (nT)
    BzIMF::Float64 = -9999.0 # GSM z component of interplanetary magnetic field (nT)
    G1::Float64 = -9999.0
    G2::Float64 = -9999.0
    G3::Float64 = -9999.0
    W1::Float64 = -9999.0
    W2::Float64 = -9999.0
    W3::Float64 = -9999.0
    W4::Float64 = -9999.0
    W5::Float64 = -9999.0
    W6::Float64 = -9999.0
    AL::Float64 = -9999.0 # auroral index
end

MagInput(nt) = MagInput(; nt...)
MagInput(d::AbstractDict) = MagInput(; (Symbol(k) => v for (k, v) in d)...)

abstract type AbstractCoordinateSystem end

struct CoordinateVector{T, C} <: FieldVector{3, T}
    x::T
    y::T
    z::T
    sym::C
end

# Arithmetic results (e.g. `a - b`) are not positions in a coordinate system
StaticArrays.similar_type(::Type{<:CoordinateVector{T}}, ::Type{T}, s::Size) where {T} = similar_type(SVector{3, T}, T, s)

# Order gives IRBEM's `sysaxes` code (0-based)
const COORD_SYSTEMS = (:GDZ, :GEO, :GSM, :GSE, :SM, :GEI, :MAG, :SPH, :RLL, :HEE, :HAE, :HEEQ, :TOD, :J2000, :TEME)

for sys in COORD_SYSTEMS
    @eval struct $sys <: AbstractCoordinateSystem end
    @eval $sys(x, y, z) = CoordinateVector(promote(x, y, z)..., $sys())
    @eval $sys(x) = (length(x) == 3 || throw(DimensionMismatch("expected 3 components, got $(length(x))")); $sys(x[1], x[2], x[3]))
    @eval export $sys
end

@doc "Geodetic (altitude, latitude, east longitude - km, deg, deg)" GDZ
@doc "Geocentric Geographic (Earth radii)" GEO
@doc "Geocentric Solar Magnetospheric (GSM)\n\nX points sunward from Earth's center. The X-Z plane is defined to contain Earth's dipole axis (positive North)." GSM
@doc "Geocentric Solar Ecliptic (GSE) (Earth radii)" GSE
@doc "Solar Magnetic (SM)" SM
@doc "Geocentric Equatorial Inertial (GEI), true of date (Earth radii)" GEI
@doc "Geomagnetic (MAG)" MAG
@doc "Spherical GEO (SPH) (radial distance, latitude, east longitude - Earth radii, deg, deg)" SPH
@doc "Geodetic (radial distance, latitude, East longitude - Earth radii, deg, deg)\n\nA re-expression of [`GDZ`](@ref) coordinates using radial distance instead of altitude above the reference ellipsoid." RLL
@doc "Heliocentric Earth Ecliptic (HEE) (Earth radii)" HEE
@doc "Heliocentric Aries Ecliptic (HAE) (Earth radii)" HAE
@doc "Heliocentric Earth Equatorial (HEEQ) (Earth radii)" HEEQ
@doc "True of Date (TOD), same as [`GEI`](@ref) (Earth radii)" TOD
@doc "[`GEI`](@ref) at J2000: mean equator and mean equinox of J2000 (Earth radii)" J2000
@doc "True Equator Mean Equinox (TEME), the inertial system of SGP4 (Earth radii)" TEME
