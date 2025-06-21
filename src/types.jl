import Base: String

mutable struct MagneticField
    kext::Int32
    options::Vector{Int32}
    sysaxes::Int32
end

function MagneticField(; options=[0, 0, 0, 0, 0], kext="OPQ77", sysaxes="GDZ")
    kext_val = parse_kext(kext)
    sysaxes_val = coord_sys(sysaxes)
    MagneticField(kext_val, options, sysaxes_val)
end

abstract type AbstractCoordinateSystem end

struct CoordinateVector{T,C} <: FieldVector{3,T}
    x::T
    y::T
    z::T
    sym::C
end

for sys in (:GDZ, :GEO, :GSM, :GSE, :SM, :GEI, :MAG, :SPH, :RLL, :HEE, :HAE, :HEEQ, :J2000)
    @eval struct $sys <: AbstractCoordinateSystem end
    @eval $sys(x, y, z) = CoordinateVector(x, y, z, $sys())
    @eval $sys(x) = CoordinateVector(x..., $sys())
    @eval export $sys
    @eval Base.String(::Type{$sys}) = $(String(sys))
end

@doc """Geocentric Solar Magnetospheric (GSM)\n\nX points sunward from Earth's center. The X-Z plane is defined to contain Earth's dipole axis (positive North).
""" GSM

coord(v::CoordinateVector) = v.sym
Base.String(::S) where {S<:AbstractCoordinateSystem} = String(S)
