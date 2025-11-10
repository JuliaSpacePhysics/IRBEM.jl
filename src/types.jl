import Base: String

mutable struct MagneticField
    kext::Int32
    options::Vector{Int32}
    sysaxes::Int32
end

"""
Magnetic field inputs

See [IRBEM Documentation](https://prbem.github.io/IRBEM/api/general_information.html#magnetic-field-inputs)
"""
@kwdef mutable struct MagInput
    Kp::Float64 = 0.0
    Dst::Float64 = 0.0
    Dsw::Float64 = 0.0 # solar wind density (cm-3)
    Vsw::Float64 = 0.0 # solar wind velocity (km/s)
    Pdyn::Float64 = 0.0 # solar wind dynamic pressure (nPa)
    ByIMF::Float64 = 0.0 # GSM y component of interplanetary magnetic field (nT)
    BzIMF::Float64 = 0.0 # GSM z component of interplanetary magnetic field (nT)
    G1::Float64 = 0.0
    G2::Float64 = 0.0
    G3::Float64 = 0.0
    W1::Float64 = 0.0
    W2::Float64 = 0.0
    W3::Float64 = 0.0
    W4::Float64 = 0.0
    W5::Float64 = 0.0
    W6::Float64 = 0.0
    AL::Float64 = 0.0 # auroral index
end

MagInput(nt) = MagInput(; nt...)
function MagInput(d::AbstractDict)
    out = MagInput()
    for (key, param) in d
        sym = Symbol(key)
        setfield!(out, sym, convert(Float64, param))
    end
    out
end

Base.unsafe_convert(::Type{Ptr{Float64}}, a::MagInput) =
    Base.unsafe_convert(Ptr{Float64}, pointer_from_objref(a))

const param_indices = fieldnames(MagInput)

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
    @eval $sys(x, y, z) = CoordinateVector(promote(x, y, z)..., $sys())
    @eval $sys(x) = (@assert length(x) == 3; CoordinateVector(x[1], x[2], x[3], $sys()))
    @eval export $sys
    @eval Base.String(::Type{$sys}) = $(String(sys))
end

@doc """Geocentric Solar Magnetospheric (GSM)\n\nX points sunward from Earth's center. The X-Z plane is defined to contain Earth's dipole axis (positive North).
""" GSM

coord(v::CoordinateVector) = v.sym
coord_sys(v::CoordinateVector) = coord_sys(v.sym)
Base.String(::S) where {S<:AbstractCoordinateSystem} = String(S)
