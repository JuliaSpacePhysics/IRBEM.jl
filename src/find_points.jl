"""
https://prbem.github.io/IRBEM/api/magnetic_coordinates.html#points-of-interest-on-the-field-line
"""

"""
    find_mirror_point(time, x, alpha, coord="GDZ", maginput=Dict(); kext=KEXT[], options=OPTIONS[])
    find_mirror_point(model::MagneticField, X, alpha, maginput=Dict())

Find the magnitude and location of the mirror point along a field line traced from any given location and local pitch-angle.

# Arguments
$SIG_DOC
- `alpha`: Local pitch angle in degrees

# Outputs
- Blocal: magnitude of magnetic field at point (nT)
- Bmirr: magnitude of the magnetic field at the mirror point (nT)
- posit (array of 3 double): GEO coordinates of the mirror point (Re)

References: [IRBEM API](https://prbem.github.io/IRBEM/api/magnetic_coordinates.html#routine-FIND_MIRROR_POINT)
"""
find_mirror_point(arg1, arg2, alpha, args...; kw_args...) =
    find_mirror_point1(prepare_irbem(arg1, arg2, args...; kw_args...)[2:end]..., Float64(alpha))

function find_mirror_point1(args...; kw_args...)
    Blocal, Bmirr, posit = Ref{Float64}(), Ref{Float64}(), zeros(Float64, 3)
    find_mirror_point1!(args..., Blocal, Bmirr, posit)
    return (; Blocal=Blocal[], Bmirr=Bmirr[], posit)
end

"""
    find_foot_point(time, x, stop_alt, hemi_flag, coord="GDZ", maginput=Dict(); kext=KEXT[], options=OPTIONS[])
    find_foot_point(model::MagneticField, X, stop_alt, hemi_flag, maginput=Dict())

Find the footprint of a field line that passes through location X in a given hemisphere.

# Arguments
$SIG_DOC
- `stop_alt`: Altitude in km where to stop field line tracing
- `hemi_flag`: Hemisphere flag (0: same as SM z, +1: northern, -1: southern)

# Outputs
- XFOOT: GDZ coordinates of the foot point (Re)
- BFOOT: magnetic field vector (GEO) at the foot point (nT)
- BFOOTMAG: magnitude of the magnetic field at the foot point (nT)

References: [IRBEM API](https://prbem.github.io/IRBEM/api/magnetic_coordinates.html#routine-FIND_FOOT_POINT)
"""
find_foot_point(arg1, arg2, stop_alt, hemi_flag, args...; kw_args...) =
    find_foot_point1(prepare_irbem(arg1, arg2, args...; kw_args...)[2:end]..., Float64(stop_alt), Int32(hemi_flag))

function find_foot_point1(args...; kw_args...)
    BFOOTMAG, XFOOT, BFOOT = Ref{Float64}(), zeros(Float64, 3), zeros(Float64, 3)
    find_foot_point1!(args..., BFOOTMAG, XFOOT, BFOOT)
    return (; XFOOT, BFOOT, BFOOTMAG=BFOOTMAG[])
end


"""
    find_magequator($SIG1)
    find_magequator($SIG2)

Find the coordinates of the magnetic equator from tracing the magnetic field line from the input location.
Returns a named tuple with fields Bmin and XGEO (location of magnetic equator in GEO coordinates).

# Arguments
$SIG_DOC
"""
find_magequator(args...; kwargs...) = find_magequator1(prepare_irbem(args...; kwargs...)[2:end]...)

function find_magequator1(args...; kw_args...)
    Bmin, XGEO = Ref{Float64}(), zeros(Float64, 3)
    find_magequator1!(args..., Bmin, XGEO)
    return (; Bmin=Bmin[], XGEO)
end