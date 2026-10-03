# https://prbem.github.io/IRBEM/api/magnetic_coordinates.html#points-of-interest-on-the-field-line
# IRBEM only has single-point versions of these routines, so we loop over points.

"""
    find_mirror_point(time, x, alpha, [coord="GDZ",] maginput=(; ); kext=KEXT[], options=OPTIONS[])
    find_mirror_point(model::MagneticField, X, alpha, maginput=(; ))

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
function find_mirror_point(arg1, arg2, alpha, args...; kw...)
    p, single = prepare_irbem(arg1, arg2, args...; kw...)
    nt = (; Blocal = _out(p), Bmirr = _out(p), posit = _out(p, 3))
    return _output(_each_point!(find_mirror_point1!, p, (Float64(alpha),), nt), single)
end

"""
    find_foot_point(time, x, stop_alt, hemi_flag, [coord="GDZ",] maginput=(; ); kext=KEXT[], options=OPTIONS[])
    find_foot_point(model::MagneticField, X, stop_alt, hemi_flag, maginput=(; ))

Find the footprint of a field line that passes through location X in a given hemisphere.

# Arguments
$SIG_DOC
- `stop_alt`: Altitude in km where to stop field line tracing
- `hemi_flag`: Hemisphere flag (0: same as SM z, +1: northern, -1: southern)

# Outputs
- XFOOT: GDZ coordinates of the foot point (km, deg, deg)
- BFOOT: magnetic field vector (GEO) at the foot point (nT)
- BFOOTMAG: magnitude of the magnetic field at the foot point (nT)

References: [IRBEM API](https://prbem.github.io/IRBEM/api/magnetic_coordinates.html#routine-FIND_FOOT_POINT)
"""
function find_foot_point(arg1, arg2, stop_alt, hemi_flag, args...; kw...)
    p, single = prepare_irbem(arg1, arg2, args...; kw...)
    nt = (; XFOOT = _out(p, 3), BFOOT = _out(p, 3), BFOOTMAG = _out(p))
    return _output(_each_point!(find_foot_point1!, p, (Float64(stop_alt), Int32(hemi_flag)), nt), single)
end

"""
    find_magequator($SIG1)
    find_magequator($SIG2)

Find the coordinates of the magnetic equator from tracing the magnetic field line from the input location.
Returns a named tuple with fields Bmin and XGEO (location of magnetic equator).

# Arguments
$SIG_DOC
"""
function find_magequator(args...; kw...)
    p, single = prepare_irbem(args...; kw...)
    nt = (; Bmin = _out(p), XGEO = _out(p, 3))
    return _output(_each_point!(find_magequator1!, p, (), nt), single)
end
