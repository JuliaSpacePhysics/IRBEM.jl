# https://prbem.github.io/IRBEM/api/magnetic_coordinates.html#field-tracing

"""
    trace_field_line($SIG1, R0=1.0)
    trace_field_line($SIG2; R0=1.0)

Trace a full field line which crosses the input position until radial distance `R0=1.0` (Re).

# Outputs
- Lm: L McIlwain
- Blocal (array of 3000 double): magnitude of magnetic field at point (nT)
- Bmin: magnitude of magnetic field at equator (nT)
- XJ: I, related to second adiabatic invariant (Re)
- posit (array of (3, 3000) double): Cartesian coordinates in GEO along the field line
- Nposit: number of points in posit

Reference: [IRBEM API](https://prbem.github.io/IRBEM/api/magnetic_coordinates.html#routine-TRACE_FIELD_LINE)
"""
function trace_field_line(args...; R0=1.0, kw...)
    # Prepare arguments
    max_points = 3000
    posit = zeros(Float64, 3, max_points)
    Nposit = Ref{Int32}(0)
    Blocal = zeros(Float64, max_points)
    nt = (; Lm = RF64(), Blocal, Bmin = RF64(), XJ = RF64(), posit, Nposit)

    trace_field_line2_1_!(prepare_irbem(args...; kw...)[2:end]..., Float64(R0), nt...)
    # Extract valid positions
    valid_posit = posit[:, 1:Nposit[]]
    return (map(_deref, nt)..., posit=valid_posit)
end

"""
    drift_bounce_orbit($SIG1, alpha=90, R0=1)
    drift_bounce_orbit($SIG2; alpha=90, R0=1)

Trace a full drift-bounce orbit for particles with a specified pitch angle `alpha=90` at the input location until radial distance `R0=1.0` (Re). 
Returns only positions between mirror points, with 25 azimuths.

# Outputs:
- Lm: L McIlwain
- Lstar: L Roederer or Φ=2π Bo/L* (nT Re2), depending on the options value
- Blocal (array of (1000, 25) double): magnitude of magnetic field at point (nT)
- Bmin: magnitude of magnetic field at equator (nT)
- XJ: I, related to second adiabatic invariant (Re)
- posit (array of (3, 1000, 25) double): Cartesian coordinates in GEO along the drift shell
- Nposit (array of 25 integer): number of points in posit along each traced field line

Reference: [IRBEM API](https://prbem.github.io/IRBEM/api/magnetic_coordinates.html#routine-DRIFT_BOUNCE_ORBIT)
"""
function drift_bounce_orbit(args...; alpha=90, R0=1, kw...)
    # Prepare arguments
    max_points = 1000
    n_azimuth = 25
    posit = Array{Float64,3}(undef, 3, max_points, n_azimuth)
    Nposit = zeros(Int32, n_azimuth)
    Blocal = zeros(Float64, max_points, n_azimuth)
    nt = (; Lm = RF64(), Lstar = RF64(), Blocal, Bmin = RF64(), Bmirr = RF64(), XJ = RF64(), posit, Nposit, hmin = RF64(), hmin_lon = RF64())

    drift_bounce_orbit2_1_!(prepare_irbem(args...; kw...)[2:end]..., Float64(alpha), Float64(R0), nt...)
    clean_posit!(posit, Nposit)
    return map(_deref, nt)
end


"""
    drift_shell($SIG1)
    drift_shell($SIG2)

Trace a full drift shell for particles that have their mirror point at the input location.

# Outputs
- `Lm`: L McIlwain
- `Lstar`: L Roederer or Φ=2π Bo/L* (nT Re2), depending on the options value
- `Blocal` (array of (1000, 48)): magnitude of magnetic field at point (nT)
- `Bmin`: magnitude of magnetic field at equator (nT)
- `XJ`: I, related to second adiabatic invariant (Re)
- `posit` (array of (3, 1000, 48)): Cartesian coordinates in GEO along the drift shell
- `Nposit` (array of 48 integer): number of points in posit along each traced field line

Reference: [IRBEM API](https://prbem.github.io/IRBEM/api/magnetic_coordinates.html#routine-DRIFT_SHELL)
"""
function drift_shell(args...; kw...)
    # Prepare arguments
    max_points = 1000
    n_azimuth = 48
    posit = zeros(Float64, 3, max_points, n_azimuth)
    Nposit = zeros(Int32, n_azimuth)
    Blocal = zeros(Float64, max_points, n_azimuth)
    nt = (; Lm = RF64(), Lstar = RF64(), Blocal, Bmin = RF64(), XJ = RF64(), posit, Nposit)

    drift_shell1_!(prepare_irbem(args...; kw...)[2:end]..., nt...)
    clean_posit!(posit, Nposit)
    return map(_deref, nt)
end