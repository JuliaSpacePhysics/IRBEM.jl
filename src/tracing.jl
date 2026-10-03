# https://prbem.github.io/IRBEM/api/magnetic_coordinates.html#field-tracing
# These routines take a single point; output sizes are fixed by IRBEM.

"""
    trace_field_line($SIG1, R0=1.0)
    trace_field_line($SIG2; R0=1.0)

Trace a full field line which crosses the input position until radial distance `R0=1.0` (Re).

# Outputs
- Lm: L McIlwain
- Blocal (array of Nposit double): magnitude of magnetic field at point (nT)
- Bmin: magnitude of magnetic field at equator (nT)
- XJ: I, related to second adiabatic invariant (Re)
- posit (array of (3, Nposit) double): Cartesian coordinates in GEO along the field line
- Nposit: number of points in posit

Reference: [IRBEM API](https://prbem.github.io/IRBEM/api/magnetic_coordinates.html#routine-TRACE_FIELD_LINE)
"""
function trace_field_line(args...; R0 = 1.0, kw...)
    max_points = 3000
    posit = Matrix{Float64}(undef, 3, max_points)
    Blocal = Vector{Float64}(undef, max_points)
    Lm, Bmin, XJ, Nposit = Ref{Float64}(), Ref{Float64}(), Ref{Float64}(), Ref{Int32}()
    trace_field_line2_1_!(_single(prepare_irbem(args...; kw...))..., Float64(R0), Lm, Blocal, Bmin, XJ, posit, Nposit)
    N = Nposit[]
    nt = (; Lm, Blocal = Blocal[1:N], Bmin, XJ, posit = posit[:, 1:N], Nposit = N)
    return map(_nan, nt)
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
- Bmirr: magnitude of magnetic field at mirror point (nT)
- XJ: I, related to second adiabatic invariant (Re)
- posit (array of (3, 1000, 25) double): Cartesian coordinates in GEO along the drift shell
- Nposit (array of 25 integer): number of points in posit along each traced field line
- hmin, hmin_lon: GDZ altitude (km) and longitude (deg) of the lowest point of the drift shell, among all traced points

Entries of `posit` and `Blocal` beyond `Nposit` are `NaN`.

Reference: [IRBEM API](https://prbem.github.io/IRBEM/api/magnetic_coordinates.html#routine-DRIFT_BOUNCE_ORBIT)
"""
function drift_bounce_orbit(args...; alpha = 90, R0 = 1, kw...)
    max_points, n_azimuth = 1000, 25
    posit = Array{Float64, 3}(undef, 3, max_points, n_azimuth)
    Blocal = Matrix{Float64}(undef, max_points, n_azimuth)
    Nposit = zeros(Int32, n_azimuth)
    nt = (;
        Lm = Ref{Float64}(), Lstar = Ref{Float64}(), Blocal, Bmin = Ref{Float64}(), Bmirr = Ref{Float64}(),
        XJ = Ref{Float64}(), posit, Nposit, hmin = Ref{Float64}(), hmin_lon = Ref{Float64}(),
    )
    drift_bounce_orbit2_1_!(_single(prepare_irbem(args...; kw...))..., Float64(alpha), Float64(R0), nt...)
    clean_posit!(posit, Blocal, Nposit)
    return map(_nan, nt)
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

Entries of `posit` and `Blocal` beyond `Nposit` are `NaN`.

Reference: [IRBEM API](https://prbem.github.io/IRBEM/api/magnetic_coordinates.html#routine-DRIFT_SHELL)
"""
function drift_shell(args...; kw...)
    max_points, n_azimuth = 1000, 48
    posit = Array{Float64, 3}(undef, 3, max_points, n_azimuth)
    Blocal = Matrix{Float64}(undef, max_points, n_azimuth)
    Nposit = zeros(Int32, n_azimuth)
    nt = (; Lm = Ref{Float64}(), Lstar = Ref{Float64}(), Blocal, Bmin = Ref{Float64}(), XJ = Ref{Float64}(), posit, Nposit)
    drift_shell1_!(_single(prepare_irbem(args...; kw...))..., nt...)
    clean_posit!(posit, Blocal, Nposit)
    return map(_nan, nt)
end
