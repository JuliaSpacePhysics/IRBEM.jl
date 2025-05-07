# https://prbem.github.io/IRBEM/api/magnetic_coordinates.html#field-tracing

"""
    trace_field_line(model::MagneticField, X::Dict, maginput::Dict; R0::Float64=1.0)

Trace a full field line which crosses the input position.

# Arguments
- `model::MagneticField`: The magnetic field model
- `X::Dict`: Dictionary with keys:
  - `dateTime` or `Time`: Date and time (DateTime or String)
  - `x1`, `x2`, `x3`: Position coordinates in the system specified by `sysaxes`
- `maginput::Dict`: Dictionary with magnetic field model inputs
- `R0::Float64=1.0`: Radial distance in Re units where to stop field line tracing

# Returns
- `Dict`: Contains keys Lm, Blocal, Bmin, XJ, POSIT, and Nposit
"""
function trace_field_line(model::MagneticField, X::Dict, maginput::Dict; R0=1.0)
    # Process input coordinates and time
    ntime, iyear, idoy, ut, x1, x2, x3 = process_coords_time(X)

    # Check if ntime exceeds 1 (trace_field_line only supports single time)
    if ntime > 1
        throw(ArgumentError("trace_field_line only supports a single time point"))
    end

    # Process magnetic field model inputs
    maginput_array = prepare_maginput(maginput, ntime)

    # Initialize output arrays
    max_points = 3000
    R0 = Float64(R0)
    @init_refs Float64 Lm Bmin XJ
    @init_refs Int32 Nposit
    posit = zeros(Float64, (3, max_points))
    Blocal = zeros(Float64, max_points)

    # Call IRBEM library function using @ccall
    kext = model.kext
    options = model.options
    sysaxes = model.sysaxes

    @ccall IRBEM_jll.libirbem.trace_field_line2_1_(
        kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
        iyear::Ref{Int32}, idoy::Ref{Int32}, ut::Ref{Float64},
        x1::Ref{Float64}, x2::Ref{Float64}, x3::Ref{Float64},
        maginput_array::Ptr{Float64}, R0::Ref{Float64},
        Lm::Ref{Float64}, Blocal::Ptr{Float64}, Bmin::Ref{Float64},
        XJ::Ref{Float64}, posit::Ptr{Float64}, Nposit::Ref{Int32}
    )::Cvoid

    # Extract valid positions
    Nposit = Nposit[]
    valid_posit = posit[:, 1:Nposit]

    # Return results as a dictionary
    return (;
        Lm=Lm[], Blocal, Bmin=Bmin[],
        XJ=XJ[], posit=valid_posit, Nposit
    )
end

"""
    drift_bounce_orbit(model::MagneticField, X::Dict, maginput::Dict; alpha::Float64=90.0, R0::Float64=1.0)

Trace a full drift-bounce orbit for particles with a specified pitch angle at the input location. Returns only positions between mirror points, with 25 azimuths.

# Arguments
- `model::MagneticField`: The magnetic field model
- `X::Dict`: Dictionary with keys:
  - `dateTime` or `Time`: Date and time (DateTime or String)
  - `x1`, `x2`, `x3`: Position coordinates in the system specified by `sysaxes`
- `maginput::Dict`: Dictionary with magnetic field model inputs
- `alpha::Float64`: Local pitch angle in degrees (default 90)
- `R0::Float64`: Minimum radial distance allowed along the drift path (default 1.0)

# Returns
- `Dict`: Contains keys Lm, lstar, Blocal, Bmin, bmirr, XJ, POSIT, Nposit, hmin, hmin_lon

Reference: [IRBEM API](https://prbem.github.io/IRBEM/api/magnetic_coordinates.html#routine-DRIFT_BOUNCE_ORBIT)
"""
function drift_bounce_orbit(model::MagneticField, X::Dict, maginput::Dict; alpha=90, R0=1)
    # Process input coordinates and time
    ntime, iyear, idoy, ut, x1, x2, x3 = process_coords_time(X)
    if ntime > 1
        throw(ArgumentError("drift_bounce_orbit only supports a single time point"))
    end
    maginput_array = prepare_maginput(maginput, ntime)

    # Prepare arguments
    alpha = Float64(alpha)
    R0 = Float64(R0)
    @init_refs Float64 Lm Lstar Bmin Bmirr XJ hmin hmin_lon
    # Output arrays
    Blocal = zeros(Float64, 1000, 25)
    posit = Array{Float64,3}(undef, 3, 1000, 25)
    Nposit = zeros(Int32, 25)
    # Call IRBEM library function
    kext = model.kext
    options = model.options
    sysaxes = model.sysaxes
    @ccall libirbem.drift_bounce_orbit2_1_(
        kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
        iyear::Ref{Int32}, idoy::Ref{Int32}, ut::Ref{Float64},
        x1::Ref{Float64}, x2::Ref{Float64}, x3::Ref{Float64},
        alpha::Ref{Float64}, maginput_array::Ptr{Float64},
        R0::Ref{Float64}, Lm::Ref{Float64}, Lstar::Ref{Float64},
        Blocal::Ptr{Float64}, Bmin::Ref{Float64}, Bmirr::Ref{Float64},
        XJ::Ref{Float64}, posit::Ptr{Float64}, Nposit::Ptr{Int32},
        hmin::Ref{Float64}, hmin_lon::Ref{Float64}
    )::Cvoid

    clean_posit!(posit, Nposit)

    return (;
        Lm=Lm[], Lstar=Lstar[],
        Blocal, Bmin, Bmirr,
        XJ=XJ[], posit, Nposit,
        hmin=hmin[], hmin_lon=hmin_lon[]
    )
end


"""
    drift_shell(model::MagneticField, X::Dict, maginput::Dict)

Trace a full drift shell for particles that have their mirror point at the input location.

# Arguments
- `model::MagneticField`: The magnetic field model
- `X::Dict`: Dictionary with keys:
  - `dateTime` or `Time`: Date and time (DateTime or String)
  - `x1`, `x2`, `x3`: Position coordinates in the system specified by `sysaxes`
- `maginput::Dict`: Dictionary with magnetic field model inputs

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
function drift_shell(model::MagneticField, X::Dict, maginput::Dict)
    # Process input coordinates and time
    ntime, iyear, idoy, ut, x1, x2, x3 = process_coords_time(X)
    if ntime > 1
        throw(ArgumentError("drift_shell only supports a single time point"))
    end
    maginput_array = prepare_maginput(maginput, ntime)

    # Output arrays (match Python shapes)
    max_points = 1000
    n_azimuth = 48
    posit = zeros(Float64, 3, max_points, n_azimuth)
    Blocal = zeros(Float64, max_points, n_azimuth)
    Nposit = zeros(Int32, n_azimuth)
    @init_refs Float64 Lm Lstar Bmin XJ

    # Call IRBEM library function using @ccall
    kext = model.kext
    options = model.options
    sysaxes = model.sysaxes
    @ccall libirbem.drift_shell1_(
        kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
        iyear::Ref{Int32}, idoy::Ref{Int32}, ut::Ref{Float64},
        x1::Ref{Float64}, x2::Ref{Float64}, x3::Ref{Float64},
        maginput_array::Ptr{Float64},
        Lm::Ref{Float64}, Lstar::Ref{Float64},
        Blocal::Ptr{Float64}, Bmin::Ref{Float64},
        XJ::Ref{Float64}, posit::Ptr{Float64}, Nposit::Ptr{Int32}
    )::Cvoid

    clean_posit!(posit, Nposit)

    return (;
        Lm=Lm[], Lstar=Lstar[],
        Blocal, Bmin=Bmin[],
        XJ=XJ[], posit, Nposit
    )
end