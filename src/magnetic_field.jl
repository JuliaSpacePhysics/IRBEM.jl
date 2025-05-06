mutable struct MagneticField
    kext::Int32
    options::Vector{Int32}
    sysaxes::Int32
    verbose::Bool
    NTIME_MAX::Int32

    function MagneticField(;
        options::Vector{Int}=[0, 0, 0, 0, 0],
        kext::Union{String,Int}="OPQ77",
        sysaxes::Union{String,Int}="GDZ",
        verbose::Bool=false
    )
        # Set kext (external magnetic field model)
        if isa(kext, String)
            try
                kext_val = findfirst(isequal(kext), EXT_MODELS) - 1
                if kext_val === nothing
                    throw(ArgumentError("Unknown external field model: $kext"))
                end
            catch
                throw(ArgumentError("Unknown external field model: $kext"))
            end
        else
            kext_val = kext
        end

        sysaxes_val = coord_sys(sysaxes)

        # Get NTIME_MAX from library
        ntime_max = Ref{Int32}(0)
        @ccall libirbem.get_irbem_ntime_max1_(ntime_max::Ref{Int32})::Cvoid

        new(
            Int32(kext_val),
            Int32.(options),
            sysaxes_val,
            verbose,
            ntime_max[]
        )
    end
end

"""
    make_lstar(model::MagneticField, X::Dict, maginput::Dict)

Compute magnetic coordinates at a spacecraft position.

# Arguments
- `model::MagneticField`: The magnetic field model
- `X::Dict`: Dictionary with keys:
  - `dateTime` or `Time`: Date and time (DateTime or String)
  - `x1`, `x2`, `x3`: Position coordinates in the system specified by `sysaxes`
- `maginput::Dict`: Dictionary with magnetic field model inputs

# Returns
- `Dict`: Contains keys Lm, MLT, blocal, bmin, Lstar, and xj
"""
function make_lstar(model::MagneticField, X::Dict, maginput::Dict)
    # Process input coordinates and time
    ntime, iyear, idoy, ut, x1, x2, x3 = process_coords_time(X)

    # Check if ntime exceeds NTIME_MAX
    if ntime > model.NTIME_MAX
        throw(ArgumentError("Number of time steps ($ntime) exceeds NTIME_MAX ($(model.NTIME_MAX))"))
    end

    # Process magnetic field model inputs
    maginput_array = prepare_maginput(maginput, ntime)

    # Initialize output arrays
    Lm, Lstar, blocal, bmin, xj, mlt = [zeros(Float64, ntime) for _ in 1:6]

    # Call IRBEM library function using @ccall
    sysaxes = model.sysaxes
    kext = model.kext
    options = model.options
    @ccall libirbem.make_lstar1_(
        ntime::Ref{Int32}, kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        x1::Ptr{Float64}, x2::Ptr{Float64}, x3::Ptr{Float64},
        maginput_array::Ptr{Float64},
        Lm::Ptr{Float64}, Lstar::Ptr{Float64},
        blocal::Ptr{Float64}, bmin::Ptr{Float64},
        xj::Ptr{Float64}, mlt::Ptr{Float64}
    )::Cvoid

    # Return results as a dictionary
    return Dict(
        "Lm" => ntime == 1 ? Lm[1] : Lm,
        "MLT" => ntime == 1 ? mlt[1] : mlt,
        "blocal" => ntime == 1 ? blocal[1] : blocal,
        "bmin" => ntime == 1 ? bmin[1] : bmin,
        "Lstar" => ntime == 1 ? Lstar[1] : Lstar,
        "xj" => ntime == 1 ? xj[1] : xj
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
    Lm = Ref{Float64}()
    Lstar = Ref{Float64}()
    Bmin = Ref{Float64}()
    XJ = Ref{Float64}()

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

"""
    find_mirror_point(model::MagneticField, X::Dict, maginput::Dict, alpha::Float64)

Find the magnitude and location of the mirror point along a field line traced from any given location and local pitch-angle.

# Arguments
- `model::MagneticField`: The magnetic field model
- `X::Dict`: Dictionary with keys:
  - `dateTime` or `Time`: Date and time (DateTime or String)
  - `x1`, `x2`, `x3`: Position coordinates in the system specified by `sysaxes`
- `maginput::Dict`: Dictionary with magnetic field model inputs
- `alpha::Float64`: Local pitch angle in degrees

# Returns
- `Dict`: Contains keys blocal, bmin, and POSIT (GEO coordinates of mirror point)
"""
function find_mirror_point(model::MagneticField, X::Dict, maginput::Dict, alpha::Float64)
    # Process input coordinates and time
    ntime, iyear, idoy, ut, x1, x2, x3 = process_coords_time(X)

    # Check if ntime exceeds 1 (find_mirror_point only supports single time)
    if ntime > 1
        throw(ArgumentError("find_mirror_point only supports a single time point"))
    end

    # Process magnetic field model inputs
    maginput_array = prepare_maginput(maginput, ntime)

    # Initialize output arrays
    blocal = Ref{Float64}()
    bmirror = Ref{Float64}()
    POSIT = zeros(Float64, 3)

    # Call IRBEM library function using @ccall
    kext = Ref{Int32}(model.kext)
    options = model.options
    sysaxes = Ref{Int32}(model.sysaxes)
    alpha = Ref{Float64}(alpha)

    @ccall libirbem.find_mirror_point1_(
        kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        x1::Ptr{Float64}, x2::Ptr{Float64}, x3::Ptr{Float64},
        alpha::Ref{Float64}, maginput_array::Ptr{Float64},
        blocal::Ref{Float64}, bmirror::Ref{Float64}, POSIT::Ptr{Float64}
    )::Cvoid

    (; blocal=blocal[], bmin=bmirror[], POSIT)
end

"""
    find_foot_point(model::MagneticField, X::Dict, maginput::Dict, stop_alt::Float64, hemi_flag::Int)

Find the footprint of a field line that passes through location X in a given hemisphere.

# Arguments
- `model::MagneticField`: The magnetic field model
- `X::Dict`: Dictionary with keys:
  - `dateTime` or `Time`: Date and time (DateTime or String)
  - `x1`, `x2`, `x3`: Position coordinates in the system specified by `sysaxes`
- `maginput::Dict`: Dictionary with magnetic field model inputs
- `stop_alt::Float64`: Altitude in km where to stop field line tracing
- `hemi_flag::Int`: Hemisphere flag (0: same as SM z, +1: northern, -1: southern)

# Returns
- `Dict`: Contains keys XFOOT, BFOOT, and BFOOTMAG
"""
function find_foot_point(model::MagneticField, X::Dict, maginput::Dict, stop_alt, hemi_flag::Int)
    # Process input coordinates and time
    ntime, iyear, idoy, ut, x1, x2, x3 = process_coords_time(X)

    # Check if ntime exceeds 1 (find_foot_point only supports single time)
    if ntime > 1
        throw(ArgumentError("find_foot_point only supports a single time point"))
    end

    # Process magnetic field model inputs
    maginput_array = prepare_maginput(maginput, ntime)

    # Initialize output arrays
    XFOOT = zeros(Float64, 3)
    BFOOT = zeros(Float64, 3)
    BFOOTMAG = Ref{Float64}(0.0)

    # Call IRBEM library function using @ccall
    kext = model.kext
    options = model.options
    sysaxes = model.sysaxes
    stop_alt = Ref{Float64}(stop_alt)
    hemi_flag = Ref{Int32}(hemi_flag)

    @ccall libirbem.find_foot_point1_(
        kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        x1::Ptr{Float64}, x2::Ptr{Float64}, x3::Ptr{Float64},
        stop_alt::Ref{Float64}, hemi_flag::Ref{Int32}, maginput_array::Ptr{Float64},
        XFOOT::Ptr{Float64}, BFOOT::Ptr{Float64}, BFOOTMAG::Ptr{Float64}
    )::Cvoid

    return (; XFOOT, BFOOT, BFOOTMAG=BFOOTMAG[])
end

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
- `Dict`: Contains keys Lm, blocal, bmin, xj, POSIT, and Nposit
"""
function trace_field_line(model::MagneticField, X::Dict, maginput::Dict; R0::Float64=1.0)
    # Process input coordinates and time
    ntime, iyear, idoy, ut, x1, x2, x3 = process_coords_time(X)

    # Check if ntime exceeds 1 (trace_field_line only supports single time)
    if ntime > 1
        throw(ArgumentError("trace_field_line only supports a single time point"))
    end

    # Process magnetic field model inputs
    maginput_array = prepare_maginput(maginput, ntime)

    # Initialize output arrays
    R0 = Ref{Float64}(R0)
    Lm = Ref{Float64}()
    blocal = Ref{Float64}()
    bmin = Ref{Float64}()
    xj = Ref{Float64}()

    # Maximum number of points in field line
    max_points = 1000
    posit = zeros(Float64, (3, max_points))
    nposit = Ref{Int32}(0)

    # Call IRBEM library function using @ccall
    kext = model.kext
    options = model.options
    sysaxes = model.sysaxes

    @ccall IRBEM_jll.libirbem.trace_field_line2_1_(
        kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        x1::Ptr{Float64}, x2::Ptr{Float64}, x3::Ptr{Float64},
        maginput_array::Ptr{Float64},
        R0::Ref{Float64}, Lm::Ref{Float64},
        blocal::Ref{Float64}, bmin::Ref{Float64},
        xj::Ref{Float64}, posit::Ptr{Float64}, nposit::Ref{Int32}
    )::Cvoid

    # Extract valid positions
    valid_posit = posit[:, 1:nposit[]]

    # Return results as a dictionary
    return Dict(
        "Lm" => Lm[],
        "blocal" => blocal[],
        "bmin" => bmin[],
        "xj" => xj[],
        "POSIT" => valid_posit,
        "Nposit" => nposit[]
    )
end

"""
    find_magequator(model::MagneticField, X::Dict, maginput::Dict)

Find the coordinates of the magnetic equator from tracing the magnetic field line from the input location.

# Arguments
- `model::MagneticField`: The magnetic field model
- `X::Dict`: Dictionary with keys:
  - `dateTime` or `Time`: Date and time (DateTime or String)
  - `x1`, `x2`, `x3`: Position coordinates in the system specified by `sysaxes`
- `maginput::Dict`: Dictionary with magnetic field model inputs

# Returns
- `Dict`: Contains keys bmin and XGEO (location of magnetic equator in GEO coordinates)
"""
function find_magequator(model::MagneticField, X::Dict, maginput::Dict)
    # Process input coordinates and time
    ntime, iyear, idoy, ut, x1, x2, x3 = process_coords_time(X)

    # Check if ntime exceeds 1 (find_magequator only supports single time)
    if ntime > 1
        throw(ArgumentError("find_magequator only supports a single time point"))
    end

    # Process magnetic field model inputs
    maginput_array = prepare_maginput(maginput, ntime)

    # Initialize output arrays
    bmin = Ref{Float64}(0.0)
    XGEO = zeros(Float64, 3)

    # Call IRBEM library function using @ccall
    kext = model.kext
    options = model.options
    sysaxes = model.sysaxes
    @ccall libirbem.find_magequator1_(
        kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        x1::Ptr{Float64}, x2::Ptr{Float64}, x3::Ptr{Float64},
        maginput_array::Ptr{Float64},
        bmin::Ref{Float64}, XGEO::Ptr{Float64}
    )::Cvoid

    # Return results as a dictionary
    return (; bmin=bmin[], XGEO)
end

"""
    get_field_multi(model::MagneticField, X::Dict, maginput::Dict)

Compute the GEO vector of the magnetic field at input location for a set of internal/external magnetic field.

# Arguments
- `model::MagneticField`: The magnetic field model
- `X::Dict`: Dictionary with keys:
  - `dateTime` or `Time`: Date and time (DateTime or String)
  - `x1`, `x2`, `x3`: Position coordinates in the system specified by `sysaxes`
- `maginput::Dict`: Dictionary with magnetic field model inputs

# Returns
- `Dict`: Contains keys Bgeo (GEO components of B field) and Bmag (magnitude of B field)
"""
function get_field_multi(model::MagneticField, X::Dict, maginput::Dict)
    # Process input coordinates and time
    ntime, iyear, idoy, ut, x1, x2, x3 = process_coords_time(X)

    # Process magnetic field model inputs
    maginput_array = prepare_maginput(maginput, ntime)

    # Initialize output arrays
    Bgeo = zeros(Float64, (3, ntime))
    Bmag = zeros(Float64, ntime)

    # Call IRBEM library function using @ccall
    ntime_ref = Ref{Int32}(ntime)
    kext = model.kext
    options = model.options
    sysaxes = model.sysaxes

    @ccall libirbem.get_field_multi_(
        ntime_ref::Ref{Int32},
        kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        x1::Ptr{Float64}, x2::Ptr{Float64}, x3::Ptr{Float64},
        maginput_array::Ptr{Float64}, Bgeo::Ptr{Float64}, Bmag::Ptr{Float64}
    )::Cvoid

    if ntime == 1
        (; Bgeo=Bgeo[:, 1], Bmag=Bmag[1])
    else
        (; Bgeo, Bmag)
    end
end

"""
    get_mlt(model::MagneticField, X::Dict)

Get Magnetic Local Time (MLT) from a Cartesian GEO position and date.

# Arguments
- `model::MagneticField`: The magnetic field model
- `X::Dict`: Dictionary with keys:
  - `dateTime` or `Time`: Date and time (DateTime or String)
  - `x1`, `x2`, `x3`: Position coordinates in GEO system

# Returns
- `Float64`: The MLT value (hours)
"""
function get_mlt(model::MagneticField, X::Dict)
    if model.verbose
        println("Running get_mlt")
    end

    # Process time
    dt = get_datetime(X)
    iyr = Int32(year(dt))
    idoy = Int32(dayofyear(dt))
    ut = Float64(hour(dt) * 3600 + minute(dt) * 60 + second(dt))

    # Get GEO coordinates
    xgeo = [Float64(X["x1"]), Float64(X["x2"]), Float64(X["x3"])]

    # Initialize output
    mlt = Ref{Float64}(0.0)

    # Call IRBEM library function using @ccall
    @ccall libirbem.get_mlt1_(
        iyr::Ref{Int32}, idoy::Ref{Int32}, ut::Ref{Float64},
        xgeo::Ptr{Float64},
        mlt::Ref{Float64}
    )::Cvoid

    return mlt[]
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
- `Dict`: Contains keys Lm, lstar, blocal, bmin, bmirr, xj, POSIT, Nposit, hmin, hmin_lon

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
    # Output arrays
    Lm = Ref{Float64}(0.0)
    Lstar = Ref{Float64}(0.0)
    Bmin = Ref{Float64}(0.0)
    Bmirr = Ref{Float64}(0.0)
    XJ = Ref{Float64}(0.0)
    hmin = Ref{Float64}(0.0)
    hmin_lon = Ref{Float64}(0.0)
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
