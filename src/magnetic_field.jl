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
    @ccall libirbem.make_lstar1_(
        Ref(Int32(ntime))::Ref{Int32},
        Ref(model.kext)::Ref{Int32},
        model.options::Ptr{Int32},
        Ref(model.sysaxes)::Ref{Int32},
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

# Returns
- `Dict`: Contains keys Lm, Lstar, blocal, bmin, xj, POSIT, and Nposit
"""
function drift_shell(model::MagneticField, X::Dict, maginput::Dict)
    # Process input coordinates and time
    ntime, iyear, idoy, ut, x1, x2, x3 = process_coords_time(X)

    # Check if ntime exceeds 1 (drift_shell only supports single time)
    if ntime > 1
        throw(ArgumentError("drift_shell only supports a single time point"))
    end

    # Process magnetic field model inputs
    maginput_array = prepare_maginput(maginput, ntime)

    # Initialize output arrays
    Lm = Ref{Float64}(0.0)
    Lstar = Ref{Float64}(0.0)
    blocal = Ref{Float64}(0.0)
    bmin = Ref{Float64}(0.0)
    xj = Ref{Float64}(0.0)

    # Maximum number of points in drift shell
    max_points = 1000
    posit = zeros(Float64, (3, max_points))
    nposit = Ref{Int32}(0)

    # Call IRBEM library function using @ccall
    @ccall libirbem.drift_shell1_(
        Ref(model.kext)::Ref{Int32},
        model.options::Ptr{Int32},
        Ref(model.sysaxes)::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        x1::Ptr{Float64}, x2::Ptr{Float64}, x3::Ptr{Float64},
        maginput_array::Ptr{Float64},
        Lm::Ref{Float64}, Lstar::Ref{Float64},
        blocal::Ref{Float64}, bmin::Ref{Float64},
        xj::Ref{Float64}, posit::Ptr{Float64}, nposit::Ref{Int32}
    )::Cvoid

    # Extract valid positions
    valid_posit = posit[:, 1:nposit[]]

    # Return results as a dictionary
    return Dict(
        "Lm" => Lm[],
        "Lstar" => Lstar[],
        "blocal" => blocal[],
        "bmin" => bmin[],
        "xj" => xj[],
        "POSIT" => valid_posit,
        "Nposit" => nposit[]
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
    blocal = Ref{Float64}(0.0)
    bmirror = Ref{Float64}(0.0)
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

    kext = model.kext
    options = model.options
    sysaxes = model.sysaxes

    # Call IRBEM library function using @ccall
    kext = Ref{Int32}(kext)
    sysaxes = Ref{Int32}(sysaxes)
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
    kext = Ref{Int32}(model.kext)
    options = model.options
    sysaxes = Ref{Int32}(model.sysaxes)

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
    xgeo = zeros(Float64, 3)

    # Call IRBEM library function using @ccall
    @ccall libirbem.find_magequator1_(
        Ref(model.kext)::Ref{Int32},
        model.options::Ptr{Int32},
        Ref(model.sysaxes)::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        x1::Ptr{Float64}, x2::Ptr{Float64}, x3::Ptr{Float64},
        maginput_array::Ptr{Float64},
        bmin::Ref{Float64},
        xgeo::Ptr{Float64}
    )::Cvoid

    # Return results as a dictionary
    return Dict(
        "bmin" => bmin[],
        "XGEO" => xgeo
    )
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
    kext = Ref{Int32}(model.kext)
    options = model.options
    sysaxes = Ref{Int32}(model.sysaxes)

    @ccall libirbem.get_field_multi_(
        ntime_ref::Ref{Int32}, kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
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
        Ref(iyr)::Ref{Int32},
        Ref(idoy)::Ref{Int32},
        Ref(ut)::Ref{Float64},
        xgeo::Ptr{Float64},
        mlt::Ref{Float64}
    )::Cvoid

    return mlt[]
end
