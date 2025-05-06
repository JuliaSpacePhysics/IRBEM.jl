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
- `Dict`: Contains keys Lm, MLT, Blocal, Bmin, Lstar, and XJ
"""
function make_lstar(model::MagneticField, X::Dict, maginput::Dict)
    # Process input coordinates and time
    ntime, iyear, idoy, ut, x1, x2, x3 = process_coords_time(X)

    # Check if ntime exceeds NTIME_MAX
    if ntime > NTIME_MAX[]
        throw(ArgumentError("Number of time steps ($ntime) exceeds NTIME_MAX ($NTIME_MAX[])"))
    end

    # Process magnetic field model inputs
    maginput_array = prepare_maginput(maginput, ntime)

    # Initialize output arrays
    Lm, Lstar, Blocal, Bmin, XJ, mlt = [zeros(Float64, ntime) for _ in 1:6]

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
        Blocal::Ptr{Float64}, Bmin::Ptr{Float64},
        XJ::Ptr{Float64}, mlt::Ptr{Float64}
    )::Cvoid

    # Return results as a dictionary
    return Dict(
        "Lm" => ntime == 1 ? Lm[1] : Lm,
        "MLT" => ntime == 1 ? mlt[1] : mlt,
        "Blocal" => ntime == 1 ? Blocal[1] : Blocal,
        "Bmin" => ntime == 1 ? Bmin[1] : Bmin,
        "Lstar" => ntime == 1 ? Lstar[1] : Lstar,
        "XJ" => ntime == 1 ? XJ[1] : XJ
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
    kext = model.kext
    options = model.options
    sysaxes = model.sysaxes

    @ccall libirbem.get_field_multi_(
        ntime::Ref{Int32},
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
