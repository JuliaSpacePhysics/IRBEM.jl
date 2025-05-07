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
- `NamedTuple`: Contains fields Lm, MLT, Blocal, Bmin, Lstar, and XJ
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
    Lm, Lstar, Blocal, Bmin, XJ, MLT = [zeros(Float64, ntime) for _ in 1:6]

    # Call IRBEM library function using @ccall
    sysaxes = model.sysaxes
    kext = model.kext
    options = model.options
    make_lstar1!(
        ntime, kext, options, sysaxes, iyear, idoy, ut, x1, x2, x3, maginput_array, # inputs
        Lm, Lstar, Blocal, Bmin, XJ, MLT # outputs
    )
    return map(_only, (; Lm, MLT, Blocal, Bmin, Lstar, XJ))
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
- `NamedTuple`: Contains fields Bgeo (GEO components of B field) and Bmag (magnitude of B field)
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
    get_field_multi!(
        ntime, kext, options, sysaxes, iyear, idoy, ut, x1, x2, x3, maginput_array, # inputs
        Bgeo, Bmag # outputs
    )
    return map(_only, (; Bgeo, Bmag))
end

"""
    get_mlt(model::MagneticField, X::Dict)

Get Magnetic Local Time (MLT) from a Cartesian GEO position and date.

# Arguments
- `model::MagneticField`: The magnetic field model
- `X::Dict`: Dictionary with keys:
  - `dateTime` or `Time`: Date and time (DateTime or String)
  - `x1`, `x2`, `x3`: Position coordinates in GEO system
"""
function get_mlt(X::Dict)
    _, iyear, idoy, ut, _, _, _ = process_coords_time(X)
    xgeo = Float64[X["x1"], X["x2"], X["x3"]]

    # Initialize output
    mlt = Ref{Float64}(0.0)
    get_mlt1!(iyear, idoy, ut, xgeo, mlt)
    return mlt[]
end
