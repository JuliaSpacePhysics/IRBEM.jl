"""
https://prbem.github.io/IRBEM/api/magnetic_coordinates.html#points-of-interest-on-the-field-line
"""

"""
    find_mirror_point(model::MagneticField, X::Dict, maginput::Dict, alpha)

Find the magnitude and location of the mirror point along a field line traced from any given location and local pitch-angle.

# Arguments
- `model::MagneticField`: The magnetic field model
- `X::Dict`: Dictionary with keys:
  - `dateTime` or `Time`: Date and time (DateTime or String)
  - `x1`, `x2`, `x3`: Position coordinates in the system specified by `sysaxes`
- `maginput::Dict`: Dictionary with magnetic field model inputs
- `alpha`: Local pitch angle in degrees

# Returns
- Blocal: magnitude of magnetic field at point (nT)
- Bmirr: magnitude of the magnetic field at the mirror point (nT)
- posit (array of 3 double): GEO coordinates of the mirror point (Re)

References: [IRBEM API](https://prbem.github.io/IRBEM/api/magnetic_coordinates.html#routine-FIND_MIRROR_POINT)
"""
function find_mirror_point(model::MagneticField, X::Dict, maginput::Dict, alpha)
    # Process input coordinates and time
    ntime, iyear, idoy, ut, x1, x2, x3 = process_coords_time(X)

    # Check if ntime exceeds 1 (find_mirror_point only supports single time)
    if ntime > 1
        throw(ArgumentError("find_mirror_point only supports a single time point"))
    end

    # Process magnetic field model inputs
    maginput_array = prepare_maginput(maginput, ntime)

    # Initialize output arrays
    Blocal = Ref{Float64}()
    Bmirr = Ref{Float64}()
    POSIT = zeros(Float64, 3)

    # Call IRBEM library function using @ccall
    kext = model.kext
    options = model.options
    sysaxes = model.sysaxes
    alpha = Float64(alpha)

    @ccall libirbem.find_mirror_point1_(
        kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        x1::Ptr{Float64}, x2::Ptr{Float64}, x3::Ptr{Float64},
        alpha::Ref{Float64}, maginput_array::Ptr{Float64},
        Blocal::Ref{Float64}, Bmirr::Ref{Float64}, posit::Ptr{Float64}
    )::Cvoid

    (; Blocal=Blocal[], Bmirr=Bmirr[], posit)
end

"""
    find_foot_point(model::MagneticField, X::Dict, maginput::Dict, stop_alt, hemi_flag)

Find the footprint of a field line that passes through location X in a given hemisphere.

# Arguments
- `model::MagneticField`: The magnetic field model
- `X::Dict`: Dictionary with keys:
  - `dateTime` or `Time`: Date and time (DateTime or String)
  - `x1`, `x2`, `x3`: Position coordinates in the system specified by `sysaxes`
- `maginput::Dict`: Dictionary with magnetic field model inputs
- `stop_alt`: Altitude in km where to stop field line tracing
- `hemi_flag`: Hemisphere flag (0: same as SM z, +1: northern, -1: southern)

# Outputs
- XFOOT: GDZ coordinates of the foot point (Re)
- BFOOT: magnetic field vector (GEO) at the foot point (nT)
- BFOOTMAG: magnitude of the magnetic field at the foot point (nT)
"""
function find_foot_point(model::MagneticField, X::Dict, maginput::Dict, stop_alt, hemi_flag)
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
    stop_alt = Float64(stop_alt)
    hemi_flag = Int32(hemi_flag)

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
    find_magequator(model::MagneticField, X::Dict, maginput::Dict)

Find the coordinates of the magnetic equator from tracing the magnetic field line from the input location.

# Arguments
- `model::MagneticField`: The magnetic field model
- `X::Dict`: Dictionary with keys:
  - `dateTime` or `Time`: Date and time (DateTime or String)
  - `x1`, `x2`, `x3`: Position coordinates in the system specified by `sysaxes`
- `maginput::Dict`: Dictionary with magnetic field model inputs

# Returns
- `Dict`: Contains keys Bmin and XGEO (location of magnetic equator in GEO coordinates)
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
    Bmin = Ref{Float64}(0.0)
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
        Bmin::Ref{Float64}, XGEO::Ptr{Float64}
    )::Cvoid

    # Return results as a dictionary
    return (; Bmin=Bmin[], XGEO)
end