"""
    make_lstar($SIG1)
    make_lstar($SIG2)

Compute magnetic coordinates at a spacecraft position.

# Arguments
$SIG_DOC

# Returns
- `NamedTuple`: Contains fields Lm, MLT, Blocal, Bmin, Lstar, and XJ
"""
function make_lstar(ntime::Int32, args...)
    Lm, Lstar, Blocal, Bmin, XJ, MLT = [zeros(Float64, ntime) for _ in 1:6]
    make_lstar1!(ntime, args..., Lm, Lstar, Blocal, Bmin, XJ, MLT)
    return map(_only, (; Lm, Lstar, Blocal, Bmin, XJ, MLT))
end

"""
    get_field_multi($SIG1)
    get_field_multi($SIG2)

Compute the GEO vector of the magnetic field at input location for a set of internal/external magnetic field.

# Arguments
$SIG_DOC

# Returns
- `NamedTuple`: Contains fields Bgeo (GEO components of B field) and Bmag (magnitude of B field)
"""
function get_field_multi(ntime::Int32, args...)
    # Initialize output arrays
    Bgeo = zeros(Float64, (3, ntime))
    Bmag = zeros(Float64, ntime)
    get_field_multi!(ntime, args..., Bgeo, Bmag)
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

for f in (:make_lstar, :get_field_multi)
    @eval $f(args...; kwargs...) = $f(prepare_irbem(args...; kwargs...)...)
end