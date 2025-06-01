"""
    make_lstar($SIG1)
    make_lstar($SIG2)

Compute magnetic coordinates at a spacecraft position.

# Arguments
$SIG_DOC

# Returns
- `NamedTuple`: Contains fields Lm, MLT, Blocal, Bmin, Lstar, and XJ

# Examples
```jldoctest
julia> make_lstar("2015-02-02T06:12:43", [600.0, 60.0, 50.0], "GDZ", Dict("Kp" => 40.0))
(Lm = 3.5597242229067536, Lstar = -1.0e31, Blocal = 42271.43059990003, Bmin = 626.2258295723121, XJ = 7.020585390925573, MLT = 10.170297893176182)
```
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

function get_bderivs(ntime::Int32, args...)
    # Initialize output arrays
    Bgeo = zeros(Float64, (3, ntime))
    Bmag = zeros(Float64, ntime)
    gradBmag = zeros(Float64, (3, ntime))
    diffB = zeros(Float64, (3, 3, ntime))
    get_bderivs!(ntime, args..., Bgeo, Bmag, gradBmag, diffB)
    return map(_only, (; Bgeo, Bmag, gradBmag, diffB))
end

"""
    get_bderivs($SIG1)
    get_bderivs($SIG2)

Compute the magnetic field and its 1st-order derivatives at each input location.

# Arguments
$SIG_DOC

# Returns
- `NamedTuple`: Contains fields `Bgeo` (GEO components of B field), `Bmag` (magnitude of B field), `gradBmag` (gradients of Bmag in GEO), and `diffB` (derivatives of the magnetic field vector).

# Examples
```jldoctest
julia> get_bderivs("2015-02-02T06:12:43", [600.0, 60.0, 50.0], 0.1, "GDZ", Dict("Kp" => 40.0)) |> pprint
(Bgeo = [-21079.764883133903, -21504.21460705096, -29666.24532305791],
 Bmag = 42271.43059990003,
 gradBmag = [-49644.37271032293, -46030.37495428827, -83024.03530787815],
 diffB =
     [-13530.079906431165 31460.805163291334 53890.73134176735; 30427.464243221693 -16715.08632269888 50326.93737340687; 62620.43884602288 59981.93936448166 44395.53254933224;;;])
```
"""
function get_bderivs(arg1, arg2, dX, args...; kw...)
    return get_bderivs(prepare_irbem(arg1, arg2, args...; kw...)..., Float64(dX))
end

"""
    get_mlt(𝐫, time)
    get_mlt(x, y, z, time)

Get Magnetic Local Time (MLT) from a Cartesian GEO position `𝐫` and `time`.
"""
function get_mlt(𝐫, time)
    iyear, idoy, ut = decompose_time(time)
    xgeo = convert(Array{Float64}, 𝐫)
    mlt = Ref{Float64}(0.0)
    get_mlt1!(iyear, idoy, ut, xgeo, mlt)
    return mlt[]
end

get_mlt(x, y, z, time) = get_mlt(SA[x, y, z], time)

for f in (:make_lstar, :get_field_multi)
    @eval $f(args...; kwargs...) = $f(prepare_irbem(args...; kwargs...)...)
end