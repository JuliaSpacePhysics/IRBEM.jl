"""
    make_lstar($SIG1)
    make_lstar($SIG2)

Compute magnetic coordinates at a spacecraft position.
Returns a named tuple with fields Lm, Lstar, Blocal, Bmin, XJ, and MLT.
`Lstar` is only computed when `options[1] > 0` (`NaN` otherwise).

# Arguments
$SIG_DOC

# Examples
```jldoctest
julia> make_lstar("2015-02-02T06:12:43", GDZ(600.0, 60.0, 50.0), (; Kp = 40.0); kext = "T89")
(Lm = 3.5597242229067536, Lstar = NaN, Blocal = 42271.43059990003, Bmin = 626.2258295723121, XJ = 7.020585390925573, MLT = 10.170297893176182)
```

Reference: [IRBEM API](https://prbem.github.io/IRBEM/api/magnetic_coordinates.html#routine-MAKE_LSTAR)
"""
make_lstar(args...; kw...) = _lstar(make_lstar1_!, prepare_irbem(args...; kw...))

"""
    landi2lstar($SIG1)
    landi2lstar($SIG2)

Like [`make_lstar`](@ref), but L* is deduced empirically from Lm, I and day of year, which is much faster.
Errors in L* are below 2%. Valid for locally mirroring particles only.
IRBEM forces IGRF + Olson-Pfitzer quiet ([`OPQ77`](@ref)), whatever `kext` and `options[5]` are.

# Arguments
$SIG_DOC

Reference: [IRBEM API](https://prbem.github.io/IRBEM/api/magnetic_coordinates.html#routine-LANDI2LSTAR)
"""
landi2lstar(args...; kw...) = _lstar(landi2lstar1_!, prepare_irbem(args...; kw...))

function _lstar(f!, (p, single))
    nt = NamedTuple{(:Lm, :Lstar, :Blocal, :Bmin, :XJ, :MLT)}(ntuple(_ -> _out(p), 6))
    f!(p..., nt...)
    return _output(nt, single)
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
function get_field_multi(args...; kw...)
    p, single = prepare_irbem(args...; kw...)
    nt = (; Bgeo = _out(p, 3), Bmag = _out(p))
    get_field_multi!(p..., nt...)
    return _output(nt, single)
end

"""
    get_bderivs($SIG1)
    get_bderivs($SIG2)

Compute the magnetic field and its 1st-order derivatives at each input location.

# Arguments
$SIG_DOC
- `dX`: step size (Re) for the finite differences, placed right after the position

# Returns
- `NamedTuple`: Contains fields `Bgeo` (GEO components of B field), `Bmag` (magnitude of B field), `gradBmag` (gradients of Bmag in GEO), and `diffB` (derivatives of the magnetic field vector).

# Examples
```jldoctest
julia> get_bderivs("2015-02-02T06:12:43", GDZ(600.0, 60.0, 50.0), 0.1, (; Kp = 40.0); kext = "T89") |> pprint
(Bgeo = [-21079.764883133903, -21504.21460705096, -29666.24532305791],
 Bmag = 42271.43059990003,
 gradBmag = [-49644.37271032293, -46030.37495428827, -83024.03530787815],
 diffB =
     [-13530.079906431165 31460.805163291334 53890.73134176735; 30427.464243221693 -16715.08632269888 50326.93737340687; 62620.43884602288 59981.93936448166 44395.53254933224])
```
"""
function get_bderivs(arg1, arg2, dX, args...; kw...)
    p, single = prepare_irbem(arg1, arg2, args...; kw...)
    nt = (; Bgeo = _out(p, 3), Bmag = _out(p), gradBmag = _out(p, 3), diffB = _out(p, 3, 3))
    get_bderivs!(p..., Float64(dX), nt...)
    return _output(nt, single)
end

"""
    get_hemi($SIG1)
    get_hemi($SIG2)

Magnetic hemisphere of the input location: `+1` northern, `-1` southern, `0` invalid magnetic field.

# Arguments
$SIG_DOC

Reference: [IRBEM API](https://prbem.github.io/IRBEM/api/magnetic_coordinates.html#routine-GET_HEMI_MULTI)
"""
function get_hemi(args...; kw...)
    p, single = prepare_irbem(args...; kw...)
    xhemi = Vector{Int32}(undef, p.ntime)
    get_hemi_multi!(p..., xhemi)
    return _output((; xhemi), single).xhemi
end

"""
    get_mlt(𝐫, time)
    get_mlt(time, 𝐫::AbstractVector)
    get_mlt(x, y, z, time)
    get_mlt(X::AbstractDict)

Get Magnetic Local Time (MLT) from a Cartesian GEO position `𝐫` and `time`.
`X` is a dictionary with keys `x1`, `x2`, `x3` (GEO position) and `dateTime` or `Time`.
"""
function get_mlt(𝐫, time)
    iyear, idoy, ut = decompose_time_s(time)
    mlt = Ref{Float64}()
    get_mlt1!(iyear, idoy, ut, SVector{3, Float64}(𝐫), mlt)
    return mlt[]
end

get_mlt(x, y, z, time) = get_mlt(SA_F64[x, y, z], time)
get_mlt(time, 𝐫::AbstractVector) = get_mlt(𝐫, time)
get_mlt(X::AbstractDict) = get_mlt(X["x1"], X["x2"], X["x3"], get_datetime(X))
