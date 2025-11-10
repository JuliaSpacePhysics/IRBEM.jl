const _coord_sys_lookup = Dict(
    "GDZ" => 0,
    "GEO" => 1,
    "GSM" => 2,
    "GSE" => 3,
    "SM" => 4,
    "GEI" => 5,
    "MAG" => 6,
    "SPH" => 7,
    "RLL" => 8
)

const coord_sys_lookup = with_case_variants(_coord_sys_lookup)
# For symbols
const coord_sys_lookup_sym = NamedTuple(Symbol(k) => v for (k, v) in coord_sys_lookup)

const NTIME_MAX = let _NTIME_MAX = Ref{Int32}()
    get_irbem_ntime_max1!(_NTIME_MAX)
    _NTIME_MAX[]
end
