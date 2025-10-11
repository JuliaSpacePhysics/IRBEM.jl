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
