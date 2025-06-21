const param_indices = Dict(
    "Kp" => 1,
    "Dst" => 2,
    "dens" => 3,
    "velo" => 4,
    "Pdyn" => 5,
    "ByIMF" => 6,
    "BzIMF" => 7,
    "G1" => 8,
    "G2" => 9,
    "G3" => 10,
    "W1" => 11,
    "W2" => 12,
    "W3" => 13,
    "W4" => 14,
    "W5" => 15,
    "W6" => 16,
    "AL" => 17
)

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
const coord_sys_lookup_sym = NamedTuple(Symbol(k)=>v for (k,v) in coord_sys_lookup)
