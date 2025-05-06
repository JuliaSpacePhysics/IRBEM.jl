mutable struct MagneticField
    kext::Int32
    options::Vector{Int32}
    sysaxes::Int32
end

function MagneticField(; options=[0, 0, 0, 0, 0], kext="OPQ77", sysaxes="GDZ")
    kext_val = parse_kext(kext)
    sysaxes_val = coord_sys(sysaxes)
    MagneticField(kext_val, options, sysaxes_val)
end