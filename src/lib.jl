# Computing magnetic field coordinates
make_lstar1!(ntime, kext, options, sysaxes, iyear, idoy, ut, x1, x2, x3, maginput_array, Lm, Lstar, Blocal, Bmin, XJ, mlt) =
    @ccall libirbem.make_lstar1_(
        ntime::Ref{Int32}, kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        x1::Ptr{Float64}, x2::Ptr{Float64}, x3::Ptr{Float64},
        maginput_array::Ptr{Float64},
        Lm::Ptr{Float64}, Lstar::Ptr{Float64},
        Blocal::Ptr{Float64}, Bmin::Ptr{Float64},
        XJ::Ptr{Float64}, mlt::Ptr{Float64}
    )::Cvoid

get_mlt1!(iyear, idoy, ut, xgeo, mlt) =
    @ccall libirbem.get_mlt1_(
        iyear::Ref{Int32}, idoy::Ref{Int32}, ut::Ref{Float64}, xgeo::Ptr{Float64},
        mlt::Ref{Float64}
    )::Cvoid

# Magnetic field computation
get_field_multi!(ntime, kext, options, sysaxes, iyear, idoy, ut, x1, x2, x3, maginput_array, Bgeo, Bmag) =
    @ccall libirbem.get_field_multi_(
        ntime::Ref{Int32},
        kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        x1::Ptr{Float64}, x2::Ptr{Float64}, x3::Ptr{Float64},
        maginput_array::Ptr{Float64}, Bgeo::Ptr{Float64}, Bmag::Ptr{Float64}
    )::Cvoid

# Coordinate transformations
coord_trans_vec1!(ntime, sys_in, sys_out, iyear, idoy, ut, pos_in, pos_out) =
    @ccall libirbem.coord_trans_vec1_(
        ntime::Ref{Int32}, sys_in::Ref{Int32}, sys_out::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        pos_in::Ptr{Float64}, pos_out::Ptr{Float64}
    )::Cvoid
