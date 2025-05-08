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

# Points of interest on the field line
find_mirror_point1!(kext, options, sysaxes, iyear, idoy, ut, x1, x2, x3, maginput, alpha, Blocal, Bmirr, posit) =
    @ccall libirbem.find_mirror_point1_(
        kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        x1::Ptr{Float64}, x2::Ptr{Float64}, x3::Ptr{Float64},
        alpha::Ref{Float64}, maginput::Ptr{Float64},
        Blocal::Ref{Float64}, Bmirr::Ref{Float64}, posit::Ptr{Float64}
    )::Cvoid

find_foot_point1!(kext, options, sysaxes, iyear, idoy, ut, x1, x2, x3, maginput, stop_alt, hemi_flag, BFOOTMAG, XFOOT, BFOOT) =
    @ccall libirbem.find_foot_point1_(
        kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        x1::Ptr{Float64}, x2::Ptr{Float64}, x3::Ptr{Float64},
        stop_alt::Ref{Float64}, hemi_flag::Ref{Int32}, maginput::Ptr{Float64},
        XFOOT::Ptr{Float64}, BFOOT::Ptr{Float64}, BFOOTMAG::Ptr{Float64}
    )::Cvoid

find_magequator1!(kext, options, sysaxes, iyear, idoy, ut, x1, x2, x3, maginput, Bmin, XGEO) =
    @ccall libirbem.find_magequator1_(
        kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        x1::Ptr{Float64}, x2::Ptr{Float64}, x3::Ptr{Float64},
        maginput::Ptr{Float64},
        Bmin::Ref{Float64}, XGEO::Ptr{Float64}
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

# Library information functions
"""
Returns the size of time dimension in inputs and/or output arrays for some of the routines.

Reference: [IRBEM Documentation](https://prbem.github.io/IRBEM/api/library_infos.html#routine-GET_IRBEM_NTIME_MAX)
"""
get_irbem_ntime_max1!(NTIME_MAX) =
    @ccall libirbem.get_irbem_ntime_max1_(NTIME_MAX::Ref{Int32})::Cvoid

irbem_fortran_version1!(version) =
    @ccall libirbem.irbem_fortran_version1_(version::Ref{Int32})::Cvoid

irbem_fortran_release1!(version) =
    @ccall libirbem.irbem_fortran_release1_(version::Ptr{UInt8})::Cvoid

get_igrf_version!(version) =
    @ccall libirbem.get_igrf_version_(version::Ref{Int32})::Cvoid