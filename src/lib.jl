# Raw wrappers: common inputs first, then routine-specific inputs, then outputs.
# The `@ccall` lists arguments in Fortran order.
# IRBEM keeps state in COMMON blocks and SAVEd variables, so all calls are serialized.
const LIBIRBEM_LOCK = ReentrantLock()

# Computing magnetic field coordinates
for f in (:make_lstar1_, :landi2lstar1_)
    @eval $(Symbol(f, :!))(ntime, kext, options, sysaxes, iyear, idoy, ut, x1, x2, x3, maginput, Lm, Lstar, Blocal, Bmin, XJ, mlt) =
        @lock LIBIRBEM_LOCK @ccall libirbem.$f(
            ntime::Ref{Int32}, kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
            iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
            x1::Ptr{Float64}, x2::Ptr{Float64}, x3::Ptr{Float64},
            maginput::Ptr{Float64},
            Lm::Ptr{Float64}, Lstar::Ptr{Float64},
            Blocal::Ptr{Float64}, Bmin::Ptr{Float64},
            XJ::Ptr{Float64}, mlt::Ptr{Float64}
        )::Cvoid
end

@inline get_mlt1!(iyear, idoy, ut, xgeo, mlt) =
    @lock LIBIRBEM_LOCK @ccall libirbem.get_mlt1_(
        iyear::Ref{Int32}, idoy::Ref{Int32}, ut::Ref{Float64}, xgeo::Ptr{Float64},
        mlt::Ref{Float64}
    )::Cvoid

get_hemi_multi!(ntime, kext, options, sysaxes, iyear, idoy, ut, x1, x2, x3, maginput, xhemi) =
    @lock LIBIRBEM_LOCK @ccall libirbem.get_hemi_multi_(
        ntime::Ref{Int32}, kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        x1::Ptr{Float64}, x2::Ptr{Float64}, x3::Ptr{Float64},
        maginput::Ptr{Float64}, xhemi::Ptr{Int32}
    )::Cvoid

# Points of interest on the field line
find_mirror_point1!(kext, options, sysaxes, iyear, idoy, ut, x1, x2, x3, maginput, alpha, Blocal, Bmirr, posit) =
    @lock LIBIRBEM_LOCK @ccall libirbem.find_mirror_point1_(
        kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        x1::Ptr{Float64}, x2::Ptr{Float64}, x3::Ptr{Float64},
        alpha::Ref{Float64}, maginput::Ptr{Float64},
        Blocal::Ptr{Float64}, Bmirr::Ptr{Float64}, posit::Ptr{Float64}
    )::Cvoid

find_foot_point1!(kext, options, sysaxes, iyear, idoy, ut, x1, x2, x3, maginput, stop_alt, hemi_flag, XFOOT, BFOOT, BFOOTMAG) =
    @lock LIBIRBEM_LOCK @ccall libirbem.find_foot_point1_(
        kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        x1::Ptr{Float64}, x2::Ptr{Float64}, x3::Ptr{Float64},
        stop_alt::Ref{Float64}, hemi_flag::Ref{Int32}, maginput::Ptr{Float64},
        XFOOT::Ptr{Float64}, BFOOT::Ptr{Float64}, BFOOTMAG::Ptr{Float64}
    )::Cvoid

find_magequator1!(kext, options, sysaxes, iyear, idoy, ut, x1, x2, x3, maginput, Bmin, XGEO) =
    @lock LIBIRBEM_LOCK @ccall libirbem.find_magequator1_(
        kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        x1::Ptr{Float64}, x2::Ptr{Float64}, x3::Ptr{Float64},
        maginput::Ptr{Float64},
        Bmin::Ptr{Float64}, XGEO::Ptr{Float64}
    )::Cvoid

# Magnetic field computation
get_field_multi!(ntime, kext, options, sysaxes, iyear, idoy, ut, x1, x2, x3, maginput, Bgeo, Bmag) =
    @lock LIBIRBEM_LOCK @ccall libirbem.get_field_multi_(
        ntime::Ref{Int32},
        kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        x1::Ptr{Float64}, x2::Ptr{Float64}, x3::Ptr{Float64},
        maginput::Ptr{Float64}, Bgeo::Ptr{Float64}, Bmag::Ptr{Float64}
    )::Cvoid

get_bderivs!(ntime, kext, options, sysaxes, iyear, idoy, ut, x1, x2, x3, maginput, dX, Bgeo, Bmag, gradBmag, diffB) =
    @lock LIBIRBEM_LOCK @ccall libirbem.get_bderivs_(
        ntime::Ref{Int32},
        kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32}, dX::Ref{Float64},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        x1::Ptr{Float64}, x2::Ptr{Float64}, x3::Ptr{Float64},
        maginput::Ptr{Float64}, Bgeo::Ptr{Float64}, Bmag::Ptr{Float64},
        gradBmag::Ptr{Float64}, diffB::Ptr{Float64}
    )::Cvoid

# Coordinate transformations
coord_trans_vec1!(ntime, sys_in, sys_out, iyear, idoy, ut, pos_in, pos_out) =
    @lock LIBIRBEM_LOCK @ccall libirbem.coord_trans_vec1_(
        ntime::Ref{Int32}, sys_in::Ref{Int32}, sys_out::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        pos_in::Ptr{Float64}, pos_out::Ptr{Float64}
    )::Cvoid

# [Field tracing](https://prbem.github.io/IRBEM/api/magnetic_coordinates.html#field-tracing)
trace_field_line2_1_!(kext, options, sysaxes, iyear, idoy, ut, x1, x2, x3, maginput, R0, Lm, Blocal, Bmin, XJ, posit, Nposit) =
    @lock LIBIRBEM_LOCK @ccall libirbem.trace_field_line2_1_(
        kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        x1::Ptr{Float64}, x2::Ptr{Float64}, x3::Ptr{Float64},
        maginput::Ptr{Float64}, R0::Ref{Float64},
        Lm::Ref{Float64}, Blocal::Ptr{Float64}, Bmin::Ref{Float64},
        XJ::Ref{Float64}, posit::Ptr{Float64}, Nposit::Ref{Int32}
    )::Cvoid

drift_bounce_orbit2_1_!(kext, options, sysaxes, iyear, idoy, ut, x1, x2, x3, maginput, alpha, R0, Lm, Lstar, Blocal, Bmin, Bmirr, XJ, posit, Nposit, hmin, hmin_lon) =
    @lock LIBIRBEM_LOCK @ccall libirbem.drift_bounce_orbit2_1_(
        kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        x1::Ptr{Float64}, x2::Ptr{Float64}, x3::Ptr{Float64},
        alpha::Ref{Float64}, maginput::Ptr{Float64},
        R0::Ref{Float64}, Lm::Ref{Float64}, Lstar::Ref{Float64},
        Blocal::Ptr{Float64}, Bmin::Ref{Float64}, Bmirr::Ref{Float64},
        XJ::Ref{Float64}, posit::Ptr{Float64}, Nposit::Ptr{Int32},
        hmin::Ref{Float64}, hmin_lon::Ref{Float64}
    )::Cvoid

drift_shell1_!(kext, options, sysaxes, iyear, idoy, ut, x1, x2, x3, maginput, Lm, Lstar, Blocal, Bmin, XJ, posit, Nposit) =
    @lock LIBIRBEM_LOCK @ccall libirbem.drift_shell1_(
        kext::Ref{Int32}, options::Ptr{Int32}, sysaxes::Ref{Int32},
        iyear::Ptr{Int32}, idoy::Ptr{Int32}, ut::Ptr{Float64},
        x1::Ptr{Float64}, x2::Ptr{Float64}, x3::Ptr{Float64},
        maginput::Ptr{Float64},
        Lm::Ref{Float64}, Lstar::Ref{Float64},
        Blocal::Ptr{Float64}, Bmin::Ref{Float64},
        XJ::Ref{Float64}, posit::Ptr{Float64}, Nposit::Ptr{Int32}
    )::Cvoid

# Library information functions
irbem_fortran_version1!(version) =
    @lock LIBIRBEM_LOCK @ccall libirbem.irbem_fortran_version1_(version::Ref{Int32})::Cvoid

irbem_fortran_release1!(version) =
    @lock LIBIRBEM_LOCK @ccall libirbem.irbem_fortran_release1_(version::Ptr{UInt8})::Cvoid

get_igrf_version!(version) =
    @lock LIBIRBEM_LOCK @ccall libirbem.get_igrf_version_(version::Ref{Int32})::Cvoid
