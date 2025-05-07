"Returns the version number of the IGRF model."
get_igrf_version() = (v = Ref{Int32}(0); get_igrf_version!(v); v[])

"Provides the repository version number of the fortran source code."
irbem_fortran_version() = (v = Ref{Int32}(0); irbem_fortran_version1!(v); v[])

"Provides the repository release tag of the fortran source code."
irbem_fortran_release() = (
    v = Vector{UInt8}(undef, 80);
    irbem_fortran_release1!(v);
    strip(String(v))
)

