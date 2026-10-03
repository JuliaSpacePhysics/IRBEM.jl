"Returns the version number of the IGRF model."
get_igrf_version() = (v = Ref{Int32}(0); get_igrf_version!(v); v[])

"Provides the repository version number of the fortran source code."
irbem_fortran_version() = (v = Ref{Int32}(0); irbem_fortran_version1!(v); v[])

"Provides the repository release tag of the fortran source code."
function irbem_fortran_release()
    v = zeros(UInt8, 80)
    irbem_fortran_release1!(v)
    return strip(String(v), [' ', '\0'])
end
