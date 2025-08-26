module IRBEMDimensionalDataExt

using DimensionalData: dimnum, dims, rebuild
using DimensionalData: AbstractDimArray, TimeDim
import IRBEM
import IRBEM: get_mlt, vecf

function IRBEM.get_mlt(A::AbstractDimArray)
    dim = dimnum(A, TimeDim)
    tdim = dims(A, dim)
    out = similar(A, tdim)
    out .= get_mlt.(eachslice(parent(A), dims=dim), parent(tdim))
    return out
end

end
