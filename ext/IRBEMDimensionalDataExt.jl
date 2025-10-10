module IRBEMDimensionalDataExt

using DimensionalData: dimnum, dims, rebuild, layers
using DimensionalData: AbstractDimArray, AbstractDimStack, TimeDim
import IRBEM
import IRBEM: get_mlt, vecf

function IRBEM.get_mlt(A::AbstractDimArray)
    dim = dimnum(A, TimeDim)
    tdim = dims(A, dim)
    out = similar(A, tdim)
    odim = dim == 1 ? 2 : 1
    xyz = eachslice(parent(A), dims=odim)
    out .= get_mlt.(xyz[1], xyz[2], xyz[3], parent(tdim))
    return out
end


function IRBEM.get_mlt(A::AbstractDimStack)
    dim = dimnum(A, TimeDim)
    tdim = dims(A, dim)
    xyz = layers(A)
    @assert length(xyz) == 3
    out = similar(xyz[1], tdim)
    out .= get_mlt.(parent(xyz[1]), parent(xyz[2]), parent(xyz[3]), parent(tdim))
    return out
end

end
