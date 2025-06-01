
"""
    PythonAPI

Python-like interface to IRBEM.jl

References: [IRBEM Python API](https://github.com/PRBEM/IRBEM/blob/main/python/IRBEM/IRBEM.py)
"""
module PythonAPI

using ..IRBEM
using ..IRBEM: get_datetime
import ..IRBEM: get_mlt

"""
    get_mlt(X::Dict)

X is a dictionary specifying the time and location in GEO coordinates.

See [IRBEM Python API](https://github.com/PRBEM/IRBEM/blob/e7cecb00caf97bb6357f063d2ba1aa76d71a3705/python/IRBEM/IRBEM.py#L584)
"""
get_mlt(X::Dict) = get_mlt(X["x1"], X["x2"], X["x3"], get_datetime(X))


end
