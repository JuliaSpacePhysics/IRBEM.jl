# IRBEM.jl

[![DOI](https://zenodo.org/badge/978838161.svg)](https://doi.org/10.5281/zenodo.17577161)

Julia wrapper for the [IRBEM](https://prbem.github.io/IRBEM/) (International Radiation Belt Environment Modeling) Fortran library.

```@docs
IRBEM
```

# API Reference

## External Magnetic Field Models and Model Inputs

See [IRBEM Documentation](https://prbem.github.io/IRBEM/api/general_information.html#external-magnetic-field-model) for more details. Most routines accept a `kext` keyword argument (of integer-like type) which allows the selection of the external magnetic field model.

```@docs
IRBEM.ExternalFieldModel
MagInput
MF75
TS87
TL87
T89
OPQ77
OPD88
T96
OM97
T01
T01S
T04
A00
T07
MT
```

## Coordinate systems

See [IRBEM Documentation](https://prbem.github.io/IRBEM/api/general_information.html#coordinate-systems) for more details.

```@docs
GDZ
GEO
GSM
GSE
SM
GEI
MAG
SPH
RLL
HEE
HAE
HEEQ
TOD
J2000
TEME
```

## Computing magnetic field coordinates

```@docs
make_lstar
landi2lstar
get_mlt
get_hemi
```


## Points of interest on the field line

```@docs
find_mirror_point
find_magequator
find_foot_point
```


## Magnetic field computation

```@docs
get_field_multi
get_bderivs
```


## Field tracing

```@docs
trace_field_line
drift_shell
drift_bounce_orbit
```


## Coordinates transformations

```@docs
transform
```

## Library information

```@docs
get_igrf_version
irbem_fortran_version
irbem_fortran_release
```
