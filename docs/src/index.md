# IRBEM.jl

Julia wrapper for the IRBEM (International Radiation Belt Environment Modeling)

```@docs
IRBEM
```

# API Reference

## Computing magnetic field coordinates

```@docs
make_lstar
get_mlt
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

## Python interface

```@docs
IRBEM.PythonAPI
```

## Library information

```@autodocs
Modules = [IRBEM]
```
