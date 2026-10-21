# Changelog

## [Unreleased]

### Changed

- **Breaking**: coordinate types are GeoCotrans.jl's. `GEO`, `GSM`, `GSE`, `SM`, `GEI`, `MAG`, `J2000`, `GDZ` and `CoordinateVector` are re-exported from GeoCotrans; `HEE`, `HAE`, `HEEQ`, `TEME` subtype its frame type. `CoordinateVector` has no `sym` field.
- **Breaking**: `GDZ(lat, lon, alt)` instead of `GDZ(alt, lat, lon)`. A latitude outside [-90, 90] throws, which catches most old-order calls. IRBEM names keep IRBEM's order: `make_lstar(t, [alt, lat, lon], "GDZ")` is unchanged.
- **Breaking**: `transform(from => to, x, t; dim)`, `transform(from => to, A; dim)` and `transform(to, x::CoordinateVector, t)`, as in GeoCotrans, replace `transform(time, pos, in, out)`, `transform(time, pos, in => out)` and `transform(time, pos, "in2out")`.
- **Breaking**: `transform`, `get_mlt` and `find_magequator` are no longer exported (GeoCotrans exports its own); call them as `IRBEM.transform` etc.
- `get_mlt` converts a `CoordinateVector` in any system to GEO; it previously read its components as Cartesian GEO.

### Removed

- **Breaking**: `SPH`, `RLL`, `TOD` constructors and `S(time, pos::CoordinateVector)` conversions (e.g. `GSM(t, pos)`); use `IRBEM.transform(GSM, pos, t)`, `CoordinateVector{GEO, Spherical}`, or the names `"SPH"`, `"RLL"`, `"TOD"`.

## [0.2.0] - 2026-10-02

### Changed

- **Breaking**: unset magnetic field inputs are `-9999` (missing) instead of `0`, including `MagInput()` defaults. Pass every input your `kext` model uses; a model missing one now returns `NaN` rather than a result computed with zeros.
- **Breaking**: maginput keys must be `MagInput` field names (case-sensitive); unknown keys throw `ArgumentError`.
- **Breaking**: outputs IRBEM could not compute are `NaN`; check `isnan` instead of comparing with `-1e31`.
- **Breaking**: output shape follows the type of `time`: a scalar gives scalars, a vector gives one entry per point, even for a single point (`make_lstar([t], [x])` now returns vectors).
- **Breaking**: `find_mirror_point`, `find_foot_point` and `find_magequator` return one entry per point for vector inputs; previously only the first point was computed.
- **Breaking**: `trace_field_line` returns `Blocal` of length `Nposit`; `drift_shell` / `drift_bounce_orbit` fill `Blocal` past `Nposit` with `NaN`.
- **Breaking**: inputs that were silently misread now throw:
  - `DimensionMismatch` when time and position counts differ (`transform` previously threw `AssertionError`);
  - `ArgumentError` for several points passed to `trace_field_line`, `drift_shell` or `drift_bounce_orbit` (broadcast over points instead);
  - `ArgumentError` for `options` without 5 elements, and for unknown coordinate systems or conversion strings (previously `ErrorException`).
- Arithmetic on `CoordinateVector`s returns an `SVector`.

### Removed

- **Breaking**: `IRBEM.PythonAPI`; `get_mlt(::AbstractDict)` is now a method of `IRBEM.get_mlt`.
- **Breaking**: methods taking a leading `ntime::Int32` (`make_lstar`, `get_field_multi`, `get_bderivs`); use `(time, x, ...)` or `(model, X, ...)`.

### Added

- `landi2lstar`, `get_hemi`.
- Coordinate systems `TOD` and `TEME`; `HEE`, `HAE`, `HEEQ` and `J2000` now work in `transform` (e.g. `"geo2j2000"`).
- Calls are safe from multiple threads (serialized).

[unreleased]: https://github.com/JuliaSpacePhysics/GeoCotrans.jl/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/JuliaSpacePhysics/IRBEM.jl/compare/v0.1.7...HEAD
