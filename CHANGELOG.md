# Changelog

## [Unreleased]

### Changed

- **Breaking**: `coord` is a keyword: `make_lstar(time, x, maginput; coord = "GEO", kext)` instead of `make_lstar(time, x, "GEO", maginput; kext)`, likewise for the other `time, x` routines.
- **Breaking**: `coord` no longer accepts integer codes or coordinate types; pass the IRBEM name.
- **Breaking**: `drift_bounce_orbit(time, x, alpha, maginput; R0)`: `alpha` is positional and required, like in `find_mirror_point`.
- **Breaking**: `get_mlt(time, 𝐫)` only, like the other routines; `get_mlt(𝐫, time)` and `get_mlt(x, y, z, time)` are removed.
- **Breaking**: defaults are constants `kext = OPQ77`, `options = (0, 0, 0, 0, 0)`; the `IRBEM.KEXT` and `IRBEM.OPTIONS` refs are removed.
- **Breaking**: `transform(time, pos, in, out)` is `coord_trans(in, out, time, pos)` with IRBEM names only, and returns the shape of `pos`; the `in => out` and `"in2out"` forms are removed.

### Removed

- **Breaking**: `MagneticField` and the `f(model, X::AbstractDict, maginput)` signatures, incl. `get_mlt(X::AbstractDict)`; use `f(time, x, maginput; coord, kext, options)`.
- **Breaking**: IRBEM's own `CoordinateVector` and coordinate types (`GDZ`, `GEO`, `GSM`, ...), incl. `S(time, pos)` conversions. Pass raw components with `coord`, or GeoCotrans.jl `CoordinateVector`s (note GeoCotrans `GDZ(lat, lon, alt)`).

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

[unreleased]: https://github.com/JuliaSpacePhysics/IRBEM.jl/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/JuliaSpacePhysics/IRBEM.jl/compare/v0.1.7...v0.2.0
