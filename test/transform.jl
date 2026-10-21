@testitem "transform shapes" begin
    using Dates, DimensionalData
    using StaticArrays
    using IRBEM: transform

    t = DateTime(1996, 8, 28, 16, 46)
    ts = [t + Hour(i) for i in 0:3]
    A = [6.90274 1.0 -2.0 0.5; -1.63624 2.0 1.0 0.5; 1.91669 3.0 4.0 -6.0]
    B = transform(GEO => GSM, A, ts)
    cols = [transform(GEO => GSM, A[:, i], ts[i]) for i in 1:4]
    @test B ≈ reduce(hcat, cols)
    @test transform(GEO => GSM, SVector{3}(A[:, 1]), t) isa SVector{3, Float64}
    @test transform(GEO => GSM, A, t)[:, 3] ≈ transform(GEO => GSM, A[:, 3], t)
    @test transform(GEO => GSM, permutedims(A), ts; dim = 1) ≈ permutedims(B)
    @test transform(GEO => GSM, @view(A[:, 1:2]), ts[1:2]) ≈ B[:, 1:2]
    @test transform(GEO => GSM, SVector{3}.(eachcol(A)), ts) ≈ cols
    @test_throws DimensionMismatch transform(GEO => GSM, A, ts[1:2])

    # time-first DimArray: times from the container, dims and metadata kept
    da = DimArray(permutedims(A), (Ti(ts), Y(1:3)); metadata = Dict("COORDINATE_SYSTEM" => "GEO"))
    out = transform(GEO => GSM, da)
    @test out isa DimArray
    @test dims(out) == dims(da)
    @test metadata(out) == metadata(da)
    @test parent(out) ≈ permutedims(B)

    # names (case-insensitive) and frames are interchangeable; names may contain '2'
    @test transform("geo" => :GSM, A, ts) == B
    @test transform(:J2000 => "GEO", transform(GEO => J2000, A, ts), ts) ≈ A rtol = 1.0e-4
    @test transform(HEE => GEO, transform(GEO => HEE, A, ts), ts) ≈ A rtol = 1.0e-4
end

@testitem "transform representations" begin
    using Dates
    using GeoCotrans: Geodetic, Spherical, car2sph
    using IRBEM: transform

    t = DateTime(2015, 2, 2, 6, 12, 43)
    raw = [600.0, 60.0, 50.0] # IRBEM GDZ: alt, lat, lon
    gdz = GDZ(60.0, 50.0, 600.0)
    geo = transform("GDZ" => "GEO", raw, t)
    @test transform(GEO, gdz, t) isa CoordinateVector{GEO}
    @test transform(GEO, gdz, t) ≈ geo
    @test transform(GDZ => GEO, gdz, t) ≈ geo
    # typed representations of raw components would contradict the names' order
    @test_throws ArgumentError transform(GDZ => GEO, [60.0, 50.0, 600.0], t)
    @test_throws ArgumentError transform(GEO => (GEO(), Spherical()), raw, t)
    @test_throws ArgumentError make_lstar(t, [60.0, 50.0, 600.0], GDZ, (; Kp = 40.0); kext = "T89")
    @test transform(GDZ, GEO(geo...), t) ≈ gdz
    @test transform(GDZ, GEO(geo...), t) isa CoordinateVector{GEO, Geodetic}
    sph = CoordinateVector{GEO, Spherical}(car2sph(geo...)...)
    @test transform(GSM, sph, t) ≈ transform(GSM, GEO(geo...), t)
    @test transform((GEO(), Spherical()), GSM(1.0, 2.0, 3.0), t) ≈ car2sph(transform(GSM => GEO, [1.0, 2.0, 3.0], t)...)
    # IRBEM's SPH is (r, lat, lon) in degrees
    r, θ, φ = sph
    @test make_lstar(t, sph, (; Kp = 40.0); kext = "T89").Lm ≈ make_lstar(t, [r, 90 - rad2deg(θ), rad2deg(φ)], "SPH", (; Kp = 40.0); kext = "T89").Lm
    @test make_lstar(t, gdz, "GDZ", (; Kp = 40.0); kext = "T89").Lm == make_lstar(t, raw, "GDZ", (; Kp = 40.0); kext = "T89").Lm
    @test_throws ArgumentError make_lstar(t, gdz, "GEO")
    # the pre-0.3 GDZ(alt, lat, lon) order is caught
    @test_throws ArgumentError make_lstar(t, GDZ(600.0, 60.0, 50.0))
    @test_throws ArgumentError transform(GSE => GSM, GEO(geo...), t)
end

@testitem "agrees with GeoCotrans" begin
    using Dates, GeoCotrans

    t = DateTime(2010, 2, 25, 12)
    x = [3.0, -2.0, 1.5]
    @test IRBEM.transform(GEI => GEO, x, t) ≈ GeoCotrans.transform(GEI => GEO, x, t) rtol = 1.0e-5
    @test IRBEM.transform(GEO => GSM, x, t) ≈ GeoCotrans.transform(GEO => GSM, x, t) rtol = 1.0e-3
    @test IRBEM.transform(J2000 => GEI, x, t) ≈ GeoCotrans.transform(J2000 => GEI, x, t) rtol = 1.0e-5
end

@testitem "no clashes with GeoCotrans" begin
    import GeoCotrans
    m = Module()
    Core.eval(m, :(using GeoCotrans, IRBEM))
    for name in (:GEO, :GSM, :J2000, :GDZ, :CoordinateVector, :transform, :get_mlt, :find_magequator)
        @test getglobal(m, name) === getglobal(GeoCotrans, name)
    end
    @test getglobal(m, :HEE) === IRBEM.HEE
end
