@testitem "coord_trans" begin
    using Dates

    t = DateTime(1996, 8, 28, 16, 46)
    ts = [t + Hour(i) for i in 0:3]
    A = [6.90274 1.0 -2.0 0.5; -1.63624 2.0 1.0 0.5; 1.91669 3.0 4.0 -6.0]
    B = coord_trans("GEO", "GSM", ts, A)
    @test B[:, 3] ≈ coord_trans("geo", "gsm", ts[3], A[:, 3])
    @test_throws DimensionMismatch coord_trans("GEO", "GSM", t, A)
    @test_throws DimensionMismatch coord_trans("GEO", "GSM", ts, permutedims(A))
    @test coord_trans("J2000", "GEO", ts, coord_trans("GEO", "J2000", ts, A)) ≈ A rtol = 1.0e-4
end

@testitem "CoordinateVector positions" begin
    using Dates
    using GeoCotrans: GDZ, GEO, CoordinateVector, Spherical, car2sph

    t = DateTime(2015, 2, 2, 6, 12, 43)
    mi = (; Kp = 40.0)
    raw = [600.0, 60.0, 50.0] # IRBEM's GDZ order: alt, lat, lon
    gdz = GDZ(60.0, 50.0, 600.0)
    geo = GEO(coord_trans("GDZ", "GEO", t, raw)...)
    Lm = make_lstar(t, raw, mi; kext = "T89").Lm
    @test make_lstar(t, gdz, mi; kext = "T89").Lm == Lm
    @test make_lstar(t, geo, mi; kext = "T89").Lm ≈ Lm rtol = 1.0e-6
    @test make_lstar(t, CoordinateVector{GEO, Spherical}(car2sph(geo...)...), mi; kext = "T89").Lm ≈ Lm rtol = 1.0e-6
    @test make_lstar([t, t], [geo, geo], mi; kext = "T89").Lm ≈ [Lm, Lm] rtol = 1.0e-6
    @test_throws ArgumentError make_lstar([t, t], [geo, gdz], mi; kext = "T89")
    @test_throws ArgumentError make_lstar(t, gdz, mi; coord = "GDZ", kext = "T89")
    @test get_mlt(gdz, t) ≈ get_mlt(geo, t)
end
