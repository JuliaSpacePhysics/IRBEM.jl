@testitem "Coords transform (CoordinateVector)" begin
    using Dates, IRBEM

    time = DateTime(1996, 8, 28, 16, 46)
    pos = GEO(6.90274, -1.63624, 1.91669)
    @test GEO(time, pos) == pos
    @test GSM(time, pos) ≈ transform(time, collect(pos), "GEO", "GSM")
    @test GEO(time, GDZ(600, 60, 50)) ≈ GEO(time, GDZ(600.0, 60.0, 50.0))
    # Systems beyond RLL, incl. names containing '2'
    @test GEO(time, HEE(time, pos)) ≈ pos rtol = 1.0e-4
    @test transform(time, transform(time, collect(pos), "geo2j2000"), :J2000 => :GEO) ≈ pos rtol = 1.0e-4
end

@testitem "Coords transform time dimension validation" begin
    using Dates, IRBEM
    using IRBEM.StaticArrays

    pos = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
    t = DateTime(2024, 1, 1)

    @test transform(t, [1.0, 4.0, 7.0], "GEO", "GSM") ≈ transform(t, SA[1.0, 4.0, 7.0], "GEO", "GSM")
    @test_throws DimensionMismatch transform(t, pos, "GEO", "GSM")
    @test_throws DimensionMismatch transform([t, t], pos, "GEO", "GSM")
    @test transform([t, t, t], pos, "GEO", "GSM")[:, 2] ≈ transform(t, pos[:, 2], "GEO", "GSM")
end

@testitem "Coords transform with array views" begin
    using Dates, IRBEM

    x = rand(5, 4)
    pos_view = @view x[1:3, 1:2]
    time = [DateTime(2024, 1, 1, 12, 0, 0) for _ in 1:2]

    @test transform(time, x[1:3, 1:2], "GEO", "GSM") ≈ transform(time, pos_view, "GEO", "GSM")
end
