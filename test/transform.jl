@testitem "Coords transform (CoordinateVector)" begin
    using Dates, IRBEM

    time = DateTime(1996, 8, 28, 16, 46)
    pos = GEO(6.90274, -1.63624, 1.91669)
    @test GEO(time, pos) === pos
    @test_nowarn GSM(time, pos)
end