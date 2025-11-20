@testitem "Coords transform (CoordinateVector)" begin
    using Dates, IRBEM
    using Chairmarks

    time = DateTime(1996, 8, 28, 16, 46)
    pos = GEO(6.90274, -1.63624, 1.91669)
    @test GEO(time, pos) === pos
    @test_nowarn GSM(time, pos)
    @info @b GSM($time, $pos)
    @info GSM(time, pos)
end

@testitem "Coords transform time dimension validation" begin
    using Dates, IRBEM
    using IRBEM.StaticArrays

    pos = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]

    # Single time with single position should work
    @test_nowarn IRBEM.transform(DateTime(2024, 1, 1), [1.0, 4.0, 7.0], "GEO", "GSM")
    @test_nowarn IRBEM.transform(DateTime(2024, 1, 1), SA[1.0, 4.0, 7.0], "GEO", "GSM")

    # Single time with multiple positions should fail
    @test_throws AssertionError IRBEM.transform(DateTime(2024, 1, 1), pos, "GEO", "GSM")

    # Matching time array should work
    time_match = [DateTime(2024, 1, 1, h, 0, 0) for h in 1:3]
    @test_nowarn IRBEM.transform(time_match, pos, "GEO", "GSM")

    # Mismatched time array should fail
    time_mismatch = [DateTime(2024, 1, 1, h, 0, 0) for h in 1:2]
    @test_throws AssertionError IRBEM.transform(time_mismatch, pos, "GEO", "GSM")
end

@testitem "Coords transform with array views" begin
    using Dates, IRBEM

    # Test that views work correctly
    x = rand(5, 4)
    pos_full = x[1:3, :]
    pos_view = @view x[1:3, 1:2]
    time = [DateTime(2024, 1, 1, 12, 0, 0) for _ in 1:2]

    result_array = IRBEM.transform(time, pos_full[:, 1:2], "GEO", "GSM")
    result_view = IRBEM.transform(time, pos_view, "GEO", "GSM")

    @test result_array ≈ result_view
end
