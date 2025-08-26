@testitem "DimensionalData Extension" begin
    using DimensionalData
    using Dates

    time = DateTime("2015-02-02T06:12:43")
    # Test with 1D DimArray (single position)
    times = [time]
    data = [2.195517156287977 2.834061428571752 0.34759070278576953]
    pos = DimArray(data, (Ti(times), Y(1:3)))
    mlt_result = get_mlt(pos)
    true_MLT = 9.56999052595853
    @test mlt_result isa DimArray
    @test mlt_result == [true_MLT]

    times = [time, time]
    pos = DimArray([data; data], (Ti(times), Y(1:3)))
    mlt_result = get_mlt(pos)
    true_MLT = 9.56999052595853
    @test mlt_result isa DimArray
    @test mlt_result == [true_MLT, true_MLT]
end
