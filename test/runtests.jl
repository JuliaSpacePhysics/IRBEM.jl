using TestItems, TestItemRunner

@testsnippet Share begin
    using Dates
    model = MagneticField(options=[0, 0, 0, 0, 0], kext="T89")
    X = Dict(
        "dateTime" => DateTime("2015-02-02T06:12:43"),
        "x1" => 600.0,  # km
        "x2" => 60.0,   # lat
        "x3" => 50.0    # lon
    )
    maginput = Dict("Kp" => 40.0)

    n = 3
    X_array = Dict(
        # "dateTime" => fill(DateTime("2015-02-02T06:12:43"), n),
        "dateTime" => DateTime("2015-02-02T06:12:43"),
        "x1" => fill(600.0, n),  # km
        "x2" => fill(60.0, n),   # lat
        "x3" => fill(50.0, n)    # lon
    )
    maginput_array = Dict("Kp" => fill(40.0, n))
end

@testitem "make_lstar" setup = [Share] begin
    l_star_true_dict = Dict(
        "Lm" => 3.5597242229067536, "MLT" => 10.170297893176182,
        "blocal" => 42271.43059990003, "bmin" => 626.2258295723121,
        "Lstar" => -1e+31, "xj" => 7.020585390925573
    )
    result = make_lstar(model, X, maginput)
    @test result == l_star_true_dict
end

@testitem "find_foot_point" setup = [Share] begin
    _foot_point_true = (;
        XFOOT=[99.99412846343064, 61.113869939535036, 50.55633537632344],
        BFOOT=[-25644.012241653385, -25370.689449132995, -38649.994779664776],
        BFOOTMAG=52868.793663583165
    )
    stopAlt = 100
    hemiFlag = 0
    @test find_foot_point(model, X, maginput, stopAlt, hemiFlag) == _foot_point_true
end

@testitem "get_field_multi" setup = [Share] begin
    true_Bgeo = [
        -21079.764883133903 -21078.12221121894 -21078.12221121894;
        -21504.21460705096 -21508.430942943523 -21508.430942943523;
        -29666.24532305791 -29637.46273232981 -29637.46273232981
    ]
    true_Bl = [42271.43059990003, 42252.56246417121, 42252.56246417121]

    result = get_field_multi(model, X_array, maginput_array)
    @test result[1] ≈ true_Bgeo
    @test result[2] ≈ true_Bl
end

@testitem "find_mirror_point" setup = [Share] begin
    true_blocal = 42271.43059990003
    true_bmin = 42271.43059990003
    true_POSIT = [0.35282136776620165, 0.4204761325793738, 0.9448914452448274]

    alpha = 90.0  # Local pitch angle in degrees
    result = find_mirror_point(model, X, maginput, alpha)

    @test result[1] == true_blocal
    @test result[2] == true_bmin
    @test result[3] ≈ true_POSIT
end

@testitem "find_magequator" setup = [Share] begin
    true_bmin = 626.2258295723121
    true_XGEO = [2.1962220856733894, 2.8360222891612192, 0.3472455620354017]
    result = find_magequator(model, X, maginput)

    @test result[1] == true_bmin
    @test result[2] ≈ true_XGEO
end


@testitem "get_mlt" begin
    model = MagneticField(verbose=false, kext="T89")

    # Test with GEO coordinates
    X = Dict(
        "dateTime" => DateTime("2015-02-02T06:12:43"),
        "x1" => 2.195517156287977,
        "x2" => 2.834061428571752,
        "x3" => 0.34759070278576953
    )

    mlt = get_mlt(model, X)

    # MLT should be a scalar value between 0 and 24
    @test 0 <= mlt <= 24
end


@testitem "Utility functions" begin
    using Dates
    @test IRBEM.beta(100.0) ≈ 0.5482 atol = 1e-4
    @test IRBEM.gamma(100.0) ≈ 1.1956 atol = 1e-4

    # Test coordinate system lookup
    @test IRBEM.coord_sys("GDZ") == 0
    @test IRBEM.coord_sys("GEO") == 1
    @test IRBEM.coord_sys("GSM") == 2

    # Test datetime extraction
    dt = DateTime("2015-02-02T06:12:43")
    X = Dict("dateTime" => dt)
    @test IRBEM.get_datetime(X) == dt
    X_str = Dict("Time" => "2015-02-02T06:12:43")
    @test IRBEM.get_datetime(X_str) == dt
end