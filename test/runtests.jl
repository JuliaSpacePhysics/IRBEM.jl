# Reference: https://github.com/PRBEM/IRBEM/blob/main/python/IRBEM/test_IRBEM.py

using TestItems, TestItemRunner
@run_package_tests

@testitem "Aqua" begin
    using Aqua
    Aqua.test_all(IRBEM)
end

@testitem "JET" begin
    using JET
    IRBEM.workload()
    @test_call IRBEM.workload()
    test_package(IRBEM)
end

@testsnippet Share begin
    using Dates
    kext = "T89"
    model = MagneticField(options = [0, 0, 0, 0, 0], kext = kext)
    dipol_model = MagneticField(options = [0, 0, 5, 0, 5], kext = 0)
    t = DateTime("2015-02-02T06:12:43")
    x = [600.0, 60.0, 50.0]
    X = Dict(
        "dateTime" => t,
        "x1" => x[1],  # km
        "x2" => x[2],   # lat
        "x3" => x[3]    # lon
    )
    maginput = Dict(:Kp => 40.0)
    maginput_nt = (; Kp = 40.0)

    n = 3
    X_array = Dict(
        "dateTime" => fill(t, n),
        "x1" => fill(x[1], n),  # km
        "x2" => fill(x[2], n),   # lat
        "x3" => fill(x[3], n)    # lon
    )
    maginput_array = Dict("Kp" => fill(40.0, n))
    maginput_array_nt = (; Kp = fill(40.0, n))

    approx(a, b; kw...) = all(map((u, v) -> isapprox(u, v; nans = true, kw...), a, b))

    function _compute_dipole_L_shell(posit)
        x = posit[1, :, :]
        y = posit[2, :, :]
        z = posit[3, :, :]
        r = sqrt.(x .^ 2 .+ y .^ 2 .+ z .^ 2)
        theta = atan.(sqrt.(x .^ 2 .+ y .^ 2), z)
        return r ./ sin.(theta) .^ 2
    end
end

@testitem "make_lstar" setup = [Share] begin
    l_star_true = (
        Lm = 3.5597242229067536, Lstar = NaN,
        Blocal = 42271.43059990003, Bmin = 626.2258295723121,
        XJ = 7.020585390925573, MLT = 10.170297893176182,
    )
    @test approx(make_lstar(model, X, maginput), l_star_true)
    @test approx(make_lstar(t, x, "GDZ", maginput; kext), l_star_true)
    @test approx(make_lstar(t, x, :GDZ, maginput_nt; kext), l_star_true)
    @test approx(make_lstar(t, GDZ(x), maginput_nt; kext), l_star_true)
    @test approx(make_lstar([t, t], [x, x], "GDZ", (; Kp = [40.0, 40.0]); kext), map(v -> fill(v, 2), l_star_true))
    # A scalar maginput is shared by all points
    @test make_lstar([t, t, t], [x, x, x], "GDZ", maginput; kext).Lm ≈ fill(l_star_true.Lm, 3)
    # The coordinate system of a vector of CoordinateVector comes from its elements
    xgeo = GEO(transform(t, x, "GDZ", "GEO"))
    @test make_lstar([t, t], [xgeo, xgeo], maginput_nt; kext).Lm ≈ fill(l_star_true.Lm, 2) rtol = 1.0e-6
    @test_throws DimensionMismatch make_lstar([t, t], x, "GDZ", maginput_nt; kext)
    # Unset model inputs are missing for IRBEM, not zero
    @test isnan(make_lstar(t, GDZ(x); kext).Lm)
end

@testitem "landi2lstar" setup = [Share] begin
    pos = GEO(4.0, 0.0, 0.5)
    exact = make_lstar(t, pos; kext = OPQ77, options = [1, 0, 0, 0, 0])
    fast = landi2lstar(t, pos)
    @test fast.Lm == exact.Lm
    @test fast.Lstar ≈ exact.Lstar rtol = 0.02
    # IRBEM overwrites options in place; user options must not change
    m = MagneticField(options = [0, 0, 0, 0, 1], kext = 5, sysaxes = "GEO")
    landi2lstar(m, Dict("dateTime" => t, "x1" => pos[1], "x2" => pos[2], "x3" => pos[3]))
    @test m.options == [0, 0, 0, 0, 1]
end

@testitem "get_hemi" setup = [Share] begin
    @test get_hemi([t, t], [GDZ(x), GDZ(600.0, -60.0, 50.0)], maginput_nt; kext) == [1, -1]
    @test get_hemi(model, X, maginput) == 1
end

@testitem "get_field_multi" setup = [Share] begin
    true_Bgeo = [-21079.764883133903, -21504.21460705096, -29666.24532305791]
    true_Bl = 42271.43059990003

    result = get_field_multi(model, X_array, maginput_array)
    @test result.Bgeo ≈ repeat(true_Bgeo, 1, n)
    @test result.Bmag ≈ fill(true_Bl, n)
    @test approx(get_field_multi(model, X, maginput), (true_Bgeo, true_Bl))
    @test approx(result, get_field_multi(model, X_array, maginput_array_nt))
end

@testitem "get_bderivs" setup = [Share] begin
    res = get_bderivs(model, X, 0.1, maginput)
    @test res.Bmag ≈ 42271.43059990003
    @test res.gradBmag ≈ [-49644.37271032293, -46030.37495428827, -83024.03530787815]
    @test size(res.diffB) == (3, 3)
    @test approx(res, get_bderivs(t, GDZ(x), 0.1, maginput_nt; kext))
end

@testitem "get_mlt" setup = [Share] begin
    # Corresponds to test_get_mlt in Python
    input_dict = Dict(
        "dateTime" => t,
        "x1" => 2.195517156287977,
        "x2" => 2.834061428571752,
        "x3" => 0.34759070278576953
    )
    r = [2.195517156287977, 2.834061428571752, 0.34759070278576953]

    true_MLT = 9.56999052595853
    @test get_mlt(input_dict) ≈ true_MLT
    @test get_mlt(r, t) == get_mlt(t, r) ≈ true_MLT
end

@testitem "trace_field_line" setup = [Share] begin
    res = trace_field_line(model, X, maginput)
    @test length(res.Blocal) == size(res.posit, 2) == res.Nposit
    @test approx(res, trace_field_line(t, GDZ(x), maginput_nt; kext))
    @test_throws ArgumentError trace_field_line([t, t], [x, x], "GDZ", maginput_nt; kext)
end

@testitem "drift_shell" setup = [Share] begin
    using NaNStatistics
    res = drift_shell(dipol_model, X, maginput)
    Lm = res.Lm
    @test Lm ≈ 4.326679 atol = 1.0e-2
    L_posit = _compute_dipole_L_shell(res.posit)
    @test nanmaximum(abs.(L_posit .- Lm)) / Lm <= 1.0e-2
    @test count(!isnan, res.Blocal) == sum(res.Nposit)
end

@testitem "drift_bounce_orbit" setup = [Share] begin
    using NaNStatistics
    res = drift_bounce_orbit(dipol_model, X, maginput)
    Lm = res.Lm
    @test Lm ≈ 4.326679 atol = 1.0e-2
    alt = X["x1"]
    @test abs((res.hmin - alt) / alt) <= 1.0e-2
    L_posit = _compute_dipole_L_shell(res.posit)
    @test nanmaximum(abs.(L_posit .- Lm)) / Lm <= 1.0e-2
end

@testitem "Thread safety" begin
    code = """
    using IRBEM, Dates
    t = DateTime(2015, 2, 2, 6, 12, 43)
    f(i) = make_lstar(t, GDZ(600.0 + 10i, 60.0, 50.0), (; Kp = 40.0); kext = isodd(i) ? "T89" : "OPQ77").Lm
    serial = f.(1:100)
    threaded = similar(serial)
    Threads.@threads for i in 1:100
        threaded[i] = f(i)
    end
    exit(isequal(serial, threaded) ? 0 : 1)
    """
    @test success(`$(Base.julia_cmd()) -t 4 --project=$(Base.active_project()) -e $code`)
end

@testitem "Coords transform (single/multi entry)" begin
    using Dates, IRBEM

    # Single entry: GEO→GEO
    time = DateTime(1996, 8, 28, 16, 46)
    pos = [6.90274, -1.63624, 1.91669]
    @test transform(time, pos, "GEO", "GEO") == pos
    @test transform(time, pos, "GEO" => "GEO") == pos
    @test transform(time, pos, "geo2geo") == pos

    # Multi entry: GEO→MAG
    times = [DateTime(1996, 8, 28, 16, 46), DateTime(1996, 8, 28, 16, 46)]
    poses = [[6.90274, -1.63624, 1.91669] [6.90274, -1.63624, 1.91669]]
    multi_result = transform(times, poses, "GEO", "MAG")
    @test multi_result[:, 1] == multi_result[:, 2] == transform(time, pos, :GEO, :MAG)
    @test size(multi_result) == size(poses)
end

@testitem "Utility functions" begin
    using Dates

    @test IRBEM.parse_kext("None") == 0
    @test IRBEM.parse_kext("OPQ77") == IRBEM.parse_kext(:opq77) == 5
    @test IRBEM.parse_kext(5) == 5
    @test_throws ArgumentError IRBEM.parse_kext("T99")

    @test IRBEM.coord_sys("GDZ") == IRBEM.coord_sys(GDZ) == IRBEM.coord_sys(:gdz) == 0
    @test IRBEM.coord_sys("gsm") == 2
    @test IRBEM.coord_sys(TEME) == 14
    @test IRBEM.parse_coord_transform("J20002geo") == (:J2000, :GEO)
    @test IRBEM.parse_coord_transform("gsm_to_sm") == (:GSM, :SM)

    @test_throws ArgumentError make_lstar(DateTime(2020), GDZ(600.0, 60.0, 50.0), (; kp = 40.0))
    @test_throws ArgumentError make_lstar(DateTime(2020), GDZ(600.0, 60.0, 50.0); options = [0, 0, 0])

    dt = DateTime("2015-02-02T06:12:43")
    @test IRBEM.get_datetime(Dict("dateTime" => dt)) == dt
    @test IRBEM.get_datetime(Dict("Time" => "2015-02-02T06:12:43")) == dt
end
