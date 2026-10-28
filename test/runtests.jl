# Reference: https://github.com/PRBEM/IRBEM/blob/main/python/IRBEM/test_IRBEM.py

using TestItems, TestItemRunner
@run_package_tests filter = ti -> isempty(VERSION.prerelease) || !(:jet in ti.tags)

@testitem "Aqua" begin
    using Aqua
    Aqua.test_all(IRBEM)
end

@testitem "JET" tags = [:jet] begin
    using JET
    IRBEM.workload()
    @test_call IRBEM.workload()
    test_package(IRBEM; target_modules = (IRBEM,))
end

@testsnippet Share begin
    using Dates
    kext = "T89"
    dipole = (; kext = 0, options = [0, 0, 5, 0, 5])
    t = DateTime("2015-02-02T06:12:43")
    x = [600.0, 60.0, 50.0] # GDZ: alt, lat, lon
    maginput = Dict(:Kp => 40.0)
    maginput_nt = (; Kp = 40.0)

    n = 3
    ts, xs = fill(t, n), fill(x, n)
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
    @test approx(make_lstar(t, x, maginput; coord = "GDZ", kext), l_star_true)
    @test approx(make_lstar(t, x, maginput_nt; kext), l_star_true)
    @test approx(make_lstar([t, t], [x, x], (; Kp = [40.0, 40.0]); kext), map(v -> fill(v, 2), l_star_true))
    # A scalar maginput is shared by all points
    @test make_lstar([t, t, t], [x, x, x], maginput; kext).Lm ≈ fill(l_star_true.Lm, 3)
    @test_throws DimensionMismatch make_lstar([t, t], x, maginput_nt; kext)
    # Output shape follows the time input, not the number of points
    @test make_lstar([t], [x], maginput_nt; kext).Lm isa Vector
    @inferred make_lstar(t, x, maginput_nt; kext)
    # Unset model inputs are missing for IRBEM, not zero
    @test isnan(make_lstar(t, x; kext).Lm)
end

@testitem "landi2lstar" setup = [Share] begin
    pos = [4.0, 0.0, 0.5]
    exact = make_lstar(t, pos; coord = "GEO", kext = OPQ77, options = [1, 0, 0, 0, 0])
    fast = landi2lstar(t, pos; coord = "GEO")
    @test fast.Lm == exact.Lm
    @test fast.Lstar ≈ exact.Lstar rtol = 0.02
    # IRBEM overwrites options in place; user options must not change
    options = [0, 0, 0, 0, 1]
    landi2lstar(t, pos; coord = "GEO", kext = 5, options)
    @test options == [0, 0, 0, 0, 1]
end

@testitem "get_hemi" setup = [Share] begin
    @test get_hemi([t, t], [x, [600.0, -60.0, 50.0]], maginput_nt; kext) == [1, -1]
end

@testitem "get_field_multi" setup = [Share] begin
    true_Bgeo = [-21079.764883133903, -21504.21460705096, -29666.24532305791]
    true_Bl = 42271.43059990003

    result = get_field_multi(ts, xs, maginput_array; kext)
    @test result.Bgeo ≈ repeat(true_Bgeo, 1, n)
    @test result.Bmag ≈ fill(true_Bl, n)
    @test approx(get_field_multi(t, x, maginput; kext), (true_Bgeo, true_Bl))
    @test approx(result, get_field_multi(ts, xs, maginput_array_nt; kext))
end

@testitem "get_bderivs" setup = [Share] begin
    res = get_bderivs(t, x, 0.1, maginput_nt; kext)
    @test res.Bmag ≈ 42271.43059990003
    @test res.gradBmag ≈ [-49644.37271032293, -46030.37495428827, -83024.03530787815]
    @test size(res.diffB) == (3, 3)
end

@testitem "get_mlt" setup = [Share] begin
    # Corresponds to test_get_mlt in Python
    r = [2.195517156287977, 2.834061428571752, 0.34759070278576953]

    true_MLT = 9.56999052595853
    @test get_mlt(t, r) ≈ true_MLT

    using DimensionalData, TimeseriesUtilities
    t2 = [t, t + Hour(6)]
    A = [r r]
    mlt2 = get_mlt.(t2, eachcol(A))
    @test get_mlt(Timed(A, t2; dim = 2)) == Timed(mlt2, t2)
    pos = DimArray(permutedims(A), (Ti(t2), Y(1:3)))
    @test get_mlt(pos) == mlt2
end

@testitem "trace_field_line" setup = [Share] begin
    res = trace_field_line(t, x, maginput_nt; kext)
    @test length(res.Blocal) == size(res.posit, 2) == res.Nposit
    @test_throws ArgumentError trace_field_line([t, t], [x, x], maginput_nt; kext)
end

@testitem "drift_shell" setup = [Share] begin
    res = drift_shell(t, x, maginput; dipole...)
    Lm = res.Lm
    @test Lm ≈ 4.326679 atol = 1.0e-2
    L_posit = _compute_dipole_L_shell(res.posit)
    @test maximum(abs, filter(!isnan, L_posit .- Lm)) / Lm <= 1.0e-2
    @test count(!isnan, res.Blocal) == sum(res.Nposit)
end

@testitem "drift_bounce_orbit" setup = [Share] begin
    res = drift_bounce_orbit(t, x, 90, maginput; dipole...)
    Lm = res.Lm
    @test Lm ≈ 4.326679 atol = 1.0e-2
    alt = x[1]
    @test abs((res.hmin - alt) / alt) <= 1.0e-2
    L_posit = _compute_dipole_L_shell(res.posit)
    @test maximum(abs, filter(!isnan, L_posit .- Lm)) / Lm <= 1.0e-2
end

@testitem "Thread safety" begin
    code = """
    using IRBEM, Dates
    t = DateTime(2015, 2, 2, 6, 12, 43)
    f(i) = make_lstar(t, [600.0 + 10i, 60.0, 50.0], (; Kp = 40.0); kext = isodd(i) ? "T89" : "OPQ77").Lm
    serial = f.(1:100)
    threaded = similar(serial)
    Threads.@threads for i in 1:100
        threaded[i] = f(i)
    end
    exit(isequal(serial, threaded) ? 0 : 1)
    """
    @test success(`$(Base.julia_cmd()) -t 4 --project=$(Base.active_project()) -e $code`)
end

@testitem "Utility functions" begin
    using Dates

    @test IRBEM.parse_kext("None") == 0
    @test IRBEM.parse_kext("OPQ77") == IRBEM.parse_kext(:opq77) == 5
    @test IRBEM.parse_kext(5) == 5
    @test_throws ArgumentError IRBEM.parse_kext("T99")

    @test IRBEM.coord_sys("GDZ") == 0
    @test IRBEM.coord_sys("gsm") == IRBEM.coord_sys(:GSM) == 2
    @test IRBEM.coord_sys("TEME") == 14
    @test_throws ArgumentError IRBEM.coord_sys("RTN")

    @test_throws ArgumentError make_lstar(DateTime(2020), [600.0, 60.0, 50.0], (; kp = 40.0))
    @test_throws ArgumentError make_lstar(DateTime(2020), [600.0, 60.0, 50.0]; options = [0, 0, 0])
end
