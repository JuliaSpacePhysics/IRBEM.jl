@testitem "find_mirror_point" setup = [Share] begin
    true_Blocal = 42271.43059990003
    true_Bmin = 42271.43059990003
    true_POSIT = [0.35282136776620165, 0.4204761325793738, 0.9448914452448274]

    alpha = 90.0  # Local pitch angle in degrees
    @test approx(find_mirror_point(t, x, alpha, maginput; kext), (true_Blocal, true_Bmin, true_POSIT))
    multi = find_mirror_point([t, t], [x, x], alpha, maginput; coord = "GDZ", kext)
    @test multi.posit ≈ [true_POSIT true_POSIT]
end

@testitem "find_magequator" setup = [Share] begin
    true_Bmin = 626.2258295723121
    true_XGEO = [2.1962220856733894, 2.8360222891612192, 0.3472455620354017]
    @test approx(find_magequator(t, x, maginput; kext), (true_Bmin, true_XGEO))
    multi = find_magequator(ts, xs, maginput_array; kext)
    @test approx(multi, (fill(true_Bmin, n), repeat(true_XGEO, 1, n)))
end

@testitem "find_foot_point" setup = [Share] begin
    _foot_point_true = (;
        XFOOT = [99.99412846343064, 61.113869939535036, 50.55633537632344],
        BFOOT = [-25644.012241653385, -25370.689449132995, -38649.994779664776],
        BFOOTMAG = 52868.793663583165,
    )
    stopAlt = 100
    hemiFlag = 0
    @test approx(find_foot_point(t, x, stopAlt, hemiFlag, maginput; kext), _foot_point_true)
end
