@testitem "find_mirror_point" setup = [Share] begin
    true_Blocal = 42271.43059990003
    true_Bmin = 42271.43059990003
    true_POSIT = [0.35282136776620165, 0.4204761325793738, 0.9448914452448274]

    alpha = 90.0  # Local pitch angle in degrees
    result = find_mirror_point(model, X, alpha, maginput)

    @test result[1] == true_Blocal
    @test result[2] == true_Bmin
    @test result[3] == true_POSIT
    @test result == find_mirror_point(t, x, alpha, "GDZ", maginput; kext = "T89")
end

@testitem "find_magequator" setup = [Share] begin
    true_Bmin = 626.2258295723121
    true_XGEO = [2.1962220856733894, 2.8360222891612192, 0.3472455620354017]
    result = find_magequator(model, X, maginput)

    @test result[1] == true_Bmin
    @test result[2] == true_XGEO
end

@testitem "find_foot_point" setup = [Share] begin
    _foot_point_true = (;
        XFOOT = [99.99412846343064, 61.113869939535036, 50.55633537632344],
        BFOOT = [-25644.012241653385, -25370.689449132995, -38649.994779664776],
        BFOOTMAG = 52868.793663583165,
    )
    stopAlt = 100
    hemiFlag = 0
    @test find_foot_point(model, X, stopAlt, hemiFlag, maginput) == _foot_point_true
    @test find_foot_point(t, x, stopAlt, hemiFlag, "GDZ", maginput; kext = "T89") == _foot_point_true
end
