@testitem "Library information functions" begin
    @test IRBEM.get_igrf_version() == 14
    @test_nowarn IRBEM.irbem_fortran_version()
    @test_nowarn IRBEM.irbem_fortran_release()
end
