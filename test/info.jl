@testitem "Library information functions" begin
    @test IRBEM.get_igrf_version() == 14
    @test IRBEM.irbem_fortran_version() > 0
    @test !isempty(IRBEM.irbem_fortran_release())
end
