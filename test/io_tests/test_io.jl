#=
Tests for IO related codes
=#

using CellBase
using Test

@testset "IO" begin
    this_dir = splitpath(@__FILE__)[1:end-1]
    @testset "SHELX" begin
        cell = CellBase.read_res(joinpath(this_dir..., "lno.res"))
        @test species(cell) == [:Li, :Li, :O, :O, :O, :O, :Ni, :Ni]
        @test cellpar(cell)[1] == 2.73724
        @test cellpar(cell)[2] == 6.02753
    end

    @testset "CELL" begin
        cell = CellBase.read_cell(joinpath(this_dir..., "Si2.cell"))
        @test species(cell) == [:Si, :Si]
        @test cellmat(cell)[1] ≈ 2.69546455  atol=1e-5
    end

    @testset "CASTEP" begin
        snapshots = CellBase.read_castep(joinpath(this_dir..., "Fe.castep"), only_first=false)
        @test length(snapshots) == 11
        snap = snapshots[1]
        @test snap.forces[1] ≈ -3.20036 
    end
end