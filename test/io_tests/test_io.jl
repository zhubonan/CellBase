#=
Tests for IO related codes
=#

using CellBase
using Test

@testset "IO" begin
    this_dir = splitpath(@__FILE__)[1:end-1]
    @testset "SHELX" begin
        fpath = joinpath(this_dir..., "lno.res")
        cell = CellBase.read_res(fpath)
        @test species(cell) == [:Li, :Li, :O, :O, :O, :O, :Ni, :Ni]
        @test cellpar(cell)[1] == 2.73724
        @test cellpar(cell)[2] == 6.02753

        # Read packed res
        fpath = joinpath(this_dir..., "lno.pack.res")
        cells = CellBase.read_res_many(fpath)
        @test length(cells) == 2
        @test cellpar(cells[1])[1] == 2.73724
        @test cellpar(cells[1])[2] == 6.02753
    end

    @testset "CELL" begin
        cell = CellBase.read_cell(joinpath(this_dir..., "Si2.cell"))
        @test species(cell) == [:Si, :Si]
        @test cellmat(cell)[1] ≈ 2.69546455 atol = 1e-5

        # Test writing CEll file
        pos = rand(3, 3)
        latt = rand(3, 3)
        tmp, f = mktemp()
        CellBase.CellIO.write_cell(f, latt, pos, [:H, :H, :H])
        close(f)
        celltmp = CellBase.read_cell(tmp)
        @test all(cellmat(celltmp) .== latt)
        @test all(positions(celltmp) .== pos)
        rm(tmp)

        CellBase.write_cell(tmp, celltmp)
        celltmp2 = CellBase.read_cell(tmp)
        @test all(cellmat(celltmp2) .== cellmat(celltmp))
    end

    @testset "CASTEP" begin
        snapshots =
            CellBase.read_castep(joinpath(this_dir..., "Fe.castep"), only_first=false)
        @test length(snapshots) == 11
        snap = snapshots[1]
        @test snap.forces[1] ≈ -3.20036
    end

end
