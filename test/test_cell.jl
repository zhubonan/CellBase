using Test
using CellBase
import Spglib

@testset "Cell" begin
    mat = Float64[
        2 0 0
        0 2 0
        0 0 3
    ]
    example_cell = Cell(Lattice(mat), [1, 1, 1, 1], rand(3, 4) .* 10 .- 5)
    example_cell2 = Cell(Lattice(mat), [1, 2, 2, 1], rand(3, 4) .* 10 .- 5)
    example_cell2.arrays[:forces] = rand(3, 3, 4)

    @testset "Construct" begin
        @test begin
            Cell(Lattice(mat), [1, 1, 1, 1], rand(3, 4))
            true
        end
        @test begin
            Cell(Lattice(mat), [:H, :H, :H, :H], rand(3, 4))
            true
        end
    end

    @testset "Sort" begin
        c = Cell(Lattice(mat), [:H, :O, :O, :H], rand(3, 4))
        c.arrays[:forces] = rand(3, 3, 4)

        c_new = sort(c)
        @test species(c_new) == [:H, :H, :O, :O]
        sort!(c)
        @test species(c) == [:H, :H, :O, :O]
    end

    @testset "Methods" begin
        @test nions(example_cell) == 4
        @test size(positions(example_cell)) == (3, 4)
        @test volume(example_cell) > 0
        @test all(x -> x == :H, species(example_cell))
        @test all(x -> x == 1, atomic_numbers(example_cell))
        ss = CellBase.make_supercell(example_cell, 2, 2, 2)
        @test CellBase.natoms(ss) == CellBase.natoms(example_cell) * 8
        @test CellBase.cellpar(ss)[1:3] == CellBase.cellpar(example_cell)[1:3] .* 2
        @test CellBase.cellpar(ss)[4:6] == CellBase.cellpar(example_cell)[4:6]

        @test get_cellmat(example_cell) !== cellmat(example_cell)
        @test get_positions(example_cell) !== positions(example_cell)
        @test get_lattice(example_cell) !== lattice(example_cell)

        # Test getting wrapped array
        @test length(CellBase.sposarray(example_cell)) == length(example_cell)
        @test begin
            pos = CellBase.wrapped_spos(example_cell)
            all(all(x .< 3) for x in pos)
        end

        wrap!(example_cell)
        @test all(all(x .< 3) for x in CellBase.sposarray(example_cell))
    end

    @testset "interface" begin
        @test length(example_cell) == 4
        @test length(example_cell[[1, 2]]) == 2
        @test isa(example_cell[1], CellBase.Site)
        stmp = deepcopy(example_cell)
        site = stmp[1]
        site.position[1] = -10.0
        @test stmp.positions[1, 1] == -10.0

        # Test clipping
        tmp = example_cell2[[1, 3]]
        @test length(tmp) == 2
        @test size(tmp.arrays[:forces]) == (3, 3, 2)
        @test tmp.arrays[:forces][:, :, 1] == example_cell2.arrays[:forces][:, :, 1]
        @test tmp.arrays[:forces][:, :, 2] == example_cell2.arrays[:forces][:, :, 3]
    end
end

@testset "Site" begin
    @testset "Construct" begin
        @test (Site(rand(3), 1, :H); true)
    end
    @testset "methods" begin
        s1 = Site([0.0, 0.0, 0.0], 1, :H)
        s2 = Site([1.0, 1.0, -1.0], 1, :H)
        @test distance_between(s1, s2) ≈ sqrt(3)
        @test distance_squared_between(s1, s2) ≈ 3
        @test distance_between(s1, s2, [0, 0, 1]) ≈ sqrt(2)
        @test distance_squared_between(s1, s2, [0, 0, 1]) ≈ 2
        @test s1.x == 0.0
        @test s2.y == 1.0
        @test s2.z == -1.0
        @test :x in propertynames(s1)
    end
end

@testset "Lattice" begin
    @testset "Construct" begin
        @test (Lattice(rand(3, 3)); true)
    end
    @testset "Methods" begin
        l = Lattice(Float64[
            2 0 0
            0 1 0
            0 0 1
        ])
        l2 = Lattice(Float64[
            2 0 0
            2 1 0
            0 0 1
        ])
        @test volume(l) == 2.0
        @test cellmat(l) === l.matrix
        @test get_cellmat(l) !== l.matrix

        @test cellvecs(l)[1] == [2, 0, 0]
        @test cellvecs(l2)[1] == [2, 2, 0]
        site = Site([2.1, 0, 0], 1, :H)
        wrap!(site, l)
        @test site.position[1] ≈ 0.1
        @test cellpar(l) == [2, 1, 1, 90, 90, 90]
    end

    @testset "Mic" begin
        l = Lattice(Float64[
            2 0 0
            0 1 0
            0 0 1
        ])
        vec, d = CellBase.mic(l, reshape([0.1, 0.1, 0.1], 3, 1))
        @test d[1] ≈ sqrt(3) / 10
        @test vec[:] == [0.1, 0.1, 0.1]

        vec, d = CellBase.mic(l, reshape([1, 1, 1], 3, 1))
        @test d[1] != sqrt(3) / 10
        @test vec[:] != [1, 1, 1]
    end


end

@testset "NL" begin
    lattice = Lattice(10, 10, 10)
    frac_pos = [
        0 0.5 0.5 0.5 0
        0 0.5 0.5 0.0 0.5
        0 0.5 0.0 0.5 0.5
    ]
    pos = cellmat(lattice) * frac_pos
    testcell = Cell(Lattice(10, 10, 10), [1, 2, 3, 4, 5], pos)

    parray = ExtendedPointArray(testcell, 6)
    @test nions_orig(parray) == 5
    @test nions_extended(parray) == 5 * 27

    # Building neighbour list
    nl = NeighbourList(parray, 8.0)
    @test num_neighbours(nl, 1) == 12
    @test num_neighbours(nl, 2) == 6
    @test num_neighbours(nl, 3) == 14

    # Distance
    list = collect(eachneighbour(nl, 1))
    @test length(list) == 12
    @test list[1][1] == 3
    @test list[1][2] == 64
    @test list[1][3] ≈ sqrt(50)

    list = collect(eachneighbour(nl, 1, unique=true))
    @test list[1][1] == 3
    @test list[1][2] == 64
    @test list[1][3] ≈ sqrt(50)
    length(list) == 3

    # Outside the cell
    positions(testcell)[1] += 101
    parray = ExtendedPointArray(testcell, 6)
    @test parray.orig_positions[1][1] != positions(testcell)[1]
    @test all(all(x .< 10) for x in parray.orig_positions)
end

@testset "Spglib" begin
    lattice = Lattice(10, 10, 10)
    frac_pos = [
        0 0.5 0.5 0.5 0
        0 0.5 0.5 0.0 0.5
        0 0.5 0.0 0.5 0.5
    ]
    pos = cellmat(lattice) * frac_pos
    testcell = Cell(Lattice(10, 10, 10), [1, 2, 3, 3, 3], pos)
    @test all(CellBase.SCell(testcell).lattice .== cellmat(lattice))
    scell = Spglib.Cell(testcell)
    cell2 = Cell(scell)
    @test all(cell2.positions .== testcell.positions)

    @test Spglib.get_international(testcell) == "Pm-3m"
    Spglib.get_dataset(testcell)
    Spglib.get_symmetry(testcell)

    testcell.metadata[:test] = 1
    std = Spglib.standardize_cell(testcell)
    prim = Spglib.find_primitive(testcell)

    @testset "niggli_reduce_cell" begin
        testcell = Cell(Lattice(5.0, 5.0, 5.0, 30.0, 30.0, 30.0), [1, 2, 3, 3, 3], pos)
        reduced_cell = niggli_reduce_cell(testcell; wrap_pos=false)
        @test all(positions(testcell) .== positions(reduced_cell))
        reduced_cell2 = niggli_reduce_cell(testcell; wrap_pos=true)
        CellBase.wrap!(reduced_cell)
        @test all(
            get_scaled_positions(reduced_cell) .== get_scaled_positions(reduced_cell2),
        )
    end

    @testset "niggli_reduce" begin
        testcell = Cell(Lattice(5.0, 5.0, 5.0, 30.0, 30.0, 30.0), [1, 2, 3, 3, 3], pos)
        testtmp = testcell
        reduced_cell = niggli_reduce(testcell)
        reduced_cell2 = niggli_reduce_cell(testcell; wrap_pos=false)
        @test all(
            isapprox.(
                get_scaled_positions(reduced_cell),
                get_scaled_positions(reduced_cell2),
                atol=1e-7,
            ),
        )
    end

    @testset "roundtrip functions" begin
        testcell = Cell(Lattice(5.0, 5.0, 5.0, 30.0, 30.0, 30.0), [1, 2, 3, 3, 3], pos)
        @test isa(standardize_cell(testcell), Cell)
        @test isa(find_primitive(testcell), Cell)
        @test isa(refine_cell(testcell), Cell)

    end

end
