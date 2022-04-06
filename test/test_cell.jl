using Test
using CellBase

@testset "Cell" begin
    mat = Float64[
        2 0 0 
        0 2 0
        0 0 3
    ]
    s = Cell(Lattice(mat), [1,1,1,1], rand(3, 4))
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
    @testset "Methods" begin
        @test nions(s) == 4
        @test size(positions(s)) == (3, 4) 
        @test volume(s) > 0
        @test all(x-> x == :H, species(s)) 
        ss = CellBase.make_supercell(s, 2, 2, 2)
        @test CellBase.natoms(ss) == CellBase.natoms(s) * 8
        @test CellBase.cellpar(ss)[1:3] == CellBase.cellpar(s)[1:3] .* 2
        @test CellBase.cellpar(ss)[4:6] == CellBase.cellpar(s)[4:6]
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
    end
end

@testset "Lattice" begin
    @testset "Construct" begin
        @test (Lattice(rand(3, 3)); true)
    end
    @testset "Methods" begin
        l = Lattice(Float64[2 0 0
                     0 1 0
                     0 0 1])
        l2 = Lattice(Float64[2 0 0
                     2 1 0
                     0 0 1])
        @test volume(l) == 2.0
        @test cellmat(l) == l.matrix
        @test cellvecs(l)[1] == [2, 0, 0]
        @test cellvecs(l2)[1] == [2, 2, 0]
        site = Site([2.1, 0, 0], 1, :H)
        wrap!(site, l)
        @test site.position[1] ≈ 0.1
        @test cellpar(l) == [2, 1, 1, 90, 90, 90]
    end

    @testset "Mic" begin
        l = Lattice(Float64[2 0 0
                     0 1 0
                     0 0 1])
        vec, d = CellBase.mic(l, reshape([0.1, 0.1, 0.1], 3, 1))
        @test d[1] ≈ sqrt(3) / 10
        @test vec[:] == [0.1, 0.1, 0.1]

        vec, d = CellBase.mic(l, reshape([1, 1, 1], 3, 1))
        @test d[1] != sqrt(3) / 10
        @test vec[:] != [1, 1, 1]
    end

end