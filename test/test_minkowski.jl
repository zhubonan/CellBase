using Test
using CellBase
@testset "Minkovski Reduce" begin
    mat = [1 2 1; 2 2 2; -1 5 3]
    reduced, H = CellBase.minkowski_reduce3d(collect(transpose(mat)))
    @test reduced == [
        1 0 -2
        0 2 1
        1 0 2
    ]
    @test H == [
        -1 2 -3
        1 -1 1
        0 0 1
    ]
    mat = [1 0 0; 0 2 0; 0 0 1]
    reduced, H = CellBase.minkowski_reduce3d(mat)
    @test reduced == [
        1 0 0
        0 0 2
        0 1 0
    ]
end
