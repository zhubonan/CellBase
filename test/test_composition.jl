using Test
using CellBase: Composition, Cell, Lattice


@testset "Composition" begin
    comp = Composition(:C => 1, :O=>2)
    @test comp[:C] == 1.0
    @test comp[:O] == 2.0

    comp = Composition("CO2")
    @test comp[:C] == 1.0
    @test comp[:O] == 2.0

    comp = Composition("(CO2)3(H2O)2")
    @test comp[:C] == 3.0
    @test comp[:O] == 8.0
    @test comp[:H] == 4.0

    cell = Cell(Lattice(10, 10, 10), [1, 1], [[0., 0., 0.], [1., 1., 1.]])
    comp = Composition(cell)
    @test comp[:H] == 2.0

    c1 = Composition("C3O4")
    c2 = Composition("O2F3")
    c3 = c1 + c2
    @test c3[:O] == 6.0
    @test c3[:F] == 3.0
    @test c3[:C] == 3.0

    c1 = Composition("C3O4") * 2
    @test c1[:C] == 6.0

    c1 = Composition("C3O4") / 2
    @test c1[:C] == 1.5
end