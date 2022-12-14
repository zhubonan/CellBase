using Test
using CellBase: Composition, Cell, Lattice, formula
using LaTeXStrings


@testset "Composition" begin
    comp = Composition(:C => 1, :O => 2)
    @test comp[:C] == 1.0
    @test comp[:O] == 2.0

    comp = Composition("CO2")
    @test comp[:C] == 1.0
    @test comp[:O] == 2.0

    comp = Composition("(CO2)3(H2O)2")
    @test comp[:C] == 3.0
    @test comp[:O] == 8.0
    @test comp[:H] == 4.0

    cell = Cell(Lattice(10, 10, 10), [1, 1], [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]])
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

    c1 = Composition("C3O4")
    @test formula(c1) == :C3O4
    show(devnull, c1)

    c1 = Composition("C3O")
    @test formula(c1) == :C3O

    c1 = Composition("C3O4")
    @test latex_formula(c1) == L"$\mathrm{C_{3}O_{4}}$"

    c1 = Composition("C3O4")
    @test latex_formula(c1) == L"$\mathrm{C_{3}O_{4}}$"
    c1 = Composition("C2O")
    @test latex_formula(c1) == L"$\mathrm{C_{2}O}$"

    # Equality test
    @test Composition(formula(c1)) == c1

    # Hash
    @test hash(Composition(:C10O10)) == hash(Composition(:C10O10))

    # Settings index
    c1 = Composition("C3O4")
    c1[:C] = 2
    @test c1 == Composition(:C2O4)
    c1[:Zn] = 2
    @test c1 == Composition(:C2O4Zn2)

    @test contains(Composition(:C2O8), Composition(:CO))

    @test Composition(:C2O2) - Composition(:O) == Composition(:C2O)

    # Obtain atomic weight
    CellBase.atomic_weight(c1)
end