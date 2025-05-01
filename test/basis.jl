@testset "Basis" begin

    # TEST ALLOCATION OF PolySet
    zp  = zeros(PolySet{Int}, 3, 7)
    zp2 = zero(zp)
    zp3 = allocate_PolySet(Int, 3,6)

    @test zp == zp2 == zp3
    @test iszero(zp.coeffs)

    @test zp[1] == [0, 0, 0, 0, 0, 0, 0]
    zp[1] = [0, 0, 1, 0, 0, 0, 0]
    zp[1,2] = 9
    zp[3,7] = 10
    @test first(zp) == [0, 9, 1, 0, 0, 0, 0]    
    @test last(zp)  == [0, 0, 0, 0, 0, 0, 10]  

    # TEST MONOMIALS WITH FLOAT64
    pm = Monomials(4)

    @test eltype(pm)    == Float64
    @test size(pm)      == (4,4)
    @test size(pm, 1)   == 4
    @test size(pm, 2)   == 4
    @test length(pm)    == 4
    @test npolys(pm)    == 4
    @test maxdeg(pm)    == 3

    pmeval = evaluate(pm, [0.0, 1.0, 2.0])

    

    @test pmeval == [1.0 1.0 1.0;
                     0.0 1.0 2.0;
                     0.0 1.0 4.0;
                     0.0 1.0 8.0]

    # TEST MONOMIALS WITH RATIONAL{INT}
    pm = Monomials(4,Rational{Int})

    @test eltype(pm)    == Rational{Int}
    @test size(pm)      == (4,4)
    @test size(pm, 1)   == 4
    @test size(pm, 2)   == 4
    @test length(pm)     == 4
    @test npolys(pm)    == 4
    @test maxdeg(pm)    == 3

    pmeval = evaluate(pm, [0//1, 1//1, 1//2])

    @test pmeval == [1//1 1//1 1//1;
                     0//1 1//1 1//2;
                     0//1 1//1 1//4;
                     0//1 1//1 1//8]

    # TESTS INTLEGENDRE

end