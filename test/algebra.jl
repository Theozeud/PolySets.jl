@testset "Algebra" begin

    # PolySet of p(x) = 1, x, 2x + 3x²
    ps = PolySet([ 1.0 0.0 0.0;
                    0.0 1.0 0.0;
                    0.0 2.0 3.0])

    c = 4.0


    # ps + c
    ps_plus_c = PolySet([5.0 0.0 0.0;
                         4.0 1.0 0.0;
                         4.0 2.0 3.0])

    @test ps + c == c + ps == ps_plus_c

    # ps - c   
    ps_minus_c = PolySet([-3.0 0.0 0.0;
                          -4.0 1.0 0.0;
                          -4.0 2.0 3.0])

    @test ps - c == ps_minus_c
    
    # c - ps   
    c_minus_ps = PolySet([3.0 0.0 0.0;
                         4.0 -1.0 0.0;
                         4.0 -2.0 -3.0])

    @test c - ps == c_minus_ps

    # ps * c
    c_times_ps = PolySet([4.0 0.0 0.0;
                          0.0 4.0 0.0;
                          0.0 8.0 12.0])

    @test c * ps == c_times_ps

    # ps / c
    ps_div_c = PolySet([0.25 0.0 0.0;
                        0.0  0.25 0.0;
                        0.0  0.5 0.75])

    @test ps / c == ps_div_c

end