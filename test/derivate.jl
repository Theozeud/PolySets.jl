
@testset "Tests for derivate" begin
    
    # PolySet of p(x) = 1, x, 2x + 3x²
    ps = PolySet([ 1.0 0.0 0.0;
                    0.0 1.0 0.0;
                    0.0 2.0 3.0])

    truedps = PolySet([0.0 0.0;
                        1.0 0.0;
                        2.0 6.0])

    dps = derivate(ps)

    @test truedps == dps

end
