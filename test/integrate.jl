
@testset "Tests for integrate" begin
    
    # PolySet of p(x) = 1, x, 2x + 3x²
    ps = PolySet([ 1.0 0.0 0.0;
                    0.0 1.0 0.0;
                    0.0 2.0 3.0])

    # Integration over [0.0 1.0]
    y = integrate(ps, 0.0, 1.0)

    @test y == [1.0, 0.5, 2.0]

    # Antiderivative of ps which cancells at 0.0
    trueips = PolySet([0.0 1.0 0.0 0.0;
                        0.0 0.0 0.5 0.0;
                        0.0 0.0 1.0 1.0])

    ips = integrate(ps, 0.0)

    @test trueips == ips

    # Test with symmetrical bounds (-b, b)
    y = integrate(ps, -1.0, 1.0)
    
    @test y == [2.0, 0.0, 2.0]
end
