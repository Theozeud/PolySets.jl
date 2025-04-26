
@testset "Tests for integrate" begin
    
    # PolySet of p(x) = 1, x, 2x + 3x²
    ps = PolysSet([ 1.0 0.0 0.0;
                    0.0 1.0 0.0;
                    0.0 2.0 3.0])

    # Integration over [0.0 1.0]
    y = integrate(ps, 0.0, 1.0)

    @test y == [1.0, 0.5, 2.0]

    # Antiderivative of ps which cancells at 0.0
    trueips = PolysSet([0.0 1.0 0.0 0.0;
                        0.0 0.0 0.5 0.0;
                        0.0 0.0 1.0 1.0])

    ips = integrate(ps, 0.0)

    @test trueips == ips

    #=
    # Test avec bornes symétriques (-b, b) spécial
    ps = PolysSet([0.0, 0.0, 1.0])  # p(x) = x²
    y = integrate(ps, -1.0, 1.0)
    # ∫(x²) dx de -1 à 1 = 2 * (1/3) = 2/3
    @test isapprox(y[1], 2/3)

    # Test de la version mutating integrate! avec cache
    ps = PolysSet([0.0 2.0 0.0 4.0]) # p(x) = 2x + 4x³
    y = zeros(Float64, npolys(ps))
    cache = zeros(Float64, maxdeg(ps))
    integrate!(y, ps, 0.0, 1.0, cache)
    # ∫(2x + 4x³) dx de 0 à 1 = [x² + x⁴]₀¹ = 1 + 1 = 2
    @test isapprox(y[1], 2.0)

    # Test de l'intégrale indéfinie
    ps = PolysSet([2.0 3.0])  # p(x) = 2 + 3x
    ips = integrate(ps, 0.0)  # devrait donner p(x) = 2x + (3/2)x²
    expected_coeffs = [0.0 2.0 (3/2)]
    @test all(isapprox.(ips.coeffs[1,:], expected_coeffs))

    # Test que l'intégrale indéfinie est correctement ajustée à 0 en a
    # (en a = 0, valeur doit être 0)
    yval = zeros(Float64, npolys(ips))
    evaluate!(yval, ips, 0.0)
    @test isapprox(yval[1], 0.0)

    # Test avec un PolysSet vide
    ps = allocate_PolysSet(Float64, 0, 2)
    y = integrate(ps, 0.0, 1.0)
    @test length(y) == 0

    # Test intégrale avec types mixtes (Int et Float)
    ps = PolysSet([1 0]) # p(x) = 1
    y = integrate(ps, 0, 2.0)
    @test isapprox(y[1], 2.0)
    =#
end
