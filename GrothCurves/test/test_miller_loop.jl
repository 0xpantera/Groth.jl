using Test
using GrothCurves
using GrothAlgebra

@testset "Miller Loop Operations" begin

    @testset "Line coefficient construction" begin
        # Test that we can create LineCoeffs with D-twist structure
        a = Fp2Element(1, 2)  # y coefficient
        b = Fp2Element(3, 4)  # x coefficient
        c = Fp2Element(5, 6)  # constant term

        coeffs = LineCoeffs(a, b, c)
        @test coeffs.a == a
        @test coeffs.b == b
        @test coeffs.c == c
    end

    @testset "Doubling step" begin
        # Test with generators
        Q = g2_generator()

        # Compute doubling step
        T_doubled, line_coeffs = doubling_step(Q)

        # Check that we got a G2Point back
        @test isa(T_doubled, G2Point)
        @test isa(line_coeffs, LineCoeffs)

        # Verify T_doubled is 2Q
        Q2_direct = Q + Q
        @test T_doubled == Q2_direct

        # Test with point at infinity
        Q_inf = zero(G2Point)
        T_inf, coeffs_inf = doubling_step(Q_inf)
        @test iszero(T_inf)

        # Check line coefficients have correct structure
        @test isa(line_coeffs.a, Fp2Element)  # y coefficient
        @test isa(line_coeffs.b, Fp2Element)  # x coefficient
        @test isa(line_coeffs.c, Fp2Element)  # constant
    end

    @testset "Addition step" begin
        # Create test points
        Q = g2_generator()
        Q2 = Q + Q  # 2Q

        # Test addition: Q + Q should give 2Q
        T_sum, line_coeffs = addition_step(Q, Q)
        @test isa(T_sum, G2Point)
        @test isa(line_coeffs, LineCoeffs)
        # Note: addition_step(Q, Q) actually calls doubling internally
        @test T_sum == Q2

        # Test Q + 2Q = 3Q
        Q3 = Q + Q + Q
        T_sum2, line_coeffs2 = addition_step(Q, Q2)
        @test T_sum2 == Q3

        # Test with identity
        Q_inf = zero(G2Point)
        T_id, coeffs_id = addition_step(Q, Q_inf)
        @test T_id == Q

        # Test identity + Q = Q
        T_id2, coeffs_id2 = addition_step(Q_inf, Q)
        @test T_id2 == Q
    end

    @testset "Line evaluation" begin
        # Test that evaluate_line produces an Fp12Element
        coeffs = LineCoeffs(
            Fp2Element(1, 2),
            Fp2Element(3, 4),
            Fp2Element(5, 6)
        )

        P = g1_generator()
        f = evaluate_line(coeffs, P)
        @test isa(f, Fp12Element)

        # Test with point at infinity
        P_inf = zero(G1Point)
        f_inf = evaluate_line(coeffs, P_inf)
        @test isa(f_inf, Fp12Element)

        # Check D-twist sparse structure (0,3,4)
        # c0 should have non-zero only at position 0
        # c1 should have non-zero at positions 0 and 1
        c0 = f[1]  # First Fp6 component
        c1 = f[2]  # Second Fp6 component

        # For D-twist with P != 0, we expect:
        # c0[1] = a*yP (non-zero), c0[2] = 0, c0[3] = 0
        # c1[1] = b*xP (non-zero), c1[2] = c (non-zero), c1[3] = 0
        @test c0[2] == zero(Fp2Element)
        @test c0[3] == zero(Fp2Element)
        @test c1[3] == zero(Fp2Element)
    end

    @testset "Miller loop basic" begin
        # Test with identity elements
        P = zero(G1Point)
        Q = zero(G2Point)

        f = miller_loop(P, Q)
        @test isa(f, Fp12Element)
        @test isone(f)  # Miller loop of identity should be 1

        # Test with generators
        P_gen = g1_generator()
        Q_gen = g2_generator()

        f_gen = miller_loop(P_gen, Q_gen)
        @test isa(f_gen, Fp12Element)
        @test !iszero(f_gen)  # Should be non-zero
        @test !isone(f_gen)   # Should not be 1 for non-zero inputs
    end

    @testset "Miller loop properties" begin
        # Test that miller_loop(0, Q) = 1
        P_zero = zero(G1Point)
        Q = g2_generator()

        f = miller_loop(P_zero, Q)
        @test isone(f)

        # Test that miller_loop(P, 0) = 1
        P = g1_generator()
        Q_zero = zero(G2Point)

        f = miller_loop(P, Q_zero)
        @test isone(f)

        # Test that miller_loop produces consistent results
        P = g1_generator()
        Q = g2_generator()

        f1 = miller_loop(P, Q)
        f2 = miller_loop(P, Q)
        @test f1 == f2  # Same inputs should give same output
    end

    @testset "NAF representation" begin
        # Test that NAF is correctly defined
        @test isa(ATE_LOOP_COUNT_NAF, Vector{Int8})
        @test length(ATE_LOOP_COUNT_NAF) == 65

        # NAF should only contain -1, 0, 1
        for digit in ATE_LOOP_COUNT_NAF
            @test digit in [-1, 0, 1]
        end

        # The last digit should be 1 (MSB of u)
        @test ATE_LOOP_COUNT_NAF[end] == 1
    end

    @testset "Frobenius endomorphism" begin
        # Test p-power endomorphism
        Q = g2_generator()

        # Apply Frobenius
        Q_pi = frobenius_g2(Q, 1)
        @test isa(Q_pi, G2Point)
        @test is_on_curve(Q_pi)

        # Apply Frobenius twice
        Q_pi2 = frobenius_g2(Q, 2)
        @test isa(Q_pi2, G2Point)
        @test is_on_curve(Q_pi2)

        # Composition should work
        Q_pi_then_pi = frobenius_g2(Q_pi, 1)
        @test Q_pi_then_pi == Q_pi2

        # Test with identity
        Q_inf = zero(G2Point)
        Q_inf_pi = frobenius_g2(Q_inf, 1)
        @test iszero(Q_inf_pi)
    end
end
