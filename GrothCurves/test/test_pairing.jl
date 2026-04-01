using Test
using GrothCurves
using GrothAlgebra

@testset "BN254 Pairing Tests" begin

    # Get generators for testing
    P = g1_generator()
    Q = g2_generator()

    @testset "Identity and Degeneracy" begin
        # Test identity elements
        O_G1 = zero(G1Point)
        O_G2 = zero(G2Point)

        @test pairing(O_G1, Q) == one(Fp12Element)
        @test pairing(P, O_G2) == one(Fp12Element)
        @test pairing(O_G1, O_G2) == one(Fp12Element)

        # Test non-degeneracy
        e_P_Q = pairing(P, Q)
        @test e_P_Q != one(Fp12Element)
    end

    @testset "Negation Properties" begin
        # Test that negation works correctly
        neg_P = -P
        neg_Q = -Q

        e_P_Q = pairing(P, Q)
        e_negP_Q = pairing(neg_P, Q)
        e_P_negQ = pairing(P, neg_Q)
        e_negP_negQ = pairing(neg_P, neg_Q)

        @test e_negP_Q == inv(e_P_Q)
        @test e_P_negQ == inv(e_P_Q)
        @test e_negP_negQ == e_P_Q
    end

    @testset "Bilinearity - Scalar Multiplication" begin
        # Test e([n]P, Q) = e(P, [n]Q) = e(P, Q)^n for various n

        # Test with n = 2
        P2 = P + P
        Q2 = Q + Q

        e_P_Q = pairing(P, Q)
        e_2P_Q = pairing(P2, Q)
        e_P_2Q = pairing(P, Q2)
        e_P_Q_squared = e_P_Q^2

        @test e_2P_Q == e_P_2Q
        @test e_2P_Q == e_P_Q_squared
        @test e_P_2Q == e_P_Q_squared

        # Test with n = 3
        P3 = P + P + P
        Q3 = Q + Q + Q

        e_3P_Q = pairing(P3, Q)
        e_P_3Q = pairing(P, Q3)
        e_P_Q_cubed = e_P_Q^3

        @test e_3P_Q == e_P_3Q
        @test e_3P_Q == e_P_Q_cubed
        @test e_P_3Q == e_P_Q_cubed

        # Test with n = 4, 5, 6, 7
        for n in 4:7
            # Compute [n]P and [n]Q
            nP = P
            for i in 2:n
                nP = nP + P
            end

            nQ = Q
            for i in 2:n
                nQ = nQ + Q
            end

            # Compute pairings
            e_nP_Q = pairing(nP, Q)
            e_P_nQ = pairing(P, nQ)
            e_P_Q_n = e_P_Q^n

            @test e_nP_Q == e_P_nQ
            @test e_nP_Q == e_P_Q_n
        end
    end

    @testset "Bilinearity - Addition Formula" begin
        # Test e(P1 + P2, Q) = e(P1, Q) * e(P2, Q)

        P2 = P + P
        P3 = P + P + P

        # P + 2P = 3P
        e_P_plus_2P_Q = pairing(P3, Q)
        e_P_Q = pairing(P, Q)
        e_2P_Q = pairing(P2, Q)
        product = e_P_Q * e_2P_Q

        @test e_P_plus_2P_Q == product

        # Test e(P, Q1 + Q2) = e(P, Q1) * e(P, Q2)
        Q2 = Q + Q
        Q3 = Q + Q + Q

        # Q + 2Q = 3Q
        e_P_Q_plus_2Q = pairing(P, Q3)
        e_P_Q = pairing(P, Q)
        e_P_2Q = pairing(P, Q2)
        product2 = e_P_Q * e_P_2Q

        @test e_P_Q_plus_2Q == product2
    end

    @testset "Cross Products" begin
        # Test e([m]P, [n]Q) = e(P, Q)^(m*n)

        P2 = P + P
        P3 = P + P + P
        Q2 = Q + Q
        Q3 = Q + Q + Q

        e_P_Q = pairing(P, Q)

        # Test e([2]P, [3]Q) = e(P, Q)^6
        e_2P_3Q = pairing(P2, Q3)
        e_P_Q_to_6 = e_P_Q^6
        @test e_2P_3Q == e_P_Q_to_6

        # Test e([3]P, [2]Q) = e(P, Q)^6
        e_3P_2Q = pairing(P3, Q2)
        @test e_3P_2Q == e_P_Q_to_6

        # Test e([2]P, [2]Q) = e(P, Q)^4
        e_2P_2Q = pairing(P2, Q2)
        e_P_Q_to_4 = e_P_Q^4
        @test e_2P_2Q == e_P_Q_to_4

        # Test e([3]P, [3]Q) = e(P, Q)^9
        e_3P_3Q = pairing(P3, Q3)
        e_P_Q_to_9 = e_P_Q^9
        @test e_3P_3Q == e_P_Q_to_9
    end

    @testset "Miller Loop and Final Exponentiation" begin
        # Test that pairing = final_exp(miller_loop)

        f = miller_loop(P, Q)
        e = final_exponentiation(f)
        e_direct = pairing(P, Q)

        @test e == e_direct

        # Test with different points
        P2 = P + P
        Q2 = Q + Q

        f2 = miller_loop(P2, Q2)
        e2 = final_exponentiation(f2)
        e2_direct = pairing(P2, Q2)

        @test e2 == e2_direct
    end

    @testset "Cyclotomic Helpers" begin
        P_ref = scalar_mul(P, BigInt(11))
        Q_ref = scalar_mul(Q, BigInt(13))
        m = final_exponentiation_easy(miller_loop(P_ref, Q_ref))

        @test GrothCurves.cyclotomic_inverse(m) == conjugate(m)
        @test GrothCurves.cyclotomic_square(m) == square(m)
        @test GrothCurves.cyclotomic_exp(m, 11) == m^11
        @test exp_by_u(m) == m^GrothCurves.BN254_U
    end

    @testset "Pairing Consistency" begin
        # Test that repeated computations give the same result

        e1 = pairing(P, Q)
        e2 = pairing(P, Q)
        e3 = pairing(P, Q)

        @test e1 == e2
        @test e2 == e3

        # Test with different points
        P_alt = P + P + P + P + P  # [5]P
        Q_alt = Q + Q + Q + Q + Q + Q + Q  # [7]Q

        e_alt1 = pairing(P_alt, Q_alt)
        e_alt2 = pairing(P_alt, Q_alt)

        @test e_alt1 == e_alt2
    end

    @testset "Edge Cases" begin
        # Test with larger scalar multiplications

        # Compute [12]P and [13]Q
        P12 = P
        for i in 2:12
            P12 = P12 + P
        end

        Q13 = Q
        for i in 2:13
            Q13 = Q13 + Q
        end

        e_P_Q = pairing(P, Q)
        e_12P_13Q = pairing(P12, Q13)
        e_P_Q_to_156 = e_P_Q^156  # 12 * 13 = 156

        @test e_12P_13Q == e_P_Q_to_156

        # Test commutativity of scalar multiplication
        # e([12]P, Q) = e(P, [12]Q)
        Q12 = Q
        for i in 2:12
            Q12 = Q12 + Q
        end

        e_12P_Q = pairing(P12, Q)
        e_P_12Q = pairing(P, Q12)

        @test e_12P_Q == e_P_12Q
    end

    @testset "Random Point Tests" begin
        # Create some pseudo-random points by scalar multiplication

        # Use different primes as scalars for variety
        P17 = P
        for i in 2:17
            P17 = P17 + P
        end

        Q23 = Q
        for i in 2:23
            Q23 = Q23 + Q
        end

        P29 = P
        for i in 2:29
            P29 = P29 + P
        end

        Q31 = Q
        for i in 2:31
            Q31 = Q31 + Q
        end

        # Test bilinearity with these points
        # e(P17 + P29, Q) = e(P17, Q) * e(P29, Q)
        P46 = P17 + P29  # [46]P = [17]P + [29]P

        e_P46_Q = pairing(P46, Q)
        e_P17_Q = pairing(P17, Q)
        e_P29_Q = pairing(P29, Q)

        @test e_P46_Q == e_P17_Q * e_P29_Q

        # e(P, Q23 + Q31) = e(P, Q23) * e(P, Q31)
        Q54 = Q23 + Q31  # [54]Q = [23]Q + [31]Q

        e_P_Q54 = pairing(P, Q54)
        e_P_Q23 = pairing(P, Q23)
        e_P_Q31 = pairing(P, Q31)

        @test e_P_Q54 == e_P_Q23 * e_P_Q31
    end
end
