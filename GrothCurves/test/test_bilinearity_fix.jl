# Test bilinearity of BN254 pairing implementation after Frobenius fix
using GrothCurves
using GrothAlgebra
using Test

@testset "BN254 Pairing Bilinearity" begin
    # Get generators
    P = g1_generator()
    Q = g2_generator()

    @testset "Basic bilinearity" begin
        # Test 1: e([2]P, Q) = e(P, [2]Q)
        P2 = P + P
        Q2 = Q + Q
        pairing1 = pairing(P2, Q)
        pairing2 = pairing(P, Q2)
        @test pairing1 == pairing2

        # Test 2: e([2]P, Q) = e(P, Q)^2
        e_P_Q = pairing(P, Q)
        e_P_Q_squared = e_P_Q^2
        e_2P_Q = pairing(P2, Q)
        @test e_2P_Q == e_P_Q_squared

        # Test 3: e(P, [2]Q) = e(P, Q)^2
        e_P_2Q = pairing(P, Q2)
        @test e_P_2Q == e_P_Q_squared

        # Test 4: e([3]P, Q) = e(P, [3]Q)
        P3 = P + P + P
        Q3 = Q + Q + Q
        e_3P_Q = pairing(P3, Q)
        e_P_3Q = pairing(P, Q3)
        @test e_3P_Q == e_P_3Q

        # Test 5: e([3]P, Q) = e(P, Q)^3
        e_P_Q_cubed = e_P_Q^3
        @test e_3P_Q == e_P_Q_cubed
    end

    @testset "Additive bilinearity" begin
        # Test 6: e(P1 + P2, Q) = e(P1, Q) * e(P2, Q)
        # Use P as P1 and [2]P as P2, so P1 + P2 = [3]P
        P2 = P + P
        P3 = P + P + P
        e_P1_plus_P2_Q = pairing(P3, Q)
        e_P1_Q = pairing(P, Q)
        e_P2_Q = pairing(P2, Q)
        product = e_P1_Q * e_P2_Q
        @test e_P1_plus_P2_Q == product

        # Test 7: e(P, Q1 + Q2) = e(P, Q1) * e(P, Q2)
        # Use Q as Q1 and [2]Q as Q2, so Q1 + Q2 = [3]Q
        Q2 = Q + Q
        Q3 = Q + Q + Q
        e_P_Q1_plus_Q2 = pairing(P, Q3)
        e_P_Q1 = pairing(P, Q)
        e_P_Q2 = pairing(P, Q2)
        product2 = e_P_Q1 * e_P_Q2
        @test e_P_Q1_plus_Q2 == product2
    end

    @testset "Non-degeneracy and identity" begin
        # Test 8: Non-degeneracy - e(P, Q) ≠ 1
        e_P_Q = pairing(P, Q)
        @test e_P_Q != one(Fp12Element)

        # Test 9: e(O, Q) = 1 (where O is point at infinity)
        O = zero(G1Point)
        e_O_Q = pairing(O, Q)
        @test e_O_Q == one(Fp12Element)

        # Test 10: e(P, O) = 1 (where O is point at infinity)
        O_G2 = zero(G2Point)
        e_P_O = pairing(P, O_G2)
        @test e_P_O == one(Fp12Element)
    end

    @testset "Scalar multiplication consistency" begin
        # Test that e([n]P, Q) = e(P, [n]Q) = e(P, Q)^n for small n
        for n in [4, 5, 6, 7]
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
            e_P_Q = pairing(P, Q)
            e_P_Q_n = e_P_Q^n

            @test e_nP_Q == e_P_nQ
            @test e_nP_Q == e_P_Q_n
        end
    end

    @testset "Negation properties" begin
        # Test that e(-P, Q) = e(P, -Q) = e(P, Q)^(-1)
        neg_P = -P
        neg_Q = -Q

        e_P_Q = pairing(P, Q)
        e_neg_P_Q = pairing(neg_P, Q)
        e_P_neg_Q = pairing(P, neg_Q)
        e_P_Q_inv = inv(e_P_Q)

        @test e_neg_P_Q == e_P_Q_inv
        @test e_P_neg_Q == e_P_Q_inv
        @test e_neg_P_Q == e_P_neg_Q
    end
end
