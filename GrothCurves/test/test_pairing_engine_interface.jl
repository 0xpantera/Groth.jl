using Test
using Random

function test_pairing_engine_interface(engine::AbstractPairingEngine)
    @testset "Pairing engine interface" begin
        # Acquire default generators
        P = g1_generator()
        Q = g2_generator()

        # Basic zero handling
        @test pairing(engine, zero(G1Point), Q) == one(GTElement)
        @test pairing(engine, P, zero(G2Point)) == one(GTElement)
        @test pairing(engine, zero(G1Point), zero(G2Point)) == one(GTElement)

        # Bilinearity across random scalars
        rng = MersenneTwister(2025)
        for _ in 1:5
            a = rand(rng, 1:50)
            b = rand(rng, 1:50)
            Pa = scalar_mul(P, a)
            Qb = scalar_mul(Q, b)
            lhs = pairing(engine, Pa, Qb)
            rhs = pairing(engine, P, Q)^(a*b)
            @test lhs == rhs
        end

        # Batch pairing equals product of individual pairings
        points_g1 = [P, scalar_mul(P, 2), scalar_mul(P, 3)]
        points_g2 = [Q, scalar_mul(Q, 2), scalar_mul(Q, 3)]
        prod_separate = pairing(engine, points_g1[1], points_g2[1]) *
                        pairing(engine, points_g1[2], points_g2[2]) *
                        pairing(engine, points_g1[3], points_g2[3])
        prod_batch = pairing_batch(engine, points_g1, points_g2)
        @test prod_batch == prod_separate
    end
end

# Execute the interface test with the default BN254 engine
@testset "BN254 pairing engine interface" begin
    test_pairing_engine_interface(BN254_ENGINE)
end
