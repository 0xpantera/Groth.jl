using GrothProofs
using GrothAlgebra: bn254_scalar, evaluate
using Test
using Random

@testset "Groth16 full pipeline" begin
    # Build example R1CS and witness for r = x*y*z*u
    r1cs = create_r1cs_example_multiplication()
    x, y, z, u = 3, 5, 7, 11
    witness = create_witness_multiplication(x, y, z, u)
    @test is_satisfied(r1cs, witness)

    # Convert to QAP
    qap = r1cs_to_qap(r1cs)

    # Setup keys
    keypair = setup_full(qap; rng=MersenneTwister(42))

    # Prove with randomized r,s (seeded for reproducibility)
    proof = prove_full(keypair.pk, qap, witness; rng=MersenneTwister(1337))

    # Verify with correct public inputs
    # Arkworks prepare_inputs uses gamma_abc[1..] with inputs excluding the leading 1.
    # Our verify_full expects the leading 1 included at index 1, then the rest inputs.
    public_inputs = witness.values[1:r1cs.num_public]
    @test verify_full(keypair.vk, proof, public_inputs)

    # Negative test: mutate a public input -> should fail
    bad_inputs = copy(public_inputs)
    bad_inputs[2] = bad_inputs[2] + one(typeof(bad_inputs[2]))
    @test !verify_full(keypair.vk, proof, bad_inputs)

    # Negative test: tweak proof element -> should fail
    tweak = keypair.pk.alpha_g1
    bad_proof = Groth16Proof(proof.A + tweak, proof.B, proof.C)
    @test !verify_full(keypair.vk, bad_proof, public_inputs)
end

@testset "Groth16 randomness and multi-inputs" begin
    # Use multiple input sets and random r,s seeds
    r1cs = create_r1cs_example_multiplication()
    qap = r1cs_to_qap(r1cs)

    # Different inputs (x,y,z,u) to exercise IC composition
    inputs_list = [(2, 3, 5, 7), (9, 1, 4, 6), (13, 17, 19, 23)]

    # Deterministic keypair for repeatable tests
    keypair = setup_full(qap; rng=MersenneTwister(2024))

    for (idx, (x, y, z, u)) in enumerate(inputs_list)
        witness = create_witness_multiplication(x, y, z, u)
        @test is_satisfied(r1cs, witness)
        public_inputs = witness.values[1:r1cs.num_public]

        # Randomized prover seeds (non-zero randomness)
        rng = MersenneTwister(10_000 + idx)
        proof = prove_full(keypair.pk, qap, witness; rng=rng, debug_no_random=false)
        @test verify_full(keypair.vk, proof, public_inputs)
    end
end

@testset "Groth16 additional negatives" begin
    r1cs = create_r1cs_example_multiplication()
    qap = r1cs_to_qap(r1cs)
    keypair = setup_full(qap; rng=MersenneTwister(7))

    x, y, z, u = 4, 8, 3, 2
    witness = create_witness_multiplication(x, y, z, u)
    public_inputs = witness.values[1:r1cs.num_public]
    proof = prove_full(keypair.pk, qap, witness; rng=MersenneTwister(99))

    # Wrong public input length (drop last)
    @test !verify_full(keypair.vk, proof, public_inputs[1:end-1])
    # Wrong public input length (add extra)
    extra = vcat(public_inputs, public_inputs[end])
    @test !verify_full(keypair.vk, proof, extra)
    # Swap two public inputs
    swapped = copy(public_inputs)
    if length(swapped) >= 4
        tmp = swapped[3]
        swapped[3] = swapped[4]
        swapped[4] = tmp
        @test !verify_full(keypair.vk, proof, swapped)
    end
    # Tamper B and C separately
    bad_proof_B = Groth16Proof(proof.A, keypair.pk.beta_g2, proof.C)
    @test !verify_full(keypair.vk, bad_proof_B, public_inputs)
    bad_proof_C = Groth16Proof(proof.A, proof.B, keypair.pk.delta_g1)
    @test !verify_full(keypair.vk, bad_proof_C, public_inputs)
end

@testset "Groth16 sum-of-products circuit" begin
    # r = x*y + z*u
    r1cs = create_r1cs_example_sum_of_products()
    qap = r1cs_to_qap(r1cs)

    keypair = setup_full(qap; rng=MersenneTwister(4242))

    # Exercise multiple input patterns (zeros/ones/large)
    inputs_list = [
        (0, 5, 7, 11),
        (1, 1, 1, 1),
        (2, 0, 3, 0),
        (12345, 6789, 101112, 131415),
    ]

    for (x, y, z, u) in inputs_list
        witness = create_witness_sum_of_products(x, y, z, u)
        @test is_satisfied(r1cs, witness)
        public_inputs = witness.values[1:r1cs.num_public]
        proof = prove_full(keypair.pk, qap, witness; rng=MersenneTwister(x + y + z + u))
        @test verify_full(keypair.vk, proof, public_inputs)
    end

    # Rerandomization invariance: different seeds produce valid proofs
    x, y, z, u = 3, 4, 5, 6
    witness = create_witness_sum_of_products(x, y, z, u)
    public_inputs = witness.values[1:r1cs.num_public]
    proof1 = prove_full(keypair.pk, qap, witness; rng=MersenneTwister(1))
    proof2 = prove_full(keypair.pk, qap, witness; rng=MersenneTwister(2))
    @test verify_full(keypair.vk, proof1, public_inputs)
    @test verify_full(keypair.vk, proof2, public_inputs)
end

@testset "QAP divisibility spot-checks" begin
    # Use both circuits to check U*V - W equals h*t at a few random points
    for builder in (create_r1cs_example_multiplication, create_r1cs_example_sum_of_products)
        r1cs = builder()
        qap = r1cs_to_qap(r1cs)

        # Simple witness
        w = builder === create_r1cs_example_multiplication ?
            create_witness_multiplication(3, 5, 7, 11) :
            create_witness_sum_of_products(3, 5, 7, 11)

        # Build combined polynomials
        u_poly = zero(qap.u[1])
        v_poly = zero(qap.v[1])
        w_poly = zero(qap.w[1])
        for i in 1:qap.num_vars
            u_poly = u_poly + w.values[i] * qap.u[i]
            v_poly = v_poly + w.values[i] * qap.v[i]
            w_poly = w_poly + w.values[i] * qap.w[i]
        end
        h_poly = compute_h_polynomial(qap, w)
        t_poly = qap.t

        # Choose a few points outside the domain [1..n]
        pts = [bn254_scalar(qap.num_constraints + k) for k in (1, 2, 3)]
        for x in pts
            lhs = evaluate(u_poly, x) * evaluate(v_poly, x) - evaluate(w_poly, x)
            rhs = evaluate(h_poly, x) * evaluate(t_poly, x)
            if lhs != rhs
                println("[QAP check] builder=$(builder) x=$(Integer(x.value)) lhs=$(Integer(lhs.value)) rhs=$(Integer(rhs.value))")
            end
            @test lhs == rhs
        end
    end
end
