using GrothProofs
using GrothAlgebra: bn254_scalar, evaluate
using Test
using Random

include("random_circuits.jl")

function public_inputs_for(r1cs::R1CS, witness::Witness)
    return r1cs.num_public > 1 ? witness.values[2:r1cs.num_public] : eltype(witness.values)[]
end

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

    # Verify with correct public inputs (arkworks-style, excludes leading 1)
    public_inputs = public_inputs_for(r1cs, witness)
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
        public_inputs = public_inputs_for(r1cs, witness)

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
    public_inputs = public_inputs_for(r1cs, witness)
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
        public_inputs = public_inputs_for(r1cs, witness)
        proof = prove_full(keypair.pk, qap, witness; rng=MersenneTwister(x + y + z + u))
        @test verify_full(keypair.vk, proof, public_inputs)

        # Prepared verifier path should agree
        pvk = prepare_verifying_key(keypair.vk)
        pin = prepare_inputs(pvk, public_inputs)
        @test verify_with_prepared(pvk, proof, pin)
    end

    # Rerandomization invariance: different seeds produce valid proofs
    x, y, z, u = 3, 4, 5, 6
    witness = create_witness_sum_of_products(x, y, z, u)
    public_inputs = public_inputs_for(r1cs, witness)
    proof1 = prove_full(keypair.pk, qap, witness; rng=MersenneTwister(1))
    proof2 = prove_full(keypair.pk, qap, witness; rng=MersenneTwister(2))
    @test verify_full(keypair.vk, proof1, public_inputs)
    @test verify_full(keypair.vk, proof2, public_inputs)
    pvk = prepare_verifying_key(keypair.vk)
    pin = prepare_inputs(pvk, public_inputs)
    @test verify_with_prepared(pvk, proof1, pin)
    @test verify_with_prepared(pvk, proof2, pin)
end

@testset "Groth16 affine product circuit" begin
    r1cs = create_r1cs_example_affine_product()
    qap = r1cs_to_qap(r1cs)
    keypair = setup_full(qap; rng=MersenneTwister(5151))

    inputs_list = [
        (0, 0, 0, 0),
        (1, 2, 3, 4),
        (7, 11, 13, 17),
        (25, 30, 5, 9),
    ]

    for (x, y, z, u) in inputs_list
        witness = create_witness_affine_product(x, y, z, u)
        @test is_satisfied(r1cs, witness)
        public_inputs = public_inputs_for(r1cs, witness)
        proof = prove_full(keypair.pk, qap, witness; rng=MersenneTwister(x + y + z + u + 1))
        @test verify_full(keypair.vk, proof, public_inputs)

        pvk = prepare_verifying_key(keypair.vk)
        pin = prepare_inputs(pvk, public_inputs)
        @test verify_with_prepared(pvk, proof, pin)
    end

    # Tamper helper witness entry to ensure detection
    witness = create_witness_affine_product(2, 3, 5, 7)
    bad = copy(witness.values)
    bad[7] += one(eltype(bad))  # corrupt s1
    bad_witness = Witness(bad)
    @test !is_satisfied(r1cs, bad_witness)
end

@testset "Groth16 square-offset circuit" begin
    r1cs = create_r1cs_example_square_offset()
    qap = r1cs_to_qap(r1cs)
    keypair = setup_full(qap; rng=MersenneTwister(8181))

    inputs_list = [
        (0, 0, 0),
        (2, 3, 5),
        (9, 1, 4),
        (123, 456, 789),
    ]

    for (x, y, c) in inputs_list
        witness = create_witness_square_offset(x, y, c)
        @test is_satisfied(r1cs, witness)
        public_inputs = public_inputs_for(r1cs, witness)
        proof = prove_full(keypair.pk, qap, witness; rng=MersenneTwister(x + y + c + 5))
        @test verify_full(keypair.vk, proof, public_inputs)

        pvk = prepare_verifying_key(keypair.vk)
        pin = prepare_inputs(pvk, public_inputs)
        @test verify_with_prepared(pvk, proof, pin)
    end

    # Wrong public input (offset) should fail verification
    witness = create_witness_square_offset(4, 6, 8)
    public_inputs = public_inputs_for(r1cs, witness)
    proof = prove_full(keypair.pk, qap, witness; rng=MersenneTwister(4040))
    bad_inputs = copy(public_inputs)
    bad_inputs[end] += one(eltype(bad_inputs))
    @test !verify_full(keypair.vk, proof, bad_inputs)
end

@testset "Groth16 randomized circuits" begin
    seeds = (101, 202, 303, 404, 505)
    for seed in seeds
        rng = MersenneTwister(seed)
        circuit = generate_small_r1cs(rng)
        r1cs = circuit.r1cs
        witness = circuit.witness
        @test is_satisfied(r1cs, witness)

        qap = r1cs_to_qap(r1cs)
        keypair = setup_full(qap; rng=MersenneTwister(seed + 1))
        proof = prove_full(keypair.pk, qap, witness; rng=MersenneTwister(seed + 2))

        public_inputs = public_inputs_for(r1cs, witness)
        @test verify_full(keypair.vk, proof, public_inputs)

        pvk = prepare_verifying_key(keypair.vk)
        pin = prepare_inputs(pvk, public_inputs)
        @test verify_with_prepared(pvk, proof, pin)

        # Mutate a random public slot and expect failure
        touched = circuit.used_public
        if !isempty(touched)
            bad_inputs = copy(public_inputs)
            idx = rand(rng, touched)
            bad_inputs[idx - 1] += one(eltype(bad_inputs))
            @test !verify_full(keypair.vk, proof, bad_inputs)
        end

        h_dense = compute_h_polynomial(qap, witness; use_coset=false)
        h_coset = compute_h_polynomial(qap, witness; use_coset=true)
        @test h_coset == h_dense
    end
end

@testset "Groth16 prepared path negatives" begin
    r1cs = create_r1cs_example_sum_of_products()
    qap = r1cs_to_qap(r1cs)
    keypair = setup_full(qap; rng=MersenneTwister(5150))

    witness = create_witness_sum_of_products(4, 6, 8, 10)
    public_inputs = public_inputs_for(r1cs, witness)
    proof = prove_full(keypair.pk, qap, witness; rng=MersenneTwister(42))

    pvk = prepare_verifying_key(keypair.vk)
    prepared = prepare_inputs(pvk, public_inputs)
    @test verify_with_prepared(pvk, proof, prepared)

    # prepare_inputs should reject empty or overlong vectors
    @test_throws ArgumentError prepare_inputs(pvk, eltype(public_inputs)[])
    overlong = vcat(public_inputs, public_inputs[end])
    @test_throws ArgumentError prepare_inputs(pvk, overlong)

    # Mutated public inputs produce distinct prepared accumulator and must fail
    bad_inputs = copy(public_inputs)
    bad_inputs[3] += one(eltype(bad_inputs))
    bad_prepared = prepare_inputs(pvk, bad_inputs)
    @test !verify_with_prepared(pvk, proof, bad_prepared)

    # Tampering proof components is detected by prepared verifier as well
    tweaked_proof = Groth16Proof(proof.A + keypair.pk.alpha_g1, proof.B, proof.C)
    @test !verify_with_prepared(pvk, tweaked_proof, prepared)
end

@testset "QAP divisibility spot-checks" begin
    # Use both circuits to check U*V - W equals h*t at a few random points
    for builder in (
        create_r1cs_example_multiplication,
        create_r1cs_example_sum_of_products,
        create_r1cs_example_affine_product,
        create_r1cs_example_square_offset,
    )
        r1cs = builder()
        qap = r1cs_to_qap(r1cs)

        # Simple witness
        w = builder === create_r1cs_example_multiplication ?
            create_witness_multiplication(3, 5, 7, 11) :
            builder === create_r1cs_example_sum_of_products ?
                create_witness_sum_of_products(3, 5, 7, 11) :
                builder === create_r1cs_example_affine_product ?
                    create_witness_affine_product(2, 3, 5, 7) :
                    create_witness_square_offset(2, 3, 5)

        # Build combined polynomials
        u_poly = zero(qap.u[1])
        v_poly = zero(qap.v[1])
        w_poly = zero(qap.w[1])
        for i in 1:qap.num_vars
            u_poly = u_poly + w.values[i] * qap.u[i]
            v_poly = v_poly + w.values[i] * qap.v[i]
            w_poly = w_poly + w.values[i] * qap.w[i]
        end
        h_dense = compute_h_polynomial(qap, w; use_coset=false)
        h_coset = compute_h_polynomial(qap, w; use_coset=true)
        @test h_coset == h_dense
        h_poly = h_coset
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

@testset "QAP coset vanishing cases" begin
    r1cs_power = create_r1cs_example_square_offset()
    qap_power = r1cs_to_qap(r1cs_power)
    @test qap_power.num_constraints == qap_power.domain.size
    witness_power = create_witness_square_offset(2, 3, 5)
    h_dense_power = compute_h_polynomial(qap_power, witness_power; use_coset=false)
    h_coset_power = compute_h_polynomial(qap_power, witness_power; use_coset=true)
    @test h_coset_power == h_dense_power

    r1cs_subset = create_r1cs_example_multiplication()
    qap_subset = r1cs_to_qap(r1cs_subset)
    @test qap_subset.num_constraints < qap_subset.domain.size
    witness_subset = create_witness_multiplication(3, 5, 7, 11)
    h_dense_subset = compute_h_polynomial(qap_subset, witness_subset; use_coset=false)
    h_coset_subset = compute_h_polynomial(qap_subset, witness_subset; use_coset=true)
    @test h_coset_subset == h_dense_subset
end

@testset "Groth16 rejects unsatisfied witness" begin
    r1cs = create_r1cs_example_multiplication()
    qap = r1cs_to_qap(r1cs)
    keypair = setup_full(qap; rng=MersenneTwister(6060))

    witness = create_witness_multiplication(3, 5, 7, 11)
    bad_values = copy(witness.values)
    bad_values[2] += one(eltype(bad_values))  # perturb expected output r
    bad_witness = Witness(bad_values)
    @test !is_satisfied(r1cs, bad_witness)

    @test_throws ErrorException prove_full(keypair.pk, qap, bad_witness; rng=MersenneTwister(1234))
end
