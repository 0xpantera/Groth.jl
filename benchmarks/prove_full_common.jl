using Random
using GrothAlgebra
using GrothCurves
using GrothProofs

const PROVE_FULL_PHASE_ORDER = [
    "end_to_end",
    "witness_to_scalars",
    "msm_a_g1",
    "msm_b_g1",
    "msm_b_g2",
    "compute_h_total",
    "h_poly_assembly",
    "h_dense_quotient",
    "h_coset_fft",
    "h_parity_assert",
    "h_msm",
    "l_msm",
    "final_c",
]

field_to_bigint(x) = convert(BigInt, x)

function public_inputs_for(r1cs::R1CS, witness::Witness)
    return r1cs.num_public > 1 ? witness.values[2:r1cs.num_public] : eltype(witness.values)[]
end

function sample_benchmark_field(::Type{F}, rng::AbstractRNG) where F
    F === BN254Fr || throw(ArgumentError("Benchmark fixtures only support BN254Fr today"))
    return F(rand(rng, 0:GrothCurves.BN254_ORDER_R-1))
end

function generate_benchmark_r1cs(rng::AbstractRNG; num_constraints::Int=24, max_public::Int=8, value_range=-20:20, retries::Int=100)
    num_constraints >= 2 || error("expect at least two constraints")
    for _ in 1:retries
        num_vars = 1 + 1 + num_constraints + 2
        num_public = min(max_public, num_vars)

        F = BN254Fr
        L = zeros(F, num_constraints, num_vars)
        R = zeros(F, num_constraints, num_vars)
        O = zeros(F, num_constraints, num_vars)

        values = zeros(Int, num_vars)
        values[1] = 1

        idx_r = 2
        base_values = [rand(rng, value_range) for _ in 1:num_vars]
        for i in 2:num_public
            values[i] = base_values[i]
        end

        next_var = num_public + 1
        available = collect(3:num_public)
        isempty(available) && continue

        for c in 1:num_constraints
            if c == num_constraints || next_var > num_vars
                a_idx = rand(rng, available)
                b_idx = rand(rng, available)
                L[c, a_idx] = one(F)
                R[c, b_idx] = one(F)
                O[c, idx_r] = one(F)
                values[idx_r] = values[a_idx] * values[b_idx]
                continue
            end

            gate_type = rand(rng, (:mul, :affine))
            if gate_type == :mul
                a_idx = rand(rng, available)
                b_idx = rand(rng, available)
                L[c, a_idx] = one(F)
                R[c, b_idx] = one(F)
                O[c, next_var] = one(F)
                values[next_var] = values[a_idx] * values[b_idx]
            else
                a_idx = rand(rng, available)
                b_idx = rand(rng, available)
                α = rand(rng, value_range)
                β = rand(rng, value_range)
                γ = rand(rng, value_range)
                if α == 0 && β == 0 && γ == 0
                    γ = 1
                end
                L[c, a_idx] = bn254_fr(α)
                L[c, b_idx] += bn254_fr(β)
                L[c, 1] += bn254_fr(γ)
                R[c, 1] = one(F)
                O[c, next_var] = one(F)
                values[next_var] = α * values[a_idx] + β * values[b_idx] + γ
            end

            push!(available, next_var)
            next_var += 1
        end

        witness_vals = Vector{BN254Fr}(undef, num_vars)
        for i in 1:num_vars
            witness_vals[i] = bn254_fr(values[i])
        end

        r1cs = R1CS{BN254Fr}(num_vars, num_constraints, num_public, L, R, O)
        witness = Witness{BN254Fr}(witness_vals)
        is_satisfied(r1cs, witness) || continue
        return r1cs, witness
    end

    error("failed to generate satisfying benchmark R1CS")
end

function build_fixture(name::String, description::String, r1cs::R1CS, witness::Witness; setup_seed::Int, prove_seed::Int)
    qap = r1cs_to_qap(r1cs)
    keypair = setup_full(qap; rng=MersenneTwister(setup_seed))
    public_inputs = public_inputs_for(r1cs, witness)
    return (
        name = name,
        description = description,
        r1cs = r1cs,
        qap = qap,
        witness = witness,
        public_inputs = public_inputs,
        keypair = keypair,
        setup_seed = setup_seed,
        prove_seed = prove_seed,
    )
end

function default_prove_full_fixtures()
    continuity_r1cs = create_r1cs_example_sum_of_products()
    continuity_witness = create_witness_sum_of_products(3, 5, 7, 11)

    rng = MersenneTwister(20260331)
    generated_r1cs, generated_witness = generate_benchmark_r1cs(rng; num_constraints=24, max_public=8)

    return [
        build_fixture(
            "sum_of_products_small",
            "Historical continuity fixture used by the Groth16 end-to-end benchmark",
            continuity_r1cs,
            continuity_witness;
            setup_seed = 42,
            prove_seed = 1337,
        ),
        build_fixture(
            "generated_24_constraints",
            "Deterministic generated fixture for prove_full hotspot study",
            generated_r1cs,
            generated_witness;
            setup_seed = 20260332,
            prove_seed = 20260333,
        ),
    ]
end

function fixture_metadata(fixture)
    qap = fixture.qap
    return Dict{String,Any}(
        "name" => fixture.name,
        "description" => fixture.description,
        "num_constraints" => fixture.r1cs.num_constraints,
        "num_vars" => fixture.r1cs.num_vars,
        "num_public" => fixture.r1cs.num_public,
        "domain_size" => qap.domain.size,
        "domain_log_size" => qap.domain.log_size,
        "coset_domain_size" => qap.coset_domain.size,
        "setup_seed" => fixture.setup_seed,
        "prove_seed" => fixture.prove_seed,
    )
end

function witness_to_scalars(witness::Witness)
    return copy(witness.values)
end

function witness_to_scalars_bigint(witness::Witness)
    w_vals = witness.values
    scalars = Vector{BigInt}(undef, length(w_vals))
    @inbounds for i in eachindex(w_vals)
        scalars[i] = field_to_bigint(w_vals[i])
    end
    return scalars
end

function build_combined_polynomials(qap::QAP{F}, witness::Witness{F}) where F
    u_poly = zero(qap.u[1])
    v_poly = zero(qap.v[1])
    w_poly = zero(qap.w[1])
    @inbounds for i in 1:qap.num_vars
        wi = witness.values[i]
        u_poly = u_poly + wi * qap.u[i]
        v_poly = v_poly + wi * qap.v[i]
        w_poly = w_poly + wi * qap.w[i]
    end
    return u_poly, v_poly, w_poly
end

function compute_h_dense_path(qap::QAP{F}, u_poly, v_poly, w_poly) where F
    p_poly = fft_polynomial_multiply(u_poly, v_poly) - w_poly
    return polynomial_division(p_poly, qap.t)
end

function compute_h_coset_path(qap::QAP{F}, u_poly, v_poly, w_poly, dense_len::Int) where F
    coset = qap.coset_domain
    u_eval = fft(u_poly.coeffs, coset)
    v_eval = fft(v_poly.coeffs, coset)
    w_eval = fft(w_poly.coeffs, coset)
    h_eval = similar(u_eval)
    vanishing_inv = qap.vanishing_coset_inv
    @inbounds for i in eachindex(u_eval)
        numerator = u_eval[i] * v_eval[i] - w_eval[i]
        h_eval[i] = numerator * vanishing_inv[i]
    end
    coeffs = GrothAlgebra.ifft!(h_eval, coset)
    dense_len < length(coeffs) && (coeffs = coeffs[1:dense_len])
    return Polynomial{F}(coeffs)
end

function assert_h_parity(dense_result, coset_result)
    coset_result == dense_result || error("Coset Groth16 path diverged from dense fallback")
    return coset_result
end

function sample_prove_randomizers(fixture)
    F = eltype(fixture.witness.values)
    rng = MersenneTwister(fixture.prove_seed)
    return sample_benchmark_field(F, rng), sample_benchmark_field(F, rng)
end

function prove_query_accumulators(pk::ProvingKey, scalars::AbstractVector)
    return (
        A_acc_g1 = GrothAlgebra.multi_scalar_mul(pk.A_query_g1, scalars),
        B_acc_g1 = GrothAlgebra.multi_scalar_mul(pk.B_query_g1, scalars),
        B_acc_g2 = GrothAlgebra.multi_scalar_mul(pk.B_query_g2, scalars),
    )
end

function assemble_ab_terms(pk::ProvingKey, A_acc_g1, B_acc_g1, B_acc_g2, r, s)
    A1_g1 = pk.alpha_g1 + A_acc_g1
    B1_g1 = pk.beta_g1 + B_acc_g1
    A = A1_g1 + scalar_mul(pk.delta_g1, r)
    B = pk.beta_g2 + B_acc_g2 + scalar_mul(pk.delta_g2, s)
    return (A = A, B = B, A1_g1 = A1_g1, B1_g1 = B1_g1)
end

function h_msm(pk::ProvingKey, h_poly)
    hk = length(h_poly.coeffs)
    pts_h = pk.H_query_g1[1:hk]
    return GrothAlgebra.multi_scalar_mul(pts_h, h_poly.coeffs)
end

function l_msm(pk::ProvingKey, witness::Witness)
    m = length(witness.values)
    if m <= pk.num_public
        return zero(G1Point)
    end
    priv_scalars = witness.values[(pk.num_public+1):m]
    return GrothAlgebra.multi_scalar_mul(pk.L_query_g1, priv_scalars)
end

function assemble_c(pk::ProvingKey, A, B1_g1, H, L, r, s)
    rs_delta = scalar_mul(pk.delta_g1, r * s)
    g1_b_full = B1_g1 + scalar_mul(pk.delta_g1, s)
    return H + L + scalar_mul(g1_b_full, r) + scalar_mul(A, s) - rs_delta
end
