"""
Groth16 proof system implementation (production-style scaffold).

Provides Proving/Verification keys and standard Groth16 equations.
The prover uses the coset FFT pipeline by default; the dense path survives as a
consistency check.
"""

using GrothAlgebra
using GrothCurves
using Random

# R1CS and QAP are included at the package level (GrothProofs.jl)

@inline function _rand_field(::Type{F}, rng::AbstractRNG) where F
    F === BN254Fr || throw(ArgumentError("No CSPRNG-ready sampler configured for field $F"))
    return F(rand(rng, 0:GrothCurves.BN254_ORDER_R-1))
end

@inline function _rand_field_nonzero(::Type{F}, rng::AbstractRNG) where F
    F === BN254Fr || throw(ArgumentError("No CSPRNG-ready sampler configured for field $F"))
    return F(rand(rng, 1:GrothCurves.BN254_ORDER_R-1))
end

"""
    Groth16Proof

Store the three group elements produced by Groth16.
"""
struct Groth16Proof
    A::G1Point  # [A]₁
    B::G2Point  # [B]₂
    C::G1Point  # [C]₁
end

"""
    ProvingKey

Store the precomputed group elements and metadata required by the Groth16 prover.

This follows the conventional separation:
- ProvingKey contains query vectors and delta terms for prover computation
- VerificationKey contains minimal elements for the verifier equation
"""
struct ProvingKey
    # Fixed elements
    alpha_g1::G1Point
    beta_g1::G1Point
    beta_g2::G2Point
    delta_g1::G1Point
    delta_g2::G2Point

    # Query vectors (per-variable)
    A_query_g1::Vector{G1Point}   # [u_i(τ)]₁
    B_query_g1::Vector{G1Point}   # [v_i(τ)]₁ (for cross-term in C)
    B_query_g2::Vector{G2Point}   # [v_i(τ)]₂
    C_query_g1::Vector{G1Point}   # [w_i(τ)]₁

    # h and l queries
    H_query_g1::Vector{G1Point}   # [τ^k · t(τ) / δ]₁ for k = 0..(deg(h))
    L_query_g1::Vector{G1Point}   # [(β·u_i(τ) + α·v_i(τ) + w_i(τ)) / δ]₁ for private i

    # Metadata
    num_public::Int               # includes constant 1
end

"""
    VerificationKey

Store the elements required by the Groth16 verifier, including the prepared
engine handle.
"""
struct VerificationKey{E<:AbstractPairingEngine}
    alpha_g1::G1Point
    beta_g2::G2Point
    gamma_g2::G2Point
    delta_g2::G2Point
    IC::Vector{G1Point}  # [ (β·u_i(τ) + α·v_i(τ) + w_i(τ)) / γ ]₁ for public i (including 1)
    engine::E
end

"""
    Keypair

Bundle the proving and verification keys emitted by `setup_full`.
"""
struct Keypair{E<:AbstractPairingEngine}
    pk::ProvingKey
    vk::VerificationKey{E}
end

# Export types
export Groth16Proof

"""
    setup_full(qap::QAP{F}; rng=Random.GLOBAL_RNG, engine=BN254_ENGINE) where F

Generate production-style Groth16 keys (ProvingKey, VerificationKey).
Keys work with the coset FFT prover pipeline while retaining the dense fallback
for parity assertions.

Security note: `rng` defaults to `Random.GLOBAL_RNG` (often `MersenneTwister`).
Use a CSPRNG for real deployments.
"""
function setup_full(qap::QAP{F}; rng::AbstractRNG=Random.GLOBAL_RNG, engine::AbstractPairingEngine=BN254_ENGINE) where F
    E = typeof(engine)
    # Sample toxic waste
    α = _rand_field(F, rng)
    β = _rand_field(F, rng)
    γ = _rand_field_nonzero(F, rng)
    δ = _rand_field_nonzero(F, rng)
    τ = _rand_field(F, rng)

    # Generators
    g1 = g1_generator()
    g2 = g2_generator()

    # Precompute inverses
    γ_inv = inv(γ)
    δ_inv = inv(δ)

    # Evaluate QAP polynomials at τ
    m = qap.num_vars
    A_eval = Vector{F}(undef, m)
    B_eval = Vector{F}(undef, m)
    C_eval = Vector{F}(undef, m)
    for i in 1:m
        A_eval[i] = evaluate(qap.u[i], τ)
        B_eval[i] = evaluate(qap.v[i], τ)
        C_eval[i] = evaluate(qap.w[i], τ)
    end

    # Fixed-base tables for batch generation
    g1tab = GrothAlgebra.build_fixed_table(g1)
    g2tab = GrothAlgebra.build_fixed_table(g2)

    # Map evaluations into group queries using fixed-base batch mul
    A_query_g1 = GrothAlgebra.batch_mul(g1tab, A_eval)
    B_query_g2 = GrothAlgebra.batch_mul(g2tab, B_eval)
    B_query_g1 = GrothAlgebra.batch_mul(g1tab, B_eval)
    C_query_g1 = GrothAlgebra.batch_mul(g1tab, C_eval)

    # Target polynomial at τ
    t_tau = evaluate(qap.t, τ)

    # τ^k powers for H query (match arkworks m_raw - 1 ≈ num_constraints + num_vars - 1)
    # This supports h(τ) expansion: Σ h_k τ^k
    max_k = qap.num_constraints + qap.num_vars - 1
    tau_powers = Vector{F}(undef, max_k)
    tau_powers[1] = one(F)
    for k in 2:max_k
        tau_powers[k] = tau_powers[k-1] * τ
    end

    # H_query_g1: [τ^k · t(τ) / δ]₁
    H_query_g1 = GrothAlgebra.batch_mul(g1tab, [t_tau * tau_powers[k] * δ_inv for k in 1:max_k])

    # L_query for private variables only: [(β·u_i(τ) + α·v_i(τ) + w_i(τ)) / δ]₁
    num_public = qap.num_public
    L_query_g1 = G1Point[]
    if m > num_public
        L_scalars = F[]
        sizehint!(L_scalars, m - num_public)
        for i in (num_public+1):m
            acc = β * A_eval[i] + α * B_eval[i] + C_eval[i]
            push!(L_scalars, acc * δ_inv)
        end
        L_query_g1 = GrothAlgebra.batch_mul(g1tab, L_scalars)
    end

    # IC for public inputs: [(β·u_i(τ) + α·v_i(τ) + w_i(τ)) / γ]₁, including 1
    IC_scalars = [(β * A_eval[i] + α * B_eval[i] + C_eval[i]) * γ_inv for i in 1:num_public]
    IC = GrothAlgebra.batch_mul(g1tab, IC_scalars)

    # Normalize queries to affine for efficient storage
    GrothCurves.batch_to_affine!(A_query_g1)
    GrothCurves.batch_to_affine!(B_query_g1)
    GrothCurves.batch_to_affine!(C_query_g1)
    GrothCurves.batch_to_affine!(H_query_g1)
    GrothCurves.batch_to_affine!(IC)
    GrothCurves.batch_to_affine!(L_query_g1)
    GrothCurves.batch_to_affine!(B_query_g2)

    # Fixed elements
    alpha_g1 = scalar_mul(g1, α)
    beta_g1 = scalar_mul(g1, β)
    beta_g2 = scalar_mul(g2, β)
    gamma_g2 = scalar_mul(g2, γ)
    delta_g1 = scalar_mul(g1, δ)
    delta_g2 = scalar_mul(g2, δ)

    pk = ProvingKey(
        alpha_g1,
        beta_g1,
        beta_g2,
        delta_g1,
        delta_g2,
        A_query_g1,
        B_query_g1,
        B_query_g2,
        C_query_g1,
        H_query_g1,
        L_query_g1,
        num_public,
    )

    vk = VerificationKey{E}(
        alpha_g1,
        beta_g2,
        gamma_g2,
        delta_g2,
        IC,
        engine,
    )

    return Keypair{E}(pk, vk)
end

"""
    prove_full(pk::ProvingKey, qap::QAP{F}, witness::Witness{F}; rng=Random.GLOBAL_RNG, debug_no_random=false) where F

Construct a Groth16 proof with the coset FFT pipeline.

Set `debug_no_random=true` to omit prover randomness during testing.
Security note: `rng` defaults to `Random.GLOBAL_RNG` (often `MersenneTwister`).
Use a CSPRNG for real deployments.
"""
function prove_full(pk::ProvingKey, qap::QAP{F}, witness::Witness{F}; rng::AbstractRNG=Random.GLOBAL_RNG, debug_no_random::Bool=false) where F
    w_vals = witness.values
    m = qap.num_vars

    # Randomizers
    r = debug_no_random ? zero(F) : _rand_field(F, rng)
    s = debug_no_random ? zero(F) : _rand_field(F, rng)

    # Accumulators for queries via MSM (fallbacks to scalar loops for tiny sizes)
    scalars = w_vals[1:m]
    A_acc_g1 = GrothAlgebra.multi_scalar_mul(pk.A_query_g1, scalars)
    B_acc_g2 = GrothAlgebra.multi_scalar_mul(pk.B_query_g2, scalars)
    B_acc_g1 = GrothAlgebra.multi_scalar_mul(pk.B_query_g1, scalars)

    # A and B
    A1_g1 = pk.alpha_g1 + A_acc_g1
    B1_g1 = pk.beta_g1 + B_acc_g1
    A = A1_g1 + scalar_mul(pk.delta_g1, r)
    B = pk.beta_g2 + B_acc_g2 + scalar_mul(pk.delta_g2, s)

    # h(x): compute via coset FFT path and map coefficients via H_query
    h_poly = compute_h_polynomial(qap, witness)
    hk = length(h_poly.coeffs)
    pts_h = pk.H_query_g1[1:hk]
    scalars_h = h_poly.coeffs
    H = GrothAlgebra.multi_scalar_mul(pts_h, scalars_h)

    # L for private variables
    # L for private variables via MSM
    if m > pk.num_public
        priv_scalars = w_vals[(pk.num_public+1):m]
        L = GrothAlgebra.multi_scalar_mul(pk.L_query_g1, priv_scalars)
    else
        L = zero(G1Point)
    end

    # Cross terms in C (arkworks style): C = s*g_a + r*g1_b - r*s*δ + L + H
    # where g_a = A (includes r*δ), and g1_b = (β + Σ v_i) + s*δ in G1
    rs_delta = scalar_mul(pk.delta_g1, r * s)
    g1_b_full = B1_g1 + scalar_mul(pk.delta_g1, s)
    C = H + L + scalar_mul(g1_b_full, r) + scalar_mul(A, s) - rs_delta

    return Groth16Proof(A, B, C)
end

"""
    verify_full(vk::VerificationKey, proof::Groth16Proof, public_inputs::Vector{F}) where F

Verify a Groth16 proof with the production-style verification key.

Checks the standard pairing equation `e(A,B) == e(α,β) · e(vk_x,γ) · e(C,δ)` and
performs on-curve plus subgroup checks for `A`, `B`, and `C`.
"""
function verify_full(vk::VerificationKey{E}, proof::Groth16Proof, public_inputs::Vector{F}) where {E<:AbstractPairingEngine,F}
    # On-curve checks
    if !(GrothCurves.is_on_curve(proof.A) && GrothCurves.is_on_curve(proof.B) && GrothCurves.is_on_curve(proof.C))
        return false
    end

    # Subgroup checks: multiply by r and ensure infinity
    r = GrothCurves.BN254_ORDER_R
    if !(iszero(scalar_mul(proof.A, r)) && iszero(scalar_mul(proof.B, r)) && iszero(scalar_mul(proof.C, r)))
        return false
    end

    # Compute vk_x = IC[1] + Σ input[i] * IC[i+1]
    # public_inputs excludes the leading 1 (arkworks-style)
    if isempty(vk.IC)
        return false
    end
    expected_len = length(vk.IC) - 1
    if length(public_inputs) != expected_len
        return false
    end
    vk_x = vk.IC[1]
    if expected_len > 0
        pts_ic = vk.IC[2:end]
        scal_ic = public_inputs
        vk_x += GrothAlgebra.multi_scalar_mul(pts_ic, scal_ic)
    end

    engine = vk.engine
    lhs = pairing(engine, proof.A, proof.B)
    rhs = pairing(engine, vk.alpha_g1, vk.beta_g2) * pairing(engine, vk_x, vk.gamma_g2) * pairing(engine, proof.C, vk.delta_g2)
    return lhs == rhs
end

# Export production-style API
export ProvingKey, VerificationKey, Keypair
export setup_full, prove_full, verify_full

"""
    PreparedVerificationKey

Cache `e(α,β)` and negated `γ`, `δ` elements for the prepared verifier path.
"""
struct PreparedVerificationKey{E<:AbstractPairingEngine}
    vk::VerificationKey{E}
    alpha_g1_beta_g2::GTElement
    gamma_g2_neg::G2Point
    delta_g2_neg::G2Point
end

"""
    prepare_verifying_key(vk::VerificationKey)

Compute cached `e(α,β)` and negated `γ`, `δ` for prepared verification.
"""
function prepare_verifying_key(vk::VerificationKey{E}) where {E<:AbstractPairingEngine}
    engine = vk.engine
    return PreparedVerificationKey{E}(
        vk,
        pairing(engine, vk.alpha_g1, vk.beta_g2),
        -vk.gamma_g2,
        -vk.delta_g2,
    )
end

"""
    prepare_inputs(pvk::PreparedVerificationKey, public_inputs::Vector{F}) where F

Compute ``vk_x = IC[1] + \\sum_i input[i] \\cdot IC[i+1]``.

Public inputs exclude the leading `1` (arkworks-style).
"""
function prepare_inputs(pvk::PreparedVerificationKey{E}, public_inputs::Vector{F}) where {E<:AbstractPairingEngine,F}
    vk = pvk.vk
    if isempty(vk.IC)
        throw(ArgumentError("Verification key has empty IC vector"))
    end
    expected_len = length(vk.IC) - 1
    if length(public_inputs) != expected_len
        throw(ArgumentError("Expected $(expected_len) public inputs, got $(length(public_inputs))"))
    end
    acc = vk.IC[1]
    for i in 1:expected_len
        xi = public_inputs[i]
        iszero(xi) && continue
        acc += scalar_mul(vk.IC[i+1], xi)
    end
    return acc
end

"""
    verify_with_prepared(pvk::PreparedVerificationKey, proof::Groth16Proof, prepared_inputs::G1Point)

Verify a Groth16 proof using prepared verifier artifacts.

Compares the product pairing to the cached `e(α,β)` and reuses the on-curve and
subgroup checks from `verify_full`.
"""
function verify_with_prepared(pvk::PreparedVerificationKey{E}, proof::Groth16Proof, prepared_inputs::G1Point) where {E<:AbstractPairingEngine}
    # On-curve and subgroup checks
    if !(GrothCurves.is_on_curve(proof.A) && GrothCurves.is_on_curve(proof.B) && GrothCurves.is_on_curve(proof.C))
        return false
    end
    r = GrothCurves.BN254_ORDER_R
    if !(iszero(scalar_mul(proof.A, r)) && iszero(scalar_mul(proof.B, r)) && iszero(scalar_mul(proof.C, r)))
        return false
    end

    # Multi-pairing product, single final exponentiation
    engine = pvk.vk.engine
    prod = pairing_batch(engine,
        [proof.A, prepared_inputs, proof.C],
        [proof.B, pvk.gamma_g2_neg, pvk.delta_g2_neg],
    )
    return prod == pvk.alpha_g1_beta_g2
end

export PreparedVerificationKey, prepare_verifying_key, prepare_inputs, verify_with_prepared

"""
    validate_witness_shape(r1cs::R1CS{F}, witness::Witness{F}) where F

Validate witness layout conventions:
- `length(witness.values) == r1cs.num_vars`
- `witness.values[1] == 1`

Throws `ArgumentError` with explicit diagnostics if validation fails.
"""
function validate_witness_shape(r1cs::R1CS{F}, witness::Witness{F}) where F
    n = length(witness.values)
    if n != r1cs.num_vars
        throw(ArgumentError("Expected witness length $(r1cs.num_vars), got $(n)"))
    end
    if !isone(witness.values[1])
        throw(ArgumentError("Witness convention violation: values[1] must be one($F)"))
    end
    return nothing
end

"""
    public_inputs_from_witness(r1cs::R1CS{F}, witness::Witness{F}) where F

Extract arkworks-style public inputs from a witness, excluding the leading
constant `1` at index 1.
"""
function public_inputs_from_witness(r1cs::R1CS{F}, witness::Witness{F}) where F
    validate_witness_shape(r1cs, witness)
    expected = r1cs.num_public - 1
    if expected <= 0
        return F[]
    end
    return copy(witness.values[2:r1cs.num_public])
end

"""
    setup(r1cs::R1CS{F}; rng=Random.GLOBAL_RNG, engine=BN254_ENGINE, prepare_vk::Bool=false) where F

High-level Groth16 setup wrapper:
1. Converts `r1cs` to `qap`.
2. Runs trusted setup.
3. Optionally prepares the verifying key.

Returns a named tuple with fields `pk`, `vk`, `qap`, and optionally `pvk`.
"""
function setup(r1cs::R1CS{F}; rng::AbstractRNG=Random.GLOBAL_RNG, engine::AbstractPairingEngine=BN254_ENGINE, prepare_vk::Bool=false) where F
    qap = r1cs_to_qap(r1cs)
    keypair = setup_full(qap; rng=rng, engine=engine)
    if prepare_vk
        pvk = prepare_verifying_key(keypair.vk)
        return (pk=keypair.pk, vk=keypair.vk, qap=qap, pvk=pvk)
    end
    return (pk=keypair.pk, vk=keypair.vk, qap=qap)
end

"""
    prove(pk::ProvingKey, qap::QAP{F}, witness::Witness{F}; rng=Random.GLOBAL_RNG, debug_no_random=false) where F

High-level proving wrapper over `prove_full`.
"""
function prove(pk::ProvingKey, qap::QAP{F}, witness::Witness{F}; rng::AbstractRNG=Random.GLOBAL_RNG, debug_no_random::Bool=false) where F
    if length(witness.values) != qap.num_vars
        throw(ArgumentError("Witness length $(length(witness.values)) does not match QAP num_vars $(qap.num_vars)"))
    end
    if !isone(witness.values[1])
        throw(ArgumentError("Witness convention violation: values[1] must be one($F)"))
    end
    return prove_full(pk, qap, witness; rng=rng, debug_no_random=debug_no_random)
end

"""
    process_vk(vk::VerificationKey)

Prepare a verification key for repeated verification.
"""
process_vk(vk::VerificationKey) = prepare_verifying_key(vk)

"""
    verify(vk::VerificationKey, public_inputs::Vector{F}, proof::Groth16Proof) where F

High-level verifier wrapper for non-prepared verification keys.
"""
function verify(vk::VerificationKey, public_inputs::Vector{F}, proof::Groth16Proof) where F
    return verify_full(vk, proof, public_inputs)
end

"""
    verify(pvk::PreparedVerificationKey, public_inputs::Vector{F}, proof::Groth16Proof) where F

High-level verifier wrapper for prepared verification keys.
"""
function verify(pvk::PreparedVerificationKey, public_inputs::Vector{F}, proof::Groth16Proof) where F
    prepared = prepare_inputs(pvk, public_inputs)
    return verify_with_prepared(pvk, proof, prepared)
end

"""
    verify_prepared(pvk::PreparedVerificationKey, prepared_inputs::G1Point, proof::Groth16Proof)

Verify a proof using already prepared public inputs.
"""
function verify_prepared(pvk::PreparedVerificationKey, prepared_inputs::G1Point, proof::Groth16Proof)
    return verify_with_prepared(pvk, proof, prepared_inputs)
end

export validate_witness_shape, public_inputs_from_witness
export setup, prove, process_vk, verify, verify_prepared
