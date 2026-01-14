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

"""
    _to_int(x)

Convert a field element to an `Integer` scalar for use with `scalar_mul`.
"""
_to_int(x) = Integer(x.value)

@inline function _rand_field(::Type{F}, rng::AbstractRNG) where F
    F === BN254ScalarField || throw(ArgumentError("No CSPRNG-ready sampler configured for field $F"))
    return F(rand(rng, 0:GrothCurves.BN254_ORDER_R - 1))
end

@inline function _rand_field_nonzero(::Type{F}, rng::AbstractRNG) where F
    F === BN254ScalarField || throw(ArgumentError("No CSPRNG-ready sampler configured for field $F"))
    return F(rand(rng, 1:GrothCurves.BN254_ORDER_R - 1))
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
    A_query_g1 = GrothAlgebra.batch_mul(g1tab, [_to_int(A_eval[i]) for i in 1:m])
    B_query_g2 = GrothAlgebra.batch_mul(g2tab, [_to_int(B_eval[i]) for i in 1:m])
    B_query_g1 = GrothAlgebra.batch_mul(g1tab, [_to_int(B_eval[i]) for i in 1:m])
    C_query_g1 = GrothAlgebra.batch_mul(g1tab, [_to_int(C_eval[i]) for i in 1:m])

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
    H_query_g1 = GrothAlgebra.batch_mul(g1tab, [_to_int(t_tau * tau_powers[k] * δ_inv) for k in 1:max_k])

    # L_query for private variables only: [(β·u_i(τ) + α·v_i(τ) + w_i(τ)) / δ]₁
    num_public = qap.num_public
    L_query_g1 = G1Point[]
    if m > num_public
        L_scalars = BigInt[]
        sizehint!(L_scalars, m - num_public)
        for i in (num_public+1):m
            acc = β * A_eval[i] + α * B_eval[i] + C_eval[i]
            push!(L_scalars, _to_int(acc * δ_inv))
        end
        L_query_g1 = GrothAlgebra.batch_mul(g1tab, L_scalars)
    end

    # IC for public inputs: [(β·u_i(τ) + α·v_i(τ) + w_i(τ)) / γ]₁, including 1
    IC_scalars = [_to_int((β * A_eval[i] + α * B_eval[i] + C_eval[i]) * γ_inv) for i in 1:num_public]
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
    alpha_g1 = scalar_mul(g1, _to_int(α))
    beta_g1  = scalar_mul(g1, _to_int(β))
    beta_g2  = scalar_mul(g2, _to_int(β))
    gamma_g2 = scalar_mul(g2, _to_int(γ))
    delta_g1 = scalar_mul(g1, _to_int(δ))
    delta_g2 = scalar_mul(g2, _to_int(δ))

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
    scalars = [_to_int(w_vals[i]) for i in 1:m]
    A_acc_g1 = GrothAlgebra.multi_scalar_mul(pk.A_query_g1, scalars)
    B_acc_g2 = GrothAlgebra.multi_scalar_mul(pk.B_query_g2, scalars)
    B_acc_g1 = GrothAlgebra.multi_scalar_mul(pk.B_query_g1, scalars)

    # A and B
    A1_g1 = pk.alpha_g1 + A_acc_g1
    B1_g1 = pk.beta_g1  + B_acc_g1
    A = A1_g1 + scalar_mul(pk.delta_g1, _to_int(r))
    B = pk.beta_g2  + B_acc_g2 + scalar_mul(pk.delta_g2, _to_int(s))

    # h(x): compute via coset FFT path and map coefficients via H_query
    h_poly = compute_h_polynomial(qap, witness)
    hk = length(h_poly.coeffs)
    pts_h = pk.H_query_g1[1:hk]
    scalars_h = [_to_int(c) for c in h_poly.coeffs]
    H = GrothAlgebra.multi_scalar_mul(pts_h, scalars_h)

    # L for private variables
    # L for private variables via MSM
    if m > pk.num_public
        priv_scalars = [_to_int(w_vals[i]) for i in (pk.num_public+1):m]
        L = GrothAlgebra.multi_scalar_mul(pk.L_query_g1, priv_scalars)
    else
        L = zero(G1Point)
    end

    # Cross terms in C (arkworks style): C = s*g_a + r*g1_b - r*s*δ + L + H
    # where g_a = A (includes r*δ), and g1_b = (β + Σ v_i) + s*δ in G1
    rs_delta = scalar_mul(pk.delta_g1, _to_int(r*s))
    g1_b_full = B1_g1 + scalar_mul(pk.delta_g1, _to_int(s))
    C = H + L + scalar_mul(g1_b_full, _to_int(r)) + scalar_mul(A, _to_int(s)) - rs_delta

    return Groth16Proof(A, B, C)
end

"""
    verify_full(vk::VerificationKey, proof::Groth16Proof, public_inputs::Vector{F}) where F

Verify a Groth16 proof with the production-style verification key.

Checks the standard pairing equation `e(A,B) == e(α,β) · e(vk_x,γ) · e(C,δ)` and
performs on-curve plus subgroup checks for `A`, `B`, and `C`.
"""
function verify_full(vk::VerificationKey{E}, proof::Groth16Proof, public_inputs::Vector{F}) where {E<:AbstractPairingEngine, F}
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
        scal_ic = [_to_int(public_inputs[i]) for i in 1:expected_len]
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

    Compute `vk_x = IC[1] + Σ input[i] * IC[i+1]`.

    Public inputs exclude the leading `1` (arkworks-style).
"""
function prepare_inputs(pvk::PreparedVerificationKey{E}, public_inputs::Vector{F}) where {E<:AbstractPairingEngine, F}
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
        acc += scalar_mul(vk.IC[i + 1], _to_int(xi))
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
