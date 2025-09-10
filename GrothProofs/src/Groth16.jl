"""
Groth16 proof system implementation (production-style scaffold).

Provides Proving/Verification keys and standard Groth16 equations.
The implementation currently uses dense polynomial paths (no FFT yet).
"""

using GrothAlgebra
using GrothCurves
using Random

# R1CS and QAP are included at the package level (GrothProofs.jl)

"""
Internal helper to convert a field element to an Integer scalar for scalar_mul.
"""
_to_int(x) = Integer(x.value)

"""
    Groth16Proof

A Groth16 proof consists of three group elements.
"""
struct Groth16Proof
    A::G1Point  # [A]₁
    B::G2Point  # [B]₂
    C::G1Point  # [C]₁
end

"""
Production-style key structures for Groth16.

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

struct VerificationKey
    alpha_g1::G1Point
    beta_g2::G2Point
    gamma_g2::G2Point
    delta_g2::G2Point
    IC::Vector{G1Point}  # [ (β·u_i(τ) + α·v_i(τ) + w_i(τ)) / γ ]₁ for public i (including 1)
end

struct Keypair
    pk::ProvingKey
    vk::VerificationKey
end

# Export types
export Groth16Proof

"""
    setup_full(qap::QAP{F}; rng=Random.GLOBAL_RNG) where F

Generate production-style Groth16 keys (ProvingKey, VerificationKey).
This implementation uses the current dense polynomial machinery (no FFT yet).
"""
function setup_full(qap::QAP{F}; rng::AbstractRNG=Random.GLOBAL_RNG) where F
    # Sample toxic waste
    α = F(rand(rng, 1:1000000))
    β = F(rand(rng, 1:1000000))
    γ = F(rand(rng, 1:1000000))
    δ = F(rand(rng, 1:1000000))
    τ = F(rand(rng, 1:1000000))

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

    # Map evaluations into group queries
    A_query_g1 = [scalar_mul(g1, _to_int(A_eval[i])) for i in 1:m]
    B_query_g2 = [scalar_mul(g2, _to_int(B_eval[i])) for i in 1:m]
    B_query_g1 = [scalar_mul(g1, _to_int(B_eval[i])) for i in 1:m]
    C_query_g1 = [scalar_mul(g1, _to_int(C_eval[i])) for i in 1:m]

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
    H_query_g1 = [scalar_mul(g1, _to_int(t_tau * tau_powers[k] * δ_inv)) for k in 1:max_k]

    # L_query for private variables only: [(β·u_i(τ) + α·v_i(τ) + w_i(τ)) / δ]₁
    num_public = qap.num_public
    L_query_g1 = G1Point[]
    for i in (num_public+1):m
        acc = β * A_eval[i] + α * B_eval[i] + C_eval[i]
        push!(L_query_g1, scalar_mul(g1, _to_int(acc * δ_inv)))
    end

    # IC for public inputs: [(β·u_i(τ) + α·v_i(τ) + w_i(τ)) / γ]₁, including 1
    IC = G1Point[]
    for i in 1:num_public
        acc = β * A_eval[i] + α * B_eval[i] + C_eval[i]
        push!(IC, scalar_mul(g1, _to_int(acc * γ_inv)))
    end

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

    vk = VerificationKey(
        alpha_g1,
        beta_g2,
        gamma_g2,
        delta_g2,
        IC,
    )

    return Keypair(pk, vk)
end

"""
    prove_full(pk::ProvingKey, qap::QAP{F}, witness::Witness{F}; rng=Random.GLOBAL_RNG) where F

Generate a Groth16 proof using the production-style proving key.
Uses current dense polynomial path for h and H.
"""
function prove_full(pk::ProvingKey, qap::QAP{F}, witness::Witness{F}; rng::AbstractRNG=Random.GLOBAL_RNG, debug_no_random::Bool=false) where F
    w_vals = witness.values
    m = qap.num_vars

    # Randomizers
    r = debug_no_random ? zero(F) : F(rand(rng, 1:1000000))
    s = debug_no_random ? zero(F) : F(rand(rng, 1:1000000))

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

    # h(x): compute via dense division path and map coefficients via H_query
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

Verify a Groth16 proof using the production-style verification key.
Checks the standard equation: e(A,B) == e(α,β) · e(vk_x,γ) · e(C,δ).
Performs on-curve and subgroup checks for A, B, C.
"""
function verify_full(vk::VerificationKey, proof::Groth16Proof, public_inputs::Vector{F}) where F
    # On-curve checks
    if !(GrothCurves.is_on_curve(proof.A) && GrothCurves.is_on_curve(proof.B) && GrothCurves.is_on_curve(proof.C))
        return false
    end

    # Subgroup checks: multiply by r and ensure infinity
    r = GrothCurves.BN254_ORDER_R
    if !(iszero(scalar_mul(proof.A, r)) && iszero(scalar_mul(proof.B, r)) && iszero(scalar_mul(proof.C, r)))
        return false
    end

    # Compute vk_x = IC[1] + Σ_{i=2..num_public} input[i] * IC[i]
    # public_inputs is expected to include the leading 1 at index 1
    if length(public_inputs) < 1 || length(vk.IC) < 1
        return false
    end
    if length(public_inputs) > length(vk.IC)
        return false
    end
    if length(public_inputs) > 1
        pts_ic = vk.IC[2:length(public_inputs)]
        scal_ic = [_to_int(public_inputs[i]) for i in 2:length(public_inputs)]
        vk_x = vk.IC[1] + GrothAlgebra.multi_scalar_mul(pts_ic, scal_ic)
    else
        vk_x = vk.IC[1]
    end

    lhs = pairing(proof.A, proof.B)
    rhs = pairing(vk.alpha_g1, vk.beta_g2) * pairing(vk_x, vk.gamma_g2) * pairing(proof.C, vk.delta_g2)
    return lhs == rhs
end

# Export production-style API
export ProvingKey, VerificationKey, Keypair
export setup_full, prove_full, verify_full

"""
Prepared verification key caching e(α,β) and holding γ, δ for prepared verification.
"""
struct PreparedVerificationKey
    vk::VerificationKey
    alpha_g1_beta_g2::GTElement
    gamma_g2_neg::G2Point
    delta_g2_neg::G2Point
end

"""
    prepare_verifying_key(vk::VerificationKey)

Compute cached pairing e(α,β) and negated γ, δ for prepared verification.
"""
function prepare_verifying_key(vk::VerificationKey)
    return PreparedVerificationKey(
        vk,
        pairing(vk.alpha_g1, vk.beta_g2),
        -vk.gamma_g2,
        -vk.delta_g2,
    )
end

"""
    prepare_inputs(pvk::PreparedVerificationKey, public_inputs::Vector{F}) where F

Compute vk_x = IC[1] + Σ_{i=2..} input[i] * IC[i].
public_inputs must include the leading 1 at index 1.
"""
function prepare_inputs(pvk::PreparedVerificationKey, public_inputs::Vector{F}) where F
    vk = pvk.vk
    if isempty(public_inputs) || isempty(vk.IC) || length(public_inputs) > length(vk.IC)
        throw(ArgumentError("Invalid public inputs for prepared inputs"))
    end
    acc = vk.IC[1]
    for i in 2:length(public_inputs)
        xi = public_inputs[i]
        iszero(xi) && continue
        acc += scalar_mul(vk.IC[i], _to_int(xi))
    end
    return acc
end

"""
    verify_with_prepared(pvk::PreparedVerificationKey, proof::Groth16Proof, prepared_inputs::G1Point)

Verify using prepared verifier: compare product pairing to cached e(α,β).
Includes on-curve and subgroup checks as in verify_full.
"""
function verify_with_prepared(pvk::PreparedVerificationKey, proof::Groth16Proof, prepared_inputs::G1Point)
    # On-curve and subgroup checks
    if !(GrothCurves.is_on_curve(proof.A) && GrothCurves.is_on_curve(proof.B) && GrothCurves.is_on_curve(proof.C))
        return false
    end
    r = GrothCurves.BN254_ORDER_R
    if !(iszero(scalar_mul(proof.A, r)) && iszero(scalar_mul(proof.B, r)) && iszero(scalar_mul(proof.C, r)))
        return false
    end

    # Multi-pairing product, single final exponentiation
    prod = pairing_batch(
        [proof.A, prepared_inputs, proof.C],
        [proof.B, pvk.gamma_g2_neg, pvk.delta_g2_neg],
    )
    return prod == pvk.alpha_g1_beta_g2
end

export PreparedVerificationKey, prepare_verifying_key, prepare_inputs, verify_with_prepared
