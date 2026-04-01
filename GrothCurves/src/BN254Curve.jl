# BN254 elliptic curve implementation.
#
# Defines G1/G2 groups over BN254 along with projective helpers.

# Curve tag used by the generic projective point representation
struct BN254Curve <: AbstractCurve end

# BN254 scalar field / subgroup order r (Q)
# This is the prime order of the G1/G2 subgroups used by Groth16 (aka Fr modulus).
const BN254_ORDER_R = parse(BigInt, "21888242871839275222246405745257275088548364400416034343698204186575808495617")

# Shared projective point storage (X:Y:Z) in Jacobian coordinates
const G1Point = ProjectivePoint{BN254Curve, BN254Fq}
const G2Point = ProjectivePoint{BN254Curve, Fp2Element}

const BN254_GLV_BETA = bn254_fq(parse(BigInt, "21888242871839275220042445260109153167277707414472061641714758635765020556616"))
const BN254_GLV_BETA_FP2 = Fp2Element(BN254_GLV_BETA, zero(BN254Fq))
const BN254_G1_GLV_LAMBDA = parse(BigInt, "21888242871839275217838484774961031246154997185409878258781734729429964517155")
const BN254_G2_GLV_LAMBDA = parse(BigInt, "4407920970296243842393367215006156084916469457145843978461")
const BN254_G1_GLV_DECOMP_MATRIX = (
    -parse(BigInt, "147946756881789319000765030803803410728"),
    parse(BigInt, "9931322734385697763"),
    -parse(BigInt, "9931322734385697763"),
    -parse(BigInt, "147946756881789319010696353538189108491"),
)
const BN254_G2_GLV_DECOMP_MATRIX = (
    -parse(BigInt, "147946756881789319010696353538189108491"),
    -parse(BigInt, "9931322734385697763"),
    parse(BigInt, "9931322734385697763"),
    -parse(BigInt, "147946756881789319000765030803803410728"),
)
# Keep small and medium scalars on the existing w-NAF path; the Stage 7A
# scalar sweep shows GLV only starts to pull ahead on larger bitlengths.
const BN254_G1_GLV_THRESHOLD_BITS = 192

# Stage 7 tuning on the Montgomery backend found that BN254 MSMs prefer
# smaller windows than the generic schedule at the benchmark sizes we exercise
# in the prover and benchmark suite. G2 prover query MSMs are especially
# sensitive at small sizes, where the real prove_full workload prefers w=2.
function _pippenger_window(::Type{G1Point}, size::Int)
    if size <= 32
        return 3
    elseif size <= 128
        return 5
    elseif size >= 512
        return 6
    end
    return GrothAlgebra._default_pippenger_window(size)
end

function _pippenger_window(::Type{G2Point}, size::Int)
    if size <= 32
        return 2
    elseif size <= 128
        return 5
    elseif size >= 512
        return 6
    end
    return GrothAlgebra._default_pippenger_window(size)
end

# Single-point scalar multiplication uses the same generic API, but BN254
# points benefit from a tuned w-NAF window instead of the binary fallback.
_scalar_mul_window(::Type{G1Point}, bit_length::Int) = bit_length <= 64 ? 3 : 4
_scalar_mul_window(::Type{G2Point}, ::Int) = 3

const G1_B = bn254_fq(3)
const G1_AFFINE_Z = one(BN254Fq)
const G2_AFFINE_Z = one(Fp2Element)

# Twist coefficient for the BN254 D-twist (y^2 = x^3 + b')
const G2_B_TWIST = Fp2Element(3, 0) / Fp2Element(9, 1)

# Convenient constructors for G1
G1Point(x::Integer, y::Integer, z::Integer) =
    G1Point(bn254_fq(x), bn254_fq(y), bn254_fq(z))
G1Point(x::Integer, y::Integer) = G1Point(bn254_fq(x), bn254_fq(y), G1_AFFINE_Z)

# Convenient constructors for G2
function G2Point(x0::Integer, x1::Integer, y0::Integer, y1::Integer, z0::Integer, z1::Integer)
    X = Fp2Element(bn254_fq(x0), bn254_fq(x1))
    Y = Fp2Element(bn254_fq(y0), bn254_fq(y1))
    Z = Fp2Element(bn254_fq(z0), bn254_fq(z1))
    G2Point(X, Y, Z)
end

function G2Point(x0::Integer, x1::Integer, y0::Integer, y1::Integer)
    X = Fp2Element(bn254_fq(x0), bn254_fq(x1))
    Y = Fp2Element(bn254_fq(y0), bn254_fq(y1))
    G2Point(X, Y, G2_AFFINE_Z)
end

G2Point(X::Fp2Element, Y::Fp2Element) = G2Point(X, Y, G2_AFFINE_Z)

@inline _square(x) = x * x
@inline _square(x::Fp2Element) = square(x)
@inline _glv_bit_length(k::BigInt) = iszero(k) ? 0 : ndigits(k, base=2)

@inline function _glv_decomposition_bits(k1::BigInt, k2::BigInt)
    return max(_glv_bit_length(k1), _glv_bit_length(k2))
end

function _glv_scalar_decomposition(k::Integer, coeffs::NTuple{4,BigInt})
    scalar = mod(BigInt(k), BN254_ORDER_R)
    iszero(scalar) && return ((true, BigInt(0)), (true, BigInt(0)))

    n11, n12, n21, n22 = coeffs

    beta_1, rem_1 = divrem(scalar * n22, BN254_ORDER_R)
    if rem_1 + rem_1 > BN254_ORDER_R
        beta_1 += 1
    end

    beta_2, rem_2 = divrem(scalar * (-n12), BN254_ORDER_R)
    if rem_2 + rem_2 > BN254_ORDER_R
        beta_2 += 1
    end

    b1 = beta_1 * n11 + beta_2 * n21
    b2 = beta_1 * n12 + beta_2 * n22
    k1 = scalar - b1
    k2 = -b2

    return ((k1 >= 0, abs(k1)), (k2 >= 0, abs(k2)))
end

"""
    glv_scalar_decomposition(::Type{G1Point}, k::Integer)
    glv_scalar_decomposition(::Type{G2Point}, k::Integer)

Decompose a scalar into the signed two-term representation used by the BN254
GLV scalar-multiplication path for the selected group.
"""
glv_scalar_decomposition(::Type{G1Point}, k::Integer) = _glv_scalar_decomposition(k, BN254_G1_GLV_DECOMP_MATRIX)
glv_scalar_decomposition(::Type{G2Point}, k::Integer) = _glv_scalar_decomposition(k, BN254_G2_GLV_DECOMP_MATRIX)

glv_lambda(::Type{G1Point}) = BN254_G1_GLV_LAMBDA
glv_lambda(::Type{G2Point}) = BN254_G2_GLV_LAMBDA

"""
    glv_endomorphism(p::G1Point)
    glv_endomorphism(p::G2Point)

Apply the BN254 GLV endomorphism to a projective point in the corresponding
subgroup.
"""
@inline function glv_endomorphism(p::G1Point)
    iszero(p) && return p
    return G1Point(x_coord(p) * BN254_GLV_BETA, y_coord(p), z_coord(p))
end

@inline function glv_endomorphism(p::G2Point)
    iszero(p) && return p
    return G2Point(x_coord(p) * BN254_GLV_BETA_FP2, y_coord(p), z_coord(p))
end

function _glv_joint_scalar_mul(p::P, coeffs, endo_p::P) where {P<:ProjectivePoint{BN254Curve}}
    ((sgn_k1, k1), (sgn_k2, k2)) = coeffs
    iszero(k1) && iszero(k2) && return zero(p)

    b1 = sgn_k1 ? p : -p
    b2 = sgn_k2 ? endo_p : -endo_p
    b1b2 = b1 + b2

    max_bits = _glv_decomposition_bits(k1, k2)
    acc = zero(p)
    skip_leading_zeros = true

    for bit in (max_bits-1):-1:0
        bit1 = ((k1 >> bit) & 1) == 1
        bit2 = ((k2 >> bit) & 1) == 1

        if skip_leading_zeros
            if !(bit1 || bit2)
                continue
            end
            skip_leading_zeros = false
        else
            acc = acc + acc
        end

        if bit1 && bit2
            acc = acc + b1b2
        elseif bit1
            acc = acc + b1
        elseif bit2
            acc = acc + b2
        end
    end

    return acc
end

"""
    glv_scalar_mul(p::P, k::Integer) where {P<:ProjectivePoint{BN254Curve}}
    glv_scalar_mul(p::P, k::GrothAlgebra.BN254Fr) where {P<:ProjectivePoint{BN254Curve}}

Multiply a BN254 point by `k` using the Gallant-Lambert-Vanstone decomposition
and the curve-specific endomorphism.
"""
function glv_scalar_mul(p::P, k::Integer) where {P<:ProjectivePoint{BN254Curve}}
    if k == 0
        return zero(p)
    elseif k < 0
        return -glv_scalar_mul(p, -k)
    elseif k == 1
        return p
    end

    coeffs = glv_scalar_decomposition(P, k)
    return _glv_joint_scalar_mul(p, coeffs, glv_endomorphism(p))
end

@inline _bn254fr_scalar_bigint(k::GrothAlgebra.BN254Fr) =
    GrothAlgebra.limbs_to_bigint(GrothAlgebra.canonical_limbs(k))

function glv_scalar_mul(p::P, k::GrothAlgebra.BN254Fr) where {P<:ProjectivePoint{BN254Curve}}
    if iszero(k)
        return zero(p)
    elseif isone(k)
        return p
    end

    coeffs = glv_scalar_decomposition(P, _bn254fr_scalar_bigint(k))
    return _glv_joint_scalar_mul(p, coeffs, glv_endomorphism(p))
end

"""
    GrothAlgebra.scalar_mul(p::G1Point, k::Integer)
    GrothAlgebra.scalar_mul(p::G1Point, k::GrothAlgebra.BN254Fr)

BN254 G1 scalar multiplication dispatcher.

Small and medium scalars stay on the tuned w-NAF path, while larger scalars use
the GLV path once the measured threshold is reached.
"""
function GrothAlgebra.scalar_mul(p::G1Point, k::Integer)
    if k == 0
        return zero(p)
    elseif k == 1
        return p
    elseif k == -1
        return -p
    end

    bits = GrothAlgebra._bit_length(k)
    if bits < BN254_G1_GLV_THRESHOLD_BITS
        return GrothAlgebra.scalar_mul_wnaf(p, k, _scalar_mul_window(G1Point, bits))
    end

    return glv_scalar_mul(p, k)
end

function GrothAlgebra.scalar_mul(p::G1Point, k::GrothAlgebra.BN254Fr)
    if iszero(k)
        return zero(p)
    elseif isone(k)
        return p
    end

    bits = GrothAlgebra._bit_length(k)
    if bits < BN254_G1_GLV_THRESHOLD_BITS
        return GrothAlgebra.scalar_mul_wnaf(p, k, _scalar_mul_window(G1Point, bits))
    end

    return glv_scalar_mul(p, k)
end

@inline function _to_affine_coords(p::ProjectivePoint{BN254Curve,F}) where {F}
    if iszero(p)
        return (zero(F), zero(F))
    end

    X = x_coord(p)
    Y = y_coord(p)
    Z = z_coord(p)
    if isone(Z)
        return (X, Y)
    end

    z_inv = inv(Z)
    z_inv2 = _square(z_inv)
    return (X * z_inv2, Y * z_inv2 * z_inv)
end

# Affine conversion specialised for each field type
@inline to_affine(p::G1Point) = _to_affine_coords(p)
@inline to_affine(p::G2Point) = _to_affine_coords(p)

# Curve membership checks ---------------------------------------------------

@inline function _is_on_curve_jacobian(p::ProjectivePoint{BN254Curve,F}, b::F) where {F}
    if iszero(p)
        return true
    end

    X = x_coord(p)
    Y = y_coord(p)
    Z = z_coord(p)
    X2 = _square(X)
    Z2 = _square(Z)
    Z4 = _square(Z2)
    Z6 = Z4 * Z2
    return _square(Y) == X2 * X + b * Z6
end

is_on_curve(p::G1Point) = _is_on_curve_jacobian(p, G1_B)
is_on_curve(p::G2Point) = _is_on_curve_jacobian(p, G2_B_TWIST)

# Jacobian helpers ---------------------------------------------------------

@inline function _double_jacobian(::Type{P}, X, Y, Z) where {P}
    XX = _square(X)
    YY = _square(Y)
    YYYY = _square(YY)
    S = 4 * X * YY
    M = 3 * XX
    X3 = _square(M) - 2 * S
    Y3 = M * (S - X3) - 8 * YYYY
    Z3 = 2 * Y * Z
    return P(X3, Y3, Z3)
end

@inline function _add_jacobian(::Type{P}, p, q) where {P}
    X1 = x_coord(p)
    Y1 = y_coord(p)
    Z1 = z_coord(p)
    X2 = x_coord(q)
    Y2 = y_coord(q)
    Z2 = z_coord(q)

    Z1Z1 = _square(Z1)
    Z2Z2 = _square(Z2)
    U1 = X1 * Z2Z2
    U2 = X2 * Z1Z1
    S1 = Y1 * Z2 * Z2Z2
    S2 = Y2 * Z1 * Z1Z1

    if U1 == U2
        return S1 == S2 ? _double_jacobian(P, X1, Y1, Z1) : zero(P)
    end

    H = U2 - U1
    HH = _square(H)
    I = 4 * HH
    J = H * I
    r = 2 * (S2 - S1)
    V = U1 * I
    X3 = _square(r) - J - 2 * V
    Y3 = r * (V - X3) - 2 * S1 * J
    Zsum = Z1 + Z2
    Z3 = (_square(Zsum) - Z1Z1 - Z2Z2) * H
    return P(X3, Y3, Z3)
end

# G1 arithmetic -------------------------------------------------------------

function double(p::G1Point)
    iszero(p) && return p
    return _double_jacobian(G1Point, x_coord(p), y_coord(p), z_coord(p))
end

function Base.:+(p::G1Point, q::G1Point)
    if iszero(p)
        return q
    elseif iszero(q)
        return p
    end

    return _add_jacobian(G1Point, p, q)
end

Base.:-(p::G1Point) = G1Point(x_coord(p), -y_coord(p), z_coord(p))
Base.:-(p::G1Point, q::G1Point) = p + (-q)

# G2 parameters -------------------------------------------------------------

# BN254 parameter u = 4965661367192848881 (from BN254MillerLoop)
# G2 arithmetic -------------------------------------------------------------

function double(p::G2Point)
    iszero(p) && return p
    return _double_jacobian(G2Point, x_coord(p), y_coord(p), z_coord(p))
end

function Base.:+(p::G2Point, q::G2Point)
    if iszero(p)
        return q
    elseif iszero(q)
        return p
    end

    return _add_jacobian(G2Point, p, q)
end

Base.:-(p::G2Point) = G2Point(x_coord(p), -y_coord(p), z_coord(p))
Base.:-(p::G2Point, q::G2Point) = p + (-q)

# Generators ---------------------------------------------------------------

function g1_generator()
    G1Point(1, 2)
end

function g2_generator()
    x0 = parse(BigInt, "10857046999023057135944570762232829481370756359578518086990519993285655852781")
    x1 = parse(BigInt, "11559732032986387107991004021392285783925812861821192530917403151452391805634")
    y0 = parse(BigInt, "8495653923123431417604973247489272438418190587263600148770280649306958101930")
    y1 = parse(BigInt, "4082367875863433681332203403145435568316851327593401208105741076214120093531")
    G2Point(
        Fp2Element(bn254_fq(x0), bn254_fq(x1)),
        Fp2Element(bn254_fq(y0), bn254_fq(y1)),
        G2_AFFINE_Z,
    )
end

# Export types and functions
export BN254Curve, G1Point, G2Point
export to_affine, is_on_curve, double
export g1_generator, g2_generator
export x_coord, y_coord, z_coord
export BN254_ORDER_R

# Batch normalisation ------------------------------------------------------

function _batch_to_affine!(pts::Vector{P}, ::Type{F}) where {F,P<:ProjectivePoint{BN254Curve,F}}
    n = length(pts)
    n == 0 && return pts

    prefix = Vector{F}(undef, n)
    running = one(F)
    seen_nonzero = false

    @inbounds for i in 1:n
        point = pts[i]
        if !iszero(point)
            running = running * z_coord(point)
            seen_nonzero = true
        end
        prefix[i] = running
    end

    seen_nonzero || return pts

    inv_total = inv(running)
    affine_z = one(F)

    @inbounds for i in n:-1:1
        point = pts[i]
        iszero(point) && continue
        prev = i == 1 ? affine_z : prefix[i - 1]
        z_inv = inv_total * prev
        inv_total = inv_total * z_coord(point)
        z_inv2 = _square(z_inv)
        pts[i] = P(x_coord(point) * z_inv2, y_coord(point) * z_inv2 * z_inv, affine_z)
    end

    return pts
end

batch_to_affine!(pts::Vector{G1Point}) = _batch_to_affine!(pts, BN254Fq)
batch_to_affine!(pts::Vector{G2Point}) = _batch_to_affine!(pts, Fp2Element)
