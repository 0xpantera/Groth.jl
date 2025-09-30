# BN254 elliptic curve implementation.
#
# Defines G1/G2 groups over BN254 along with projective helpers.

# Curve tag used by the generic projective point representation
struct BN254Curve <: AbstractCurve end

# BN254 scalar field / subgroup order r (Q)
# This is the prime order of the G1/G2 subgroups used by Groth16 (aka Fr modulus).
const BN254_ORDER_R = parse(BigInt, "21888242871839275222246405745257275088548364400416034343698204186575808495617")

# Shared projective point storage (X:Y:Z) in Jacobian coordinates
const G1Point = ProjectivePoint{BN254Curve, BN254Field}
const G2Point = ProjectivePoint{BN254Curve, Fp2Element}

# Twist coefficient for the BN254 D-twist (y² = x³ + b')
const G2_B_TWIST = Fp2Element(3, 0) / Fp2Element(9, 1)

# Convenient constructors for G1
G1Point(x::Integer, y::Integer, z::Integer) =
    G1Point(bn254_field(x), bn254_field(y), bn254_field(z))
G1Point(x::Integer, y::Integer) = G1Point(bn254_field(x), bn254_field(y), one(BN254Field))

# Convenient constructors for G2
function G2Point(x0::Integer, x1::Integer, y0::Integer, y1::Integer, z0::Integer, z1::Integer)
    X = Fp2Element(bn254_field(x0), bn254_field(x1))
    Y = Fp2Element(bn254_field(y0), bn254_field(y1))
    Z = Fp2Element(bn254_field(z0), bn254_field(z1))
    G2Point(X, Y, Z)
end

function G2Point(x0::Integer, x1::Integer, y0::Integer, y1::Integer)
    X = Fp2Element(bn254_field(x0), bn254_field(x1))
    Y = Fp2Element(bn254_field(y0), bn254_field(y1))
    G2Point(X, Y, one(Fp2Element))
end

G2Point(X::Fp2Element, Y::Fp2Element) = G2Point(X, Y, one(Fp2Element))

# Affine conversion specialised for each field type
function to_affine(p::G1Point)
    if iszero(p)
        return (zero(BN254Field), zero(BN254Field))
    end
    z_inv = inv(z_coord(p))
    z_inv2 = z_inv^2
    z_inv3 = z_inv2 * z_inv
    return (x_coord(p) * z_inv2, y_coord(p) * z_inv3)
end

function to_affine(p::G2Point)
    if iszero(p)
        return (zero(Fp2Element), zero(Fp2Element))
    end
    z_inv = inv(z_coord(p))
    z_inv2 = z_inv^2
    z_inv3 = z_inv2 * z_inv
    return (x_coord(p) * z_inv2, y_coord(p) * z_inv3)
end

# Curve membership checks ---------------------------------------------------

function is_on_curve(p::G1Point)
    if iszero(p)
        return true
    end
    X, Y, Z = x_coord(p), y_coord(p), z_coord(p)
    Z2 = Z^2
    Z4 = Z2^2
    Z6 = Z4 * Z2
    return Y^2 == X^3 + bn254_field(3) * Z6
end

function is_on_curve(p::G2Point)
    if iszero(p)
        return true
    end
    X, Y, Z = x_coord(p), y_coord(p), z_coord(p)
    Z2 = Z^2
    Z4 = Z2^2
    Z6 = Z4 * Z2
    return Y^2 == X^3 + G2_B_TWIST * Z6
end

# G1 arithmetic -------------------------------------------------------------

function double(p::G1Point)
    if iszero(p)
        return p
    end

    X, Y, Z = x_coord(p), y_coord(p), z_coord(p)

    S = bn254_field(4) * X * Y^2
    M = bn254_field(3) * X^2
    X3 = M^2 - bn254_field(2) * S
    Y3 = M * (S - X3) - bn254_field(8) * Y^4
    Z3 = bn254_field(2) * Y * Z

    return G1Point(X3, Y3, Z3)
end

function Base.:+(p::G1Point, q::G1Point)
    if iszero(p)
        return q
    elseif iszero(q)
        return p
    end

    X1, Y1, Z1 = x_coord(p), y_coord(p), z_coord(p)
    X2, Y2, Z2 = x_coord(q), y_coord(q), z_coord(q)

    Z1Z1 = Z1^2
    Z2Z2 = Z2^2
    U1 = X1 * Z2Z2
    U2 = X2 * Z1Z1
    S1 = Y1 * Z2 * Z2Z2
    S2 = Y2 * Z1 * Z1Z1

    if U1 == U2
        if S1 == S2
            return double(p)
        else
            return zero(G1Point)
        end
    end

    H = U2 - U1
    I = (bn254_field(2) * H)^2
    J = H * I
    r = bn254_field(2) * (S2 - S1)
    V = U1 * I
    X3 = r^2 - J - bn254_field(2) * V
    Y3 = r * (V - X3) - bn254_field(2) * S1 * J
    Z3 = ((Z1 + Z2)^2 - Z1Z1 - Z2Z2) * H

    return G1Point(X3, Y3, Z3)
end

Base.:-(p::G1Point) = G1Point(x_coord(p), -y_coord(p), z_coord(p))
Base.:-(p::G1Point, q::G1Point) = p + (-q)

# G2 parameters -------------------------------------------------------------

# BN254 parameter u = 4965661367192848881 (from BN254MillerLoop)
# G2 arithmetic -------------------------------------------------------------

function double(p::G2Point)
    if iszero(p)
        return p
    end

    X, Y, Z = x_coord(p), y_coord(p), z_coord(p)

    XX = X^2
    ZZ = Z^2
    ZZ2 = ZZ^2
    S = Fp2Element(4) * X * Y^2
    M = Fp2Element(3) * XX
    X3 = M^2 - Fp2Element(2) * S
    Y3 = M * (S - X3) - Fp2Element(8) * Y^4
    Z3 = Fp2Element(2) * Y * Z

    return G2Point(X3, Y3, Z3)
end

function Base.:+(p::G2Point, q::G2Point)
    if iszero(p)
        return q
    elseif iszero(q)
        return p
    end

    X1, Y1, Z1 = x_coord(p), y_coord(p), z_coord(p)
    X2, Y2, Z2 = x_coord(q), y_coord(q), z_coord(q)

    Z1Z1 = Z1^2
    Z2Z2 = Z2^2
    U1 = X1 * Z2Z2
    U2 = X2 * Z1Z1
    S1 = Y1 * Z2 * Z2Z2
    S2 = Y2 * Z1 * Z1Z1

    if U1 == U2
        if S1 == S2
            return double(p)
        else
            return zero(G2Point)
        end
    end

    H = U2 - U1
    I = (Fp2Element(2) * H)^2
    J = H * I
    r = Fp2Element(2) * (S2 - S1)
    V = U1 * I
    X3 = r^2 - J - Fp2Element(2) * V
    Y3 = r * (V - X3) - Fp2Element(2) * S1 * J
    Z3 = ((Z1 + Z2)^2 - Z1Z1 - Z2Z2) * H

    return G2Point(X3, Y3, Z3)
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
        Fp2Element(bn254_field(x0), bn254_field(x1)),
        Fp2Element(bn254_field(y0), bn254_field(y1)),
        one(Fp2Element)
    )
end

# Export types and functions
export BN254Curve, G1Point, G2Point
export to_affine, is_on_curve, double
export g1_generator, g2_generator
export x_coord, y_coord, z_coord
export BN254_ORDER_R

# Batch normalisation ------------------------------------------------------

function batch_to_affine!(pts::Vector{G1Point})
    n = length(pts)
    if n == 0
        return pts
    end
    idxs = [i for i in 1:n if !iszero(pts[i])]
    if isempty(idxs)
        return pts
    end
    prefix = Vector{BN254Field}(undef, length(idxs))
    prefix[1] = z_coord(pts[idxs[1]])
    for k in 2:length(idxs)
        prefix[k] = prefix[k-1] * z_coord(pts[idxs[k]])
    end
    inv_total = inv(prefix[end])
    zinv = Vector{BN254Field}(undef, length(idxs))
    for k in length(idxs):-1:1
        z_k = z_coord(pts[idxs[k]])
        prev = k == 1 ? one(BN254Field) : prefix[k-1]
        zinv_k = inv_total * prev
        zinv[k] = zinv_k
        inv_total = inv_total * z_k
    end
    for (t, i) in enumerate(idxs)
        P = pts[i]
        z_inv = zinv[t]
        z_inv2 = z_inv^2
        z_inv3 = z_inv2 * z_inv
        pts[i] = G1Point(x_coord(P) * z_inv2, y_coord(P) * z_inv3, one(BN254Field))
    end
    return pts
end

function batch_to_affine!(pts::Vector{G2Point})
    n = length(pts)
    if n == 0
        return pts
    end
    idxs = [i for i in 1:n if !iszero(pts[i])]
    if isempty(idxs)
        return pts
    end
    prefix = Vector{Fp2Element}(undef, length(idxs))
    prefix[1] = z_coord(pts[idxs[1]])
    for k in 2:length(idxs)
        prefix[k] = prefix[k-1] * z_coord(pts[idxs[k]])
    end
    inv_total = inv(prefix[end])
    zinv = Vector{Fp2Element}(undef, length(idxs))
    for k in length(idxs):-1:1
        z_k = z_coord(pts[idxs[k]])
        prev = k == 1 ? one(Fp2Element) : prefix[k-1]
        zinv_k = inv_total * prev
        zinv[k] = zinv_k
        inv_total = inv_total * z_k
    end
    for (t, i) in enumerate(idxs)
        P = pts[i]
        z_inv = zinv[t]
        z_inv2 = z_inv^2
        z_inv3 = z_inv2 * z_inv
        pts[i] = G2Point(x_coord(P) * z_inv2, y_coord(P) * z_inv3, one(Fp2Element))
    end
    return pts
end
