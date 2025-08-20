"""
BN254 elliptic curve implementation.

The BN254 curve has:
- G1: y² = x³ + 3 over Fp
- G2: y² = x³ + 3/(9+u) over Fp2
- Pairing: e: G1 × G2 → GT
"""

using GrothAlgebra
using StaticArrays

include("BN254Fields.jl")

# Curve parameters
struct BN254Curve <: AbstractCurve end

"""
    G1Point

Point on the BN254 G1 curve: y² = x³ + 3 over Fp.
Uses Jacobian coordinates (X, Y, Z) where x = X/Z² and y = Y/Z³.
"""
struct G1Point <: GroupElem{BN254Curve}
    coords::SVector{3, BN254Field}
    
    function G1Point(X::BN254Field, Y::BN254Field, Z::BN254Field)
        new(SVector(X, Y, Z))
    end
end

# Convenient constructors
G1Point(x::Integer, y::Integer, z::Integer) = G1Point(bn254_field(x), bn254_field(y), bn254_field(z))
G1Point(x::Integer, y::Integer) = G1Point(bn254_field(x), bn254_field(y), one(BN254Field))

# Access coordinates
Base.getindex(p::G1Point, i::Int) = p.coords[i]
x_coord(p::G1Point) = p[1]
y_coord(p::G1Point) = p[2]
z_coord(p::G1Point) = p[3]

# Convert to affine coordinates
function to_affine(p::G1Point)
    if iszero(p)
        return (zero(BN254Field), zero(BN254Field))
    end
    z_inv = inv(z_coord(p))
    z_inv2 = z_inv^2
    z_inv3 = z_inv2 * z_inv
    return (x_coord(p) * z_inv2, y_coord(p) * z_inv3)
end

# Identity element (point at infinity in Jacobian coordinates)
Base.zero(::Type{G1Point}) = G1Point(one(BN254Field), one(BN254Field), zero(BN254Field))
Base.zero(::G1Point) = zero(G1Point)
Base.iszero(p::G1Point) = iszero(z_coord(p))

# Check if point is on curve
function is_on_curve(p::G1Point)
    if iszero(p)
        return true
    end
    X, Y, Z = x_coord(p), y_coord(p), z_coord(p)
    Z2 = Z^2
    Z4 = Z2^2
    Z6 = Z4 * Z2
    # Check Y² = X³ + 3*Z⁶ in Jacobian coordinates
    return Y^2 == X^3 + bn254_field(3) * Z6
end

"""
    double(p::G1Point)

Point doubling in Jacobian coordinates.
"""
function double(p::G1Point)
    if iszero(p)
        return p
    end
    
    X, Y, Z = x_coord(p), y_coord(p), z_coord(p)
    
    # Using efficient doubling formulas
    S = bn254_field(4) * X * Y^2
    M = bn254_field(3) * X^2  # For curve y² = x³ + b, we have a = 0
    X3 = M^2 - bn254_field(2) * S
    Y3 = M * (S - X3) - bn254_field(8) * Y^4
    Z3 = bn254_field(2) * Y * Z
    
    return G1Point(X3, Y3, Z3)
end

"""
    +(p::G1Point, q::G1Point)

Point addition in Jacobian coordinates.
"""
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

# Scalar multiplication (use the one from GrothAlgebra.Group)
# It will automatically use our + operation

"""
    G2Point

Point on the BN254 G2 curve over Fp2.
"""
struct G2Point <: GroupElem{BN254Curve}
    coords::SVector{3, Fp2Element}
    
    function G2Point(X::Fp2Element, Y::Fp2Element, Z::Fp2Element)
        new(SVector(X, Y, Z))
    end
end

# Convenient constructors
function G2Point(x0::Integer, x1::Integer, y0::Integer, y1::Integer, z0::Integer, z1::Integer)
    X = Fp2Element(bn254_field(x0), bn254_field(x1))
    Y = Fp2Element(bn254_field(y0), bn254_field(y1))
    Z = Fp2Element(bn254_field(z0), bn254_field(z1))
    G2Point(X, Y, Z)
end

function G2Point(x0::Integer, x1::Integer, y0::Integer, y1::Integer)
    X = Fp2Element(bn254_field(x0), bn254_field(x1))
    Y = Fp2Element(bn254_field(y0), bn254_field(y1))
    Z = one(Fp2Element)
    G2Point(X, Y, Z)
end

# Access coordinates
Base.getindex(p::G2Point, i::Int) = p.coords[i]
x_coord(p::G2Point) = p[1]
y_coord(p::G2Point) = p[2]
z_coord(p::G2Point) = p[3]

# Identity element
Base.zero(::Type{G2Point}) = G2Point(one(Fp2Element), one(Fp2Element), zero(Fp2Element))
Base.zero(::G2Point) = zero(G2Point)
Base.iszero(p::G2Point) = iszero(z_coord(p))

# Convert to affine
function to_affine(p::G2Point)
    if iszero(p)
        return (zero(Fp2Element), zero(Fp2Element))
    end
    z_inv = inv(z_coord(p))
    z_inv2 = z_inv^2
    z_inv3 = z_inv2 * z_inv
    return (x_coord(p) * z_inv2, y_coord(p) * z_inv3)
end

"""
    double(p::G2Point)

Point doubling for G2.
"""
function double(p::G2Point)
    if iszero(p)
        return p
    end
    
    X, Y, Z = x_coord(p), y_coord(p), z_coord(p)
    
    S = Fp2Element(4) * X * Y^2
    M = Fp2Element(3) * X^2
    X3 = M^2 - Fp2Element(2) * S
    Y3 = M * (S - X3) - Fp2Element(8) * Y^4
    Z3 = Fp2Element(2) * Y * Z
    
    return G2Point(X3, Y3, Z3)
end

"""
    +(p::G2Point, q::G2Point)

Point addition for G2.
"""
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

# Generator points for BN254
"""
    g1_generator()

Returns the standard generator for G1.
"""
function g1_generator()
    return G1Point(1, 2)
end

"""
    g2_generator()

Returns the standard generator for G2.
"""
function g2_generator()
    # Standard generator coordinates for BN254 G2
    x0 = parse(BigInt, "10857046999023057135944570762232829481370756359578518086990519993285655852781")
    x1 = parse(BigInt, "11559732032986387107991004021392285783925812861821192530917403151452391805634")
    y0 = parse(BigInt, "8495653923123431417604973247489272438418190587263600148770280649306958101930")
    y1 = parse(BigInt, "4082367875863433681332203403145435568316851327593401208105741076214120093531")
    
    return G2Point(
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