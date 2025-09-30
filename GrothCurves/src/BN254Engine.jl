"""
    BN254Engine

Pairing engine singleton implementing the BN254 pairing interface.
"""
struct BN254Engine <: AbstractPairingEngine{BN254Curve} end

const BN254_ENGINE = BN254Engine()

export BN254Engine, BN254_ENGINE
