"""
Pairing engine singleton for BN254, implementing the generic pairing interface.
"""

struct BN254Engine <: AbstractPairingEngine{BN254Curve} end

const BN254_ENGINE = BN254Engine()

export BN254Engine, BN254_ENGINE
