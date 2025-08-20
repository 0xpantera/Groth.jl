```
                    ╔═══════════════════════════════════════════════════╗
                    ║                                                   ║
                    ║    ▄████  ██▀███   ▒█████  ▄▄▄█████▓ ██░ ██     ║
                    ║   ██▒ ▀█▒▓██ ▒ ██▒▒██▒  ██▒▓  ██▒ ▓▒▓██░ ██▒    ║
                    ║  ▒██░▄▄▄░▓██ ░▄█ ▒▒██░  ██▒▒ ▓██░ ▒░▒██▀▀██░    ║
                    ║  ░▓█  ██▓▒██▀▀█▄  ▒██   ██░░ ▓██▓ ░ ░▓█ ░██     ║
                    ║  ░▒▓███▀▒░██▓ ▒██▒░ ████▓▒░  ▒██▒ ░ ░▓█▒░██▓    ║
                    ║   ░▒   ▒ ░ ▒▓ ░▒▓░░ ▒░▒░▒░   ▒ ░░    ▒ ░░▒░▒    ║
                    ║    ░   ░   ░▒ ░ ▒░  ░ ▒ ▒░     ░     ▒ ░▒░ ░    ║
                    ║  ░ ░   ░   ░░   ░ ░ ░ ░ ▒    ░       ░  ░░ ░    ║
                    ║        ░    ░         ░ ░             ░  ░  ░    ║
                    ╚═══════════════════════════════════════════════════╝
                          ∴‥∵‥∴ zkSNARK SUPREMACY ∴‥∵‥∴
                           
                    ░▒▓█►  ρяσνє єνєяутнιηg, яєνєαℓ ησтнιηg  ◄█▓▒░
                    
         ╔═════════════════════════════════════════════════════════════╗
         ║  "In cryptography we trust, in zero-knowledge we thrive"   ║
         ║              ∞ milady privacy maximalism ∞                 ║ 
         ╚═════════════════════════════════════════════════════════════╝
                                   
                        ▓▓▓▓▓▓▓   ▓▓▓▓▓▓▓
                       ▓       ▓ ▓       ▓    ┌─────────────┐
                       ▓  ╳ ╳  ▓ ▓  ╳ ╳  ▓    │ COMMITMENT  │
                       ▓       ▓ ▓       ▓    │   HIDDEN    │
                        ▓     ▓   ▓     ▓     │  KNOWLEDGE  │
                         ▓▓▓▓▓     ▓▓▓▓▓      └─────────────┘
                           ║         ║
                        ═══╬═════════╬═══ 
                           ║         ║
                     ╔═════╩═════════╩═════╗
                     ║  TRUSTED SETUP CULT  ║
                     ╚═════════════════════╝
```

# Groth.jl

**A modular, educational implementation of the Groth16 zero-knowledge proof system**  
*Built from scratch in Julia for learning, research, and cryptographic enlightenment*

## Project Structure

This is a monorepo containing several interconnected Julia packages:

### Core Packages

- **[GrothAlgebra](./GrothAlgebra)** - Foundation package with finite field arithmetic, polynomial operations, and group theory primitives
- **[GrothCurves](./GrothCurves)** - Elliptic curve implementations, focusing on BN254 (alt-bn128) with pairing support
- **[GrothProofs](./GrothProofs)** - Zero-knowledge proof systems including R1CS, QAP conversion, and Groth16 implementation
- **[GrothCrypto](./GrothCrypto)** - High-level cryptographic protocols built on top of the primitives
- **[GrothExamples](./GrothExamples)** - Educational examples and demonstrations

### Documentation & Tools

- **[docs/](./docs)** - Overall project documentation and theory explanations
- **[benchmarks/](./benchmarks)** - Performance benchmarks across packages

## Implementation Priorities

To build a complete Groth16 implementation following the QAP approach, here are the recommended implementation priorities:

### 1. Complete BN254 Curve in GrothCurves

- Implement base field Fp and extension field Fp2
- Add G1/G2 point operations (addition, doubling, scalar multiplication)
- Implement basic pairing map functionality

### 2. Implement R1CS Structure

- Define constraint matrices (A, B, C)
- Handle variable assignments (public inputs, private witnesses)
- Add constraint satisfaction checking

### 3. QAP Transformation

- Implement Lagrange interpolation of constraint matrices
- Compute polynomials A(x), B(x), C(x)
- Generate target polynomial t(x) = ∏(x - ωⁱ)

### 4. FFT for Efficient Polynomial Operations

- Find primitive roots of unity in the field
- Implement radix-2 FFT/IFFT algorithms
- Enable polynomial multiplication via convolution

### 5. Basic Prover and Verifier

- Implement simplified trusted setup ceremony
- Add proof generation algorithm
- Implement verification equation checks

## Current Status

- **Completed**: Field arithmetic, polynomial operations, group abstractions
- **In Progress**: Elliptic curve implementations
- **TODO**: R1CS, QAP transformation, trusted setup, proof system

## Development

This project follows Julia best practices with comprehensive testing and documentation. Each module has its own test suite and documentation.

## References

Following the RareSkills ZK Book for theoretical background and implementation guidance.
