# Groth16 Implementation Cleanup Plan

## Current State Summary

### What's Working ✓
- R1CS representation with matrices L, R, O for `r = x * y * z * u`
- Witness generation and satisfaction checking
- QAP conversion via Lagrange interpolation
- Polynomial operations (add, multiply, evaluate, interpolate, divide)
- Basic finite field arithmetic using BigInt (BN254Field, Secp256k1Field)
- Example in `examples/test_r1cs_qap.jl` demonstrating R1CS → QAP workflow

### Technical Debt Issues
1. **Multiple field implementations coexist:**
   - `FieldElement.jl` - Original parametric type with UInt256 (kept for reference, not used)
   - `FiniteFields.jl` - Current active implementation using BigInt
   - `SimpleBN254Field.jl` - Orphaned implementation (kept for reference)

2. **Incomplete modules:**
   - `BN254Curve.jl` and `BN254Pairing.jl` written but commented out
   - `Groth16.jl` exists but depends on missing curve operations
   - `GrothCrypto` module is empty

3. **Test suite outdated:**
   - Tests still reference old UInt256 implementation
   - Need to be updated for BigInt-based fields

## Phase 1: Clean House ✅ (COMPLETED)

### ✅ Choose ONE field implementation
- **Decision:** Keep `FiniteFields.jl` with concrete types (BN254Field, Secp256k1Field)
- **Status:** COMPLETED - FiniteFields.jl is the active implementation
- Old implementations kept for reference but not included

### ✅ Update module structure  
- **Status:** COMPLETED - GrothAlgebra only includes FiniteFields.jl
- Added comments explaining archived files

### ✅ Fix test suite
- **Status:** COMPLETED - All 326 tests passing
  - `test_finite_fields.jl` - Tests for new field implementation (231 tests)
  - `test_groups.jl` - Group operation tests (65 tests)
  - `test_polynomials.jl` - Polynomial tests (30 tests)
  - `runtests.jl` - Main test runner (replaced with new version)
- Old test file archived as `runtests_old.jl`

### ✅ Clean up GrothCurves references
- **Status:** COMPLETED
- BN254Fields.jl updated to use GrothAlgebra's BN254Field
- Added comment about archived SimpleBN254Field.jl

## Phase 2: Complete the Curves

### Enable BN254 curve operations
1. Update `BN254Curve.jl` to use new BN254Field from FiniteFields
2. Fix G1Point and G2Point to work with BigInt-based fields
3. Re-enable in GrothCurves module
4. Test point addition, doubling, scalar multiplication

### Implement or mock pairing
1. Either implement simplified pairing or create mock for testing
2. Ensure pairing interface is correct for Groth16

## Phase 3: Complete Groth16

### Enable trusted setup
1. Uncomment and fix `Groth16.jl`
2. Update to use new field and curve types
3. Implement proper CRS generation

### Complete proof generation
1. Fix proof generation to use working components
2. Implement [A]₁, [B]₂, [C]₁ computation

### Add verification
1. Implement pairing check
2. Test with example circuit

## Phase 4: Optimization & Polish

### Performance improvements
- Consider Montgomery form for field arithmetic
- Add benchmarks using BenchmarkTools.jl
- Profile hot paths

### Code quality
- Remove or archive unused code
- Consistent naming conventions
- Add documentation
- Clean up exports

## Technical Notes

### Why we switched from UInt256 to BigInt
- Julia requires type parameters to be `isbits` types
- BigInt is not `isbits` so can't be used as type parameter
- UInt256 arithmetic caused LLVM compilation failures
- Solution: Use concrete types with BigInt internally

### Key files to remember
- `FiniteFields.jl` - Active field implementation
- `examples/test_r1cs_qap.jl` - Working example
- `FieldElement.jl` - Original design (kept for reference)

### Dependencies added
- BitIntegers.jl - For UInt256 type (still used in some places)
- StaticArrays.jl - For efficient small arrays
- BenchmarkTools.jl - For performance testing

## Next Steps

1. Run and fix the new test suite (`julia test/runtests_new.jl`)
2. Enable BN254Curve.jl with new field types
3. Complete the Groth16 proof system
4. Create full end-to-end example

## Command to continue

When returning to this project, start by:
```bash
cd GrothAlgebra
julia test/runtests_new.jl
```

Fix any test failures, then proceed with Phase 2.