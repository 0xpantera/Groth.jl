using Random
using GrothAlgebra: BN254Fr, bn254_fr

"""
    GeneratedCircuit

Simple container for randomly generated R1CS data.

- `r1cs`: the constructed `R1CS` instance
- `witness`: satisfying `Witness`
- `public_indices`: indices (1-based) of public inputs
- `description`: brief string describing constraints for debugging
"""
struct GeneratedCircuit{F}
    r1cs::R1CS{F}
    witness::Witness{F}
    public_indices::Vector{Int}
    used_public::Vector{Int}
    description::Vector{String}
end

"""
    generate_small_r1cs(rng; max_constraints=4, max_public=5, value_range=-20:20, retries=20)

Generate a compact R1CS/witness pair using basic gate templates.

The construction keeps circuits intentionally small:
- number of constraints ≤ `max_constraints`
- variables are introduced gradually, each constraint either multiplies two existing
  values or forms an affine combination.
- public inputs include the constant 1 and the first `public_count` variable slots.

Returns a `GeneratedCircuit` or throws if unable to find a valid witness within `retries` attempts.
"""
function generate_small_r1cs(rng::AbstractRNG; max_constraints::Int=4, max_public::Int=5, value_range=-20:20, retries::Int=20)
    max_constraints >= 2 || error("expect at least two constraints")
    for attempt in 1:retries

        n_constraints = rand(rng, 2:max_constraints)
        # Start with 1 (constant) plus r variable and at least one auxiliary slot
        num_vars = 1 + 1 + n_constraints + 2
        num_public = min(max_public, num_vars)

        F = BN254Fr
        L = zeros(F, n_constraints, num_vars)
        R = zeros(F, n_constraints, num_vars)
        O = zeros(F, n_constraints, num_vars)

        values = zeros(Int, num_vars)
        values[1] = 1

        idx_r = 2
        descriptions = String[]

        base_values = [rand(rng, value_range) for _ in 1:num_vars]
        for i in 2:num_public
            values[i] = base_values[i]
        end

        next_var = num_public + 1
        available = collect(3:num_public)  # exclude constant and r
        if isempty(available)
            continue
        end
        used_public = Int[]

        # Build constraints sequentially; accumulate new intermediate variables
        for c in 1:n_constraints
            if c == n_constraints || next_var > num_vars
                isempty(available) && break
                a_idx = rand(rng, available)
                b_idx = rand(rng, available)
                L[c, a_idx] = one(F)
                R[c, b_idx] = one(F)
                O[c, idx_r] = one(F)
                values[idx_r] = values[a_idx] * values[b_idx]
                push!(descriptions, "v$idx_r = v$a_idx * v$b_idx")
                if a_idx <= num_public
                    push!(used_public, a_idx)
                end
                if b_idx <= num_public
                    push!(used_public, b_idx)
                end
                break
            end

            gate_type = rand(rng, (:mul, :affine))
            if gate_type == :mul
                if isempty(available) || next_var > num_vars
                    continue
                end
                a_idx = rand(rng, available)
                b_idx = rand(rng, available)
                L[c, a_idx] = one(F)
                R[c, b_idx] = one(F)
                O[c, next_var] = one(F)
                values[next_var] = values[a_idx] * values[b_idx]
                push!(descriptions, "v$next_var = v$a_idx * v$b_idx")
                if a_idx <= num_public
                    push!(used_public, a_idx)
                end
                if b_idx <= num_public
                    push!(used_public, b_idx)
                end
            else
                if isempty(available) || next_var > num_vars
                    continue
                end
                a_idx = rand(rng, available)
                b_idx = rand(rng, available)
                α = rand(rng, value_range)
                β = rand(rng, value_range)
                γ = rand(rng, value_range)
                if α == 0 && β == 0 && γ == 0
                    γ = 1
                end
                L[c, a_idx] = bn254_fr(α)
                L[c, b_idx] += bn254_fr(β)
                L[c, 1] += bn254_fr(γ)
                R[c, 1] = one(F)
                O[c, next_var] = one(F)
                values[next_var] = α * values[a_idx] + β * values[b_idx] + γ
                push!(descriptions, "v$next_var = $α*v$a_idx + $β*v$b_idx + $γ")
                if a_idx <= num_public
                    push!(used_public, a_idx)
                end
                if b_idx <= num_public
                    push!(used_public, b_idx)
                end
            end
            push!(available, next_var)
            next_var += 1
        end

        if values[idx_r] == 0
            if length(available) >= 1
                a_idx = last(available)
                values[idx_r] = values[a_idx] * values[a_idx]
            else
                continue
            end
        end

        # Convert witness values into field elements
        witness_vals = Vector{BN254Fr}(undef, num_vars)
        for i in 1:num_vars
            witness_vals[i] = bn254_fr(values[i])
        end

        r1cs = R1CS{BN254Fr}(num_vars, n_constraints, num_public, L, R, O)
        witness = Witness{BN254Fr}(witness_vals)

        if is_satisfied(r1cs, witness)
            public_indices = collect(1:num_public)
            unique!(used_public)
            filter!(i -> 2 <= i <= num_public, used_public)
            return GeneratedCircuit(r1cs, witness, public_indices, used_public, descriptions)
        end
    end

    error("failed to generate satisfying R1CS after $(retries) attempts")
end
