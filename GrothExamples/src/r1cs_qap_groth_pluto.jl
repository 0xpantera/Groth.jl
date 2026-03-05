### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0001
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
    Pkg.instantiate()

    using PlutoUI
    using GrothAlgebra
    using GrothProofs

    TableOfContents()
end

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0002
md"""
# Toy R1CS -> QAP using Groth.jl packages

## by: francos.eth

This notebook mirrors the toy R1CS->QAP walkthrough, but uses this repository's packages:

- `GrothProofs` for R1CS, witness checks, and QAP conversion
- `GrothAlgebra` for polynomial operations and field arithmetic

We demonstrate:

1. Build an R1CS and satisfying witness
2. Convert to QAP polynomials
3. Form `A(x), B(x), C(x)` and the constraint polynomial `P(x) = A(x)B(x) - C(x)`
4. Check divisibility by the vanishing polynomial `t(x)`
5. Show failure for a bad witness
"""

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0003
md"""
## Shared environment mode

This notebook uses Pluto's "shared environment" pattern by calling:

- `Pkg.activate(joinpath(@__DIR__, ".."))`
- `Pkg.instantiate()`

So dependencies come from the repository environment, not Pluto's per-notebook manager.
"""

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0004
md"""
## Build the toy circuit and witness

We reuse the GrothProofs multiplication example:

- Variables: `[1, r, x, y, z, u, v1, v2]`
- Constraints:
  - `v1 = x * y`
  - `v2 = z * u`
  - `r = v1 * v2`
"""

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0005
begin
    x_val, y_val, z_val, u_val = 3, 5, 7, 11

    r1cs = create_r1cs_example_multiplication()
    witness = create_witness_multiplication(x_val, y_val, z_val, u_val)

    (r1cs=r1cs, witness=witness, satisfied=is_satisfied(r1cs, witness))
end

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0006
md"""
## R1CS matrices
"""

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0007
begin
    labels = ["1", "r", "x", "y", "z", "u", "v1", "v2"]

    L = [
		r1cs.L[i, j].value for i in 1:r1cs.num_constraints, j in 1:r1cs.num_vars
	]
    
	R = [
		r1cs.R[i, j].value for i in 1:r1cs.num_constraints, j in 1:r1cs.num_vars
	]
    
	O = [
		r1cs.O[i, j].value for i in 1:r1cs.num_constraints, j in 1:r1cs.num_vars
	]

    (labels=labels, L, R, O)
end

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0008
md"""
### R1CS sanity check

For each row `i`, verify `(Lw)_i * (Rw)_i == (Ow)_i`.
"""

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0009
rowdot(M, i, w) = sum(M[i, j] * w[j] for j in 1:length(w))

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0010
begin
    w = witness.values
    checks = map(1:r1cs.num_constraints) do i
        li = rowdot(r1cs.L, i, w)
        ri = rowdot(r1cs.R, i, w)
        oi = rowdot(r1cs.O, i, w)
        (
            constraint=i,
            Lw=BigInt(li.value),
            Rw=BigInt(ri.value),
            Ow=BigInt(oi.value),
            ok=(li * ri == oi),
        )
    end
    checks
end

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0011
md"""
## Convert R1CS -> QAP
"""

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0012
qap = r1cs_to_qap(r1cs)

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0013
begin
    (
        num_constraints=qap.num_constraints,
        num_vars=qap.num_vars,
        points=[BigInt(p.value) for p in qap.points],
        domain_size=qap.domain.size,
        coset_size=qap.coset_domain.size,
        coset_offset=BigInt(coset_offset(qap.coset_domain).value),
        target_degree=degree(qap.t),
    )
end

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0014
md"""
### Column-wise interpolation spot check

`u[j], v[j], w[j]` are interpolated column polynomials. Evaluating at constraint points recovers matrix entries.
"""

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0015
begin
    # j = 3 corresponds to variable x in [1, r, x, y, z, u, v1, v2]
    j = 3
    [
        (
            point=BigInt(pt.value),
            u_eval=BigInt(evaluate(qap.u[j], pt).value),
            u_matrix=BigInt(r1cs.L[i, j].value),
            v_eval=BigInt(evaluate(qap.v[j], pt).value),
            v_matrix=BigInt(r1cs.R[i, j].value),
            w_eval=BigInt(evaluate(qap.w[j], pt).value),
            w_matrix=BigInt(r1cs.O[i, j].value),
        )
        for (i, pt) in enumerate(qap.points)
    ]
end

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0016
md"""
## Build `A(x)`, `B(x)`, `C(x)` from witness
"""

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0017
begin
    F = BN254Fr
    A_poly = sum((witness.values[i] * qap.u[i] for i in 1:qap.num_vars), init=zero(Polynomial{F}))
    B_poly = sum((witness.values[i] * qap.v[i] for i in 1:qap.num_vars), init=zero(Polynomial{F}))
    C_poly = sum((witness.values[i] * qap.w[i] for i in 1:qap.num_vars), init=zero(Polynomial{F}))

    (
        A_degree=degree(A_poly),
        B_degree=degree(B_poly),
        C_degree=degree(C_poly),
        A_poly=A_poly,
        B_poly=B_poly,
        C_poly=C_poly,
    )
end

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0018
md"""
## Constraint polynomial `P(x)`

For a valid witness, `P(x) = A(x)B(x) - C(x)` vanishes on all constraint points.
"""

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0019
P_poly = fft_polynomial_multiply(A_poly, B_poly) - C_poly

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0020
begin
    [(
    	point = pt.value |> BigInt,
		u = evaluate_qap(qap, witness, pt) |> first |> x -> BigInt(x.value),
		v = evaluate_qap(qap, witness, pt) |> x -> x[2] |> x -> BigInt(x.value),
		w = evaluate_qap(qap, witness, pt) |> last |> x -> BigInt(x.value),
		P = evaluate(P_poly, pt) |> x -> BigInt(x.value),
    )
    for pt in qap.points
    ]
end

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0021
md"""
## Vanishing polynomial `t(x)` and quotient `H(x)`

`H(x)` should satisfy:

`P(x) = H(x) * t(x)`
"""

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0022
begin
    H_dense = compute_h_polynomial(qap, witness; use_coset=false)
    H_coset = compute_h_polynomial(qap, witness; use_coset=true)

    (
        t_degree = degree(qap.t),
        H_dense_degree = degree(H_dense),
        H_coset_degree = degree(H_coset),
        dense_equals_coset = (H_dense == H_coset),
        divisibility_identity_holds = (P_poly == H_dense * qap.t),
    )
end

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0023
md"""
## Bad witness demo

We perturb `r` so constraints fail, then check QAP divisibility breaks.
"""

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0024
begin
    bad_values = copy(witness.values)
    bad_values[2] = bad_values[2] + one(BN254Fr)  # corrupt public output r
    bad_witness = Witness{BN254Fr}(bad_values)

    bad_A = sum((bad_witness.values[i] * qap.u[i] for i in 1:qap.num_vars), init=zero(Polynomial{BN254Fr}))
    bad_B = sum((bad_witness.values[i] * qap.v[i] for i in 1:qap.num_vars), init=zero(Polynomial{BN254Fr}))
    bad_C = sum((bad_witness.values[i] * qap.w[i] for i in 1:qap.num_vars), init=zero(Polynomial{BN254Fr}))
    bad_P = fft_polynomial_multiply(bad_A, bad_B) - bad_C

    bad_point_eval = [
        (
            point=BigInt(pt.value),
            P=BigInt(evaluate(bad_P, pt).value),
        )
        for pt in qap.points
    ]

    bad_division_result = try
        bad_H = compute_h_polynomial(qap, bad_witness; use_coset=false)
        (ok=true, degree=degree(bad_H))
    catch err
        (ok=false, error=sprint(showerror, err))
    end

    (
        is_satisfied=is_satisfied(r1cs, bad_witness),
        point_evaluations=bad_point_eval,
        division=bad_division_result,
    )
end

# ╔═╡ 7f9fca90-1111-4f2e-93ec-52b9694d0025
md"""
## Wrap-up

This notebook demonstrates the R1CS->QAP identity entirely with Groth.jl packages:

- `GrothProofs` constructs and checks constraints
- `GrothProofs.r1cs_to_qap` builds QAP polynomials
- `GrothAlgebra` evaluates and manipulates the resulting polynomials

For the valid witness, divisibility holds; for the bad witness, it breaks.
"""

# ╔═╡ Cell order:
# ╠═7f9fca90-1111-4f2e-93ec-52b9694d0001
# ╟─7f9fca90-1111-4f2e-93ec-52b9694d0002
# ╟─7f9fca90-1111-4f2e-93ec-52b9694d0003
# ╟─7f9fca90-1111-4f2e-93ec-52b9694d0004
# ╠═7f9fca90-1111-4f2e-93ec-52b9694d0005
# ╟─7f9fca90-1111-4f2e-93ec-52b9694d0006
# ╠═7f9fca90-1111-4f2e-93ec-52b9694d0007
# ╟─7f9fca90-1111-4f2e-93ec-52b9694d0008
# ╠═7f9fca90-1111-4f2e-93ec-52b9694d0009
# ╠═7f9fca90-1111-4f2e-93ec-52b9694d0010
# ╟─7f9fca90-1111-4f2e-93ec-52b9694d0011
# ╠═7f9fca90-1111-4f2e-93ec-52b9694d0012
# ╠═7f9fca90-1111-4f2e-93ec-52b9694d0013
# ╟─7f9fca90-1111-4f2e-93ec-52b9694d0014
# ╠═7f9fca90-1111-4f2e-93ec-52b9694d0015
# ╟─7f9fca90-1111-4f2e-93ec-52b9694d0016
# ╠═7f9fca90-1111-4f2e-93ec-52b9694d0017
# ╟─7f9fca90-1111-4f2e-93ec-52b9694d0018
# ╠═7f9fca90-1111-4f2e-93ec-52b9694d0019
# ╠═7f9fca90-1111-4f2e-93ec-52b9694d0020
# ╟─7f9fca90-1111-4f2e-93ec-52b9694d0021
# ╠═7f9fca90-1111-4f2e-93ec-52b9694d0022
# ╟─7f9fca90-1111-4f2e-93ec-52b9694d0023
# ╠═7f9fca90-1111-4f2e-93ec-52b9694d0024
# ╟─7f9fca90-1111-4f2e-93ec-52b9694d0025
