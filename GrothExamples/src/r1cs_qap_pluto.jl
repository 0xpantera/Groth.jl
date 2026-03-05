### A Pluto.jl notebook ###
# v0.20.23

using Markdown
using InteractiveUtils

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300001
begin
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."))
    Pkg.instantiate()

    using PlutoUI
    using AbstractAlgebra

    TableOfContents()
end

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300002
md"""
# Toy R1CS -> QAP over 𝔽₁₇ (end-to-end, with math + code)

## by: francos.eth

This notebook computes a complete toy pipeline:

1. Define a tiny **R1CS** instance (matrices $A,B,C$) and a witness $w$
2. Convert it to a **QAP** by interpolating columns into polynomials $A_j(x),B_j(x),C_j(x)$
3. Form the witness-weighted polynomials $A(x),B(x),C(x)$
4. Define the **constraint polynomial**
   $$P(x) = A(x)B(x) - C(x)$$
5. Define the **vanishing polynomial**
   $$t(x) = \prod_{i=1}^m (x - i)$$
6. Verify the key QAP identity:
   $$\exists H(x)\ \text{s.t.}\ P(x) = H(x)\,t(x)$$

We work in the finite field $\mathbb{F}_{17}$ and use **AbstractAlgebra.jl** for:
- finite fields ($\mathrm{GF}(p)$),
- polynomial rings $\mathbb{F}_{p}[x]$,
- interpolation,
- evaluation,
- polynomial division with remainder (`divrem`).

---
"""

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300003
md"""
## Shared environment mode

This notebook uses:

- `Pkg.activate(joinpath(@__DIR__, ".."))`
- `Pkg.instantiate()`

so Pluto runs against the repository's shared `GrothExamples` environment.
"""

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300004
md"""
## Field and polynomial ring

We fix:
- a prime $p = 17$,
- the finite field $\mathbb{F}_{17}$,
- the polynomial ring $R = \mathbb{F}_{17}[x]$.

In AbstractAlgebra:
- `F = GF(p)` is the field $\mathbb{F}_p$
- `R, x = polynomial_ring(F, :x)` constructs the polynomial ring and returns the **formal indeterminate** `x`.

**Important:** `x` here is **not a number**. It is the polynomial variable living in $R$.
"""

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300005
begin
    p = 17
    F = GF(p)
    R, x = polynomial_ring(F, :x)
    (p=p, F=F, R=R, x=x)
end

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300006
md"""
## Toy witness layout

We use a witness vector with a constant term:
$$w = (1,\ X,\ Y,\ Z).$$

The leading `1` is a common convention in R1CS so constraints can include constants using linear combinations.

Pick a satisfying assignment in $\mathbb{F}_{17}$:

- `X = 3`
- `Y = 5` (since `3 + 2 = 5`)
- `Z = 15` (since `3 * 5 = 15`)
"""

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300007
w = [F(1), F(3), F(5), F(15)]

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300008
md"""
## Toy constraints

We encode two constraints ($m=2$):

1. $$X \cdot Y = Z$$
2. $$(X + 2)\cdot 1 = Y$$
"""

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300009
begin
    A = Matrix{elem_type(F)}(undef, 2, 4)
    B = Matrix{elem_type(F)}(undef, 2, 4)
    C = Matrix{elem_type(F)}(undef, 2, 4)

    A[1, :] = [F(0), F(1), F(0), F(0)]
    B[1, :] = [F(0), F(0), F(1), F(0)]
    C[1, :] = [F(0), F(0), F(0), F(1)]

    A[2, :] = [F(2), F(1), F(0), F(0)]
    B[2, :] = [F(1), F(0), F(0), F(0)]
    C[2, :] = [F(0), F(0), F(1), F(0)]

    (A=A, B=B, C=C, m=2, n=4)
end

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300010
md"""
### R1CS sanity check

Compute:
$$(Aw)_i,\ (Bw)_i,\ (Cw)_i$$
and verify:
$$(Aw)_i\cdot (Bw)_i = (Cw)_i$$
for each constraint $i$.
"""

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300011
rowdot(M, i, wv) = sum(M[i, j] * wv[j] for j in 1:length(wv))

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300012
begin
    checks = map(1:2) do i
        ai = rowdot(A, i, w)
        bi = rowdot(B, i, w)
        ci = rowdot(C, i, w)
        (
            constraint=i,
            Aw=string(ai),
            Bw=string(bi),
            Cw=string(ci),
            lhs=string(ai * bi),
            ok=ai * bi == ci,
        )
    end
    checks
end

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300013
md"""
## QAP idea: interpolate columns into polynomials

Constraint indices are:
$x \in \{1,2,\dots,m\}.$

For each witness coordinate $j\in\{1,\dots,n\}$, define polynomials:
- $A_j(x)$ such that $A_j(i) = A_{ij}$ for all constraints $i$,
- similarly $B_j(x), C_j(x)$.

Then define witness-weighted sums:
$
A(x) = \sum_{j=1}^n w_j A_j(x),\quad
B(x) = \sum_{j=1}^n w_j B_j(x),\quad
C(x) = \sum_{j=1}^n w_j C_j(x).
$

By linearity of interpolation:
$A(i) = (Aw)_i,\ \ B(i) = (Bw)_i,\ \ C(i) = (Cw)_i.$
"""

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300014
xs = [F(1), F(2)]

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300015
function lagrange_interp(xs, ys, R)
    @assert length(xs) == length(ys)
    T = base_ring(R)
    acc = zero(R)
    for i in eachindex(xs)
        num = one(R)
        den = one(T)
        for j in eachindex(xs)
            j == i && continue
            num *= R([-xs[j], one(T)])
            den *= xs[i] - xs[j]
        end
        acc += ys[i] * inv(den) * num
    end
    acc
end

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300016
begin
    nvars = size(A, 2)
    Acols = [lagrange_interp(xs, [A[i, j] for i in 1:2], R) for j in 1:nvars]
    Bcols = [lagrange_interp(xs, [B[i, j] for i in 1:2], R) for j in 1:nvars]
    Ccols = [lagrange_interp(xs, [C[i, j] for i in 1:2], R) for j in 1:nvars]
    (Acols=Acols, Bcols=Bcols, Ccols=Ccols)
end

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300017
md"""
## Combine with witness: build A(x), B(x), C(x)

We compute:
$$A(x) = \sum_j w_j A_j(x)$$
and similarly for $B(x), C(x)$.

Then we verify:
$$A(i) = (Aw)_i,\quad B(i) = (Bw)_i,\quad C(i) = (Cw)_i.$$
"""

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300018
begin
    Aofx = sum((w[j] * Acols[j] for j in eachindex(w)), init=zero(R))
    Bofx = sum((w[j] * Bcols[j] for j in eachindex(w)), init=zero(R))
    Cofx = sum((w[j] * Ccols[j] for j in eachindex(w)), init=zero(R))
    (Aofx=Aofx, Bofx=Bofx, Cofx=Cofx)
end

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300019
md"""
## Constraint polynomial P(x)

Define:
$$P(x) = A(x)B(x) - C(x).$$

At a constraint index $i$:
$$P(i) = A(i)B(i) - C(i).$$

So:
- if constraint $i$ holds, then $P(i)=0$,
- if constraint $i$ fails, then $P(i)\ne 0$.
"""

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300020
Pofx = Aofx * Bofx - Cofx

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300021
begin
    [
        (x=string(t), P=string(evaluate(Pofx, t)))
        for t in xs
    ]
end

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300022
md"""
## Vanishing polynomial t(x)

Define the **vanishing polynomial** over constraint indices:
$$t(x) = \prod_{i=1}^m (x - i).$$

It satisfies:
$$t(i)=0 \quad \forall i\in\{1,\dots,m\}.$$

Key fact:
$$P(1)=\cdots=P(m)=0 \quad\Longleftrightarrow\quad t(x)\mid P(x).$$

So if all constraints hold, there exists a polynomial $H(x)$ such that:
$$P(x) = H(x)\,t(x).$$

In practice, we compute $H$ by polynomial division:
- `divrem(P, t)` gives quotient `H` and remainder `r`,
- correctness means `r == 0`.
"""

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300023
begin
    tofx = prod((x - xi for xi in xs), init=one(R))
    q, r = divrem(Pofx, tofx)
    (tofx=tofx, quotient=q, remainder=r, divisible=iszero(r))
end

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300024
md"""
## Bad witness demo

Now we deliberately break the witness so a constraint fails.

Then we recompute $P_{\text{bad}}(x)$ and check:
- values at constraint indices are *not all zero*,
- `divrem(P_bad, t)` produces a nonzero remainder.
"""

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300025
begin
    w_bad = [F(1), F(3), F(6), F(15)]

    A_bad = sum((w_bad[j] * Acols[j] for j in eachindex(w_bad)), init=zero(R))
    B_bad = sum((w_bad[j] * Bcols[j] for j in eachindex(w_bad)), init=zero(R))
    C_bad = sum((w_bad[j] * Ccols[j] for j in eachindex(w_bad)), init=zero(R))
    P_bad = A_bad * B_bad - C_bad

    q_bad, r_bad = divrem(P_bad, tofx)

    (
        witness=string.(w_bad),
        point_values=[(x=string(t), P=string(evaluate(P_bad, t))) for t in xs],
        divisible=iszero(r_bad),
        remainder=r_bad,
        quotient=q_bad,
    )
end

# ╔═╡ 4f6c90b1-7d41-4e0a-8f7b-4bb1d1300026
md"""
## Wrap-up

We demonstrated the exact pipeline:

- **R1CS** (many constraints):
  $$(Aw)_i\cdot(Bw)_i=(Cw)_i \quad \forall i.$$

- **QAP** (one divisibility statement):
  1) interpolate columns -> $A_j(x),B_j(x),C_j(x)$
  2) combine with witness -> $A(x),B(x),C(x)$
  3) form $P(x)=A(x)B(x)-C(x)$
  4) build $t(x)=\prod_{i=1}^m(x-i)$
  5) check $t(x)\mid P(x)$ via `divrem(P,t)`.

Compare this reference notebook with `r1cs_qap_groth_pluto.jl` for the Groth.jl package-native version.
"""

# ╔═╡ Cell order:
# ╠═4f6c90b1-7d41-4e0a-8f7b-4bb1d1300001
# ╟─4f6c90b1-7d41-4e0a-8f7b-4bb1d1300002
# ╟─4f6c90b1-7d41-4e0a-8f7b-4bb1d1300003
# ╟─4f6c90b1-7d41-4e0a-8f7b-4bb1d1300004
# ╠═4f6c90b1-7d41-4e0a-8f7b-4bb1d1300005
# ╟─4f6c90b1-7d41-4e0a-8f7b-4bb1d1300006
# ╠═4f6c90b1-7d41-4e0a-8f7b-4bb1d1300007
# ╟─4f6c90b1-7d41-4e0a-8f7b-4bb1d1300008
# ╠═4f6c90b1-7d41-4e0a-8f7b-4bb1d1300009
# ╟─4f6c90b1-7d41-4e0a-8f7b-4bb1d1300010
# ╠═4f6c90b1-7d41-4e0a-8f7b-4bb1d1300011
# ╠═4f6c90b1-7d41-4e0a-8f7b-4bb1d1300012
# ╟─4f6c90b1-7d41-4e0a-8f7b-4bb1d1300013
# ╠═4f6c90b1-7d41-4e0a-8f7b-4bb1d1300014
# ╠═4f6c90b1-7d41-4e0a-8f7b-4bb1d1300015
# ╠═4f6c90b1-7d41-4e0a-8f7b-4bb1d1300016
# ╟─4f6c90b1-7d41-4e0a-8f7b-4bb1d1300017
# ╠═4f6c90b1-7d41-4e0a-8f7b-4bb1d1300018
# ╟─4f6c90b1-7d41-4e0a-8f7b-4bb1d1300019
# ╠═4f6c90b1-7d41-4e0a-8f7b-4bb1d1300020
# ╠═4f6c90b1-7d41-4e0a-8f7b-4bb1d1300021
# ╟─4f6c90b1-7d41-4e0a-8f7b-4bb1d1300022
# ╠═4f6c90b1-7d41-4e0a-8f7b-4bb1d1300023
# ╟─4f6c90b1-7d41-4e0a-8f7b-4bb1d1300024
# ╠═4f6c90b1-7d41-4e0a-8f7b-4bb1d1300025
# ╟─4f6c90b1-7d41-4e0a-8f7b-4bb1d1300026
