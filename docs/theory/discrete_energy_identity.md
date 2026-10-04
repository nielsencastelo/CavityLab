# Discrete-exact energy accounting for the 1D Yee scheme with time-varying permittivity

This note derives the energy ledger implemented in
[`src/cavitylab/solvers/fdtd1d.py`](../../src/cavitylab/solvers/fdtd1d.py). The
lossless identity is classical (see [prior art](../literature/prior_art_preprint01.md)).
The terms for time-varying permittivity (pump work), semi-implicit conduction and
soft sources are written in the same discrete form so that the full balance closes
to round-off.

## Setting

We use the Ez/Hy polarization on `[0, L]` with PEC walls (`E_0 = E_N = 0`).
`E_i^n` and `D_i^n` live on nodes `x_i = i Δx`, and `H_{i+1/2}^{n+1/2}` lives on
half nodes. The scheme is:

```
H^{n+1/2}_{i+1/2} = H^{n-1/2}_{i+1/2} + Δt/(μ Δx) (E^n_{i+1} - E^n_i)
D^{n+1}_i        = D^n_i + Δt/Δx (H^{n+1/2}_{i+1/2} - H^{n+1/2}_{i-1/2}) - Δt σ_i Ē_i - Δt J_i^{n+1/2}
E^{n+1}_i        = D^{n+1}_i / ε_i^{n+1},           Ē_i = (E^n_i + E^{n+1}_i)/2
```

Write `G` for the forward difference `(GE)_{i+1/2} = E_{i+1} - E_i`. With PEC walls,
the backward difference used in the D update is `-Gᵀ` (summation by parts with
vanishing boundary terms).

## Discrete energy

```
W^n = Δx/2 Σ_i D^n_i E^n_i  +  Δx/2 μ Σ_j H^{n+1/2}_j H^{n-1/2}_j
```

The magnetic term pairs the two half steps that surround `t_n`. Under the CFL
condition it is a positive quadratic form, which is the basis of the classical
energy-method stability proof.

## One-step balance

**Electric part.** Expanding `D'E' - DE` (primes denote step n+1) gives

```
½(D'E' - DE) = ½(D' - D)(E' + E)  +  ½(D E' - D' E)
```

With `E = D/ε`, the second term is `½ D D' (1/ε' - 1/ε)`. It vanishes when `ε` is
constant and equals the **work done by the modulation on the field**. Its continuum
limit is `-½ E² ∂ε/∂t` per unit volume, the work needed to change `ε` at fixed `D`.

**Magnetic part.**

```
½μ H^{n+1/2}(H^{n+3/2} - H^{n-1/2}) = (Δt/Δx) H^{n+1/2} · G Ē
```

**Exchange.** The Maxwell exchange terms cancel by summation by parts:
`(Δt/Δx)[(-GᵀH)·Ē + H·GĒ] = 0`.

**Result.** The exact discrete balance is

```
W^{n+1} - W^n = Δx Σ_i [ ½ D^n_i D^{n+1}_i (1/ε_i^{n+1} - 1/ε_i^n)   (pump work, signed)
                        - Δt σ_i Ē_i²                                (dissipation ≥ 0)
                        - Δt J_i^{n+1/2} Ē_i ]                       (source work, signed)
```

It holds for any `σ(x) ≥ 0`, any `ε(x, t) > 0`, any source, and any Courant number
within the CFL limit. Dissipation in nodes flagged as a port is booked as `E_out`.
The total residual is therefore pure floating-point round-off, which grows roughly
like `n_steps · ε_machine`.

## Naive estimators

Most hand-written post-processing uses continuum formulas evaluated pointwise:

| Quantity | Naive estimator | Error |
|---|---|---|
| Stored energy | `½ε(Eⁿ)² + ½μ(H̄ⁿ)²`, with `H̄ⁿ` the average of the two half steps | O((ωΔt)²), oscillating |
| Dissipation | `σ Δt (Eⁿ)²` (left point) | O(Δt) cumulative |
| Source work | `-Δt J Eⁿ` (left point) | O(Δt) cumulative |
| Pump work | `-½ (Eⁿ)² (ε^{n+1} - εⁿ)` (left point) | O(Δt) cumulative |

The residual of the naive ledger is a pure discretization artifact. It vanishes
under refinement, but at practical resolutions it is 10⁻⁵ to 10⁻² of the energy
scale and **can be positive**, which reads as apparent net energy generation.
EXP-0003, EXP-0004 and EXP-0005 quantify this.

## Sign convention for E_net

The source documents write `E_net = E_out - E_in - E_loss - ΔE_stored`. Read
literally with all quantities positive, that expression does not vanish under exact
conservation. CavityLab therefore uses

```
E_net = E_out + E_loss + ΔE_stored - E_in
```

This is the energy that appeared without a booked input. Exact conservation gives
`E_net = 0`, and the null hypothesis is `H0: E_net ≤ 0` within the numerical
uncertainty. The exit criterion of the consolidated context document,
`|E_in - E_out - ΔE_stored - E_loss| / E_scale < ε`, is exactly `|E_net| / E_scale < ε`.

## Uniform modulation reduces to a Hill equation

When `ε(t) = ε_s f(t)` with `f = 1 + δ sin(Ωt + φ)` fills the whole cavity, every
mode decouples. In both the continuum and on the grid, the spatial operator
commutes with the modulation. In dimensionless form (`τ = ω_m t`, `p = d/ε_s`,
`q = η_s h`, `Q = ω_m ε_s/σ`):

```
dp/dτ = -q - p/(Q f),      dq/dτ = p/f,      w = p²/f + q²
```

This is a damped Hill equation:

- Near `Ω = 2ω_m` the amplitude Floquet exponent is `δ/4 - 1/(2Q)`, so the
  threshold is `δ_th = 2/Q`.
- The other branch decays as `-δ/4 - 1/(2Q)`.
- A start in a single field quadrature with pump phase 0 lies exactly on the
  decaying branch. That is a degenerate case, so time-series comparisons use a
  generic pump phase.

[`cavitylab.benchmarks.parametric`](../../src/cavitylab/benchmarks/parametric.py)
measures the monodromy matrix directly from FDTD to compare both branches.
