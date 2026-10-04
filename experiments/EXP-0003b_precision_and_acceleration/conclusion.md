# EXP-0003b: Ledger precision (float64 vs float32) and CPU/GPU acceleration

**Run:** `.venv-gpu/Scripts/python.exe experiments/EXP-0003b_precision_and_acceleration/run.py`
(about 25 min, mostly the small-batch GPU accuracy runs)
**Outputs:** `results/metrics.json`, `figures/exp0003b_precision_throughput.png`
**Hardware:** NVIDIA RTX A2000 12 GB (driver 610.88, torch 2.14.1+cu126) and a laptop CPU (single-threaded numpy)

## 1. Precision of the energy ledger

Uniform modulation, Q = 200, δ ∈ {0, 0.005, 0.02, 0.04}, 400 mode periods (3.2e5 steps).

| Precision | max \|E_net\|/scale (exact ledger) | Floquet/decay rates |
|---|---|---|
| float64 (CPU = GPU, bit-for-bit) | 1.2e-13, 2.0e-14, 1.2e-14, 4.2e-15 | match theory (same values as float32 to 4 digits) |
| float32 (all float32) | 1.1e-3, 1.0e-3, 9.5e-6, 1.4e-5 | match theory |
| float32 fields + float64 ledger accumulation | 1.2e-3, 3.8e-4, 1.0e-5, 1.1e-5 | match theory |

For comparison, the naive-ledger artifacts are 3.8e-5 to 2.7e-3 (float64).

## Interpretation

- **float64.** The exact ledger floor is ~1e-13 after 3.2e5 steps, ten orders of
  magnitude below the naive artifacts. The GPU (torch, float64) reproduces the CPU
  to round-off.
- **float32.** Field dynamics remain accurate (growth and decay rates are correct),
  but the ledger floor rises to 1e-5 to 1e-3. That is **the same level as the
  naive artifacts**, so float32 cannot tell physics from artifacts. This is recorded as the finding
  `float32_floor_below_naive_artifacts = false`, an expected negative result, not a
  pass/fail gate.
- Accumulating the ledger in float64 helps only slightly, so the floor is set by
  rounding of the float32 field updates themselves, not by the accumulator. This is
  worst for decaying runs, where the stored energy falls far below the energy scale.
- **Policy.** Audit-grade runs use float64. float32 is allowed only for exploratory
  screening, and every reported number must be re-run in float64. This matches the
  project rule "re-evaluate extreme candidates in float64".

## 2. Throughput of the batched solver (N = 128, Mcell-updates/s)

| Batch B | numpy f64 | numpy f32 | GPU f64 | GPU f32 |
|---|---|---|---|---|
| 16 | 19 | 26 | 2.1 | 2.0 |
| 256 | 41 | 71 | 34 | 32 |
| 2048 | 9.1 | 38 | 251 | 245 |
| 8192 | 8.9 | 16 | **272** | **412** |

- **Small batches (≤ 256):** the CPU wins. Each time step launches about 30 small
  kernels, so the GPU is dominated by launch latency.
- **Large batches (≥ 2048):** the GPU wins. In float64 it is **~6.6× faster than the
  best CPU case** and ~30× faster than the CPU at the same batch size. float32 is
  only ~1.5× faster than float64 on the A2000 at these sizes.
- **Policy.** Parameter maps, Floquet maps, ensembles and optimizer populations go
  to the GPU in float64 (`backend="auto"`). Single long runs stay on the CPU.

## Things tried and dropped

**CUDA Graphs** (one graph launch per step) were implemented and tested. The first
version diverged from the eager path (a 3% energy difference after one step;
bug not isolated), and at B = 4096 they gave only a 1.1× speedup, because the
float64 arithmetic rate, not launch overhead, is the bottleneck. The code was
removed rather than left in the repository in a broken state.

## Next action

- Use the GPU for EXP-0005 maps and future dataset generation (EXP-0007).
- A fused custom kernel (Triton/CUDA) could remove the small-batch launch overhead.
  This is low priority.
