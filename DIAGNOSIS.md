# Diagnosis: why `SCF.Mix charge` explodes in parallel, and how we found it

This is the full reasoning chain, so anyone hitting the same symptom can follow the evidence (and avoid the same dead ends).

## Background

SIESTA can mix the self-consistent density either in the density matrix (default, Pulay) or in **reciprocal space** (`SCF.Mix charge`): the real-space density is FFT'd to ρ(G), mixed there (with optional **Kerker** preconditioning `α = w·G²/(G²+q0²)` and a ρ(G) DIIS), and FFT'd back. Kerker charge mixing is the textbook cure for **charge sloshing** in large metallic cells, so for big graphene–metal contact devices it is the path that should work where Pulay stalls.

The machinery: FFT in `rhofft.F` → mixing in `m_rhog.F90` (+ `m_diis.F90`), with array setup in `dhscf.F`.

## Symptom

Large device cell (210 atoms), FFT mesh **32 × 180 × 2430**, warm-started from the converged DM, run with MPI:

```
scf: 1   -142776.6
scf: 2  +81,972,021     <- explodes within one iteration
```

Serial (`-np 1`), identical input: scf:2 = +235,969 then recovers and converges. Small cells (e.g. mesh 32×180×720) converge fine in parallel too. So it is **parallel + large-mesh specific**.

## Step 1 — rule out the mixer (it's innocent)

A controlled 2×2 matrix on rc210 (np16, warm DM):

| DIIS | Kerker (q0²) | scf:2 |
|------|--------------|-------|
| on   | 1.0          | 81972021.677380 |
| off  | 1.0          | 81972021.677380 |
| off  | 0.0          | 81972021.677380 |

**Bit-for-bit identical.** If DIIS or Kerker damping were responsible, turning them off would move that number. They don't. The mixing math is a bystander. (Reason: with `SCF.MixCharge.SCF1=false`, the first iteration's mix is a plain copy `rhog_in = rhog`; scf:2 is built from an *unmixed* density, so the only nontrivial machinery between a sane scf:1 and the exploded scf:2 is the FFT.)

## Step 2 — isolate the FFT round-trip

Instrument `rhofft` to compute `max| FFT⁻¹(FFT(ρ)) − ρ |` on the (known-good) input density:

| case               | round-trip error |
|--------------------|------------------|
| rc070 (720), np4   | 1.6e-14 (exact)  |
| rc210 (2430), np1  | 2.8e-14 (exact)  |
| rc210 (2430), np2+ | finite, ~1–2% of max ρ, and the SCF then diverges |

So the **1-D FFTs are exact** (serial is exact for the very same large mesh), and the corruption appears **only when the mesh is split across MPI ranks**. The error magnitude **grows with the amount of pencil decomposition** (np4: ProcessorY=4/Z=1 ≈ 1.5e6 blow-up; Y=1/Z=4 ≈ 9e6), and serial (no split) is exact. This is the fingerprint of a **distributed-transpose** problem, not a 1-D-FFT or mixing problem.

## Step 3 — narrow it down (and what does NOT fix it)

- **It is the shared transpose.** SIESTA's parallel 3-D FFT (`m_fft::fft`) does pencil transposes via `redistributePencil`. With FFTW (`SIESTA_WITH_FFTW=ON`), only the 1-D kernels change — `redistributePencil` is still called unconditionally for the transposes. So **an FFTW build fails identically** (verified by reading the code paths; the transpose calls are outside every `#ifdef SIESTA__FFTW`). Switching to FFTW is not a fix.
- **`ProcessorY` workaround (GitLab #170) does not help here.** Sweeping `ProcessorY` ∈ {1,2,4} and `FFT.Processory.traditional` all still blow up (only the magnitude changes). Any `np>1` forces at least one pencil split.
- **It is not a buffer overrun** at np2 (sizes check out), nor the `ordix` "serial form only" sort (that only feeds an optional debug table), nor a global FFT-distribution-flag flip between `rhofft` (wants x-pencil) and `poison` (z-pencil) — that flip is a *separate* original-code issue that a consistent distribution choice neutralizes, but a residual round-trip error remains even back-to-back.
- A rigorous static + dynamic trace of `redistributePencil` found the *permutation* to be loss-less for these splits, and the only concrete coding defects (two `gpfa` slice typos `f(IOffSet+1:n)`→`:2*n`, and `gpfa` dummies changed to assumed-shape `A(:)`) did **not**, on their own, resolve the parallel blow-up. So the residual corruption is subtle and mesh/decomposition-dependent — exactly the under-determined failure mode documented upstream.

## Step 4 — upstream context

This matches **GitLab #170** ("Grid initialization crashes … for certain numbers of MPI tasks … the distribution of grid points over the tasks … manually setting `ProcessorY` sometimes fixes the problem … stack traces occur in the `fft` routines"; parallel-only). Related: **#169** (Y-Z process-grid ordering / non-contiguous slabs), addressed by **MR !115** (FFT distribution rework + FFTW-for-Poisson + `Poisson.Method`/`FFT.Processory.traditional`). The internal FFT also documents a hard requirement (manual §Mesh.Sizes): each mesh dimension must factor into 2/3/5 and divide `Mesh.SubDivisions` (rc210's 2430 = 2·3⁵·5 does satisfy this, so that rule is not the issue here).

## The fix — bypass the broken transpose with the exact serial FFT

Since the **serial** FFT is bit-exact even for this mesh, route the charge-mixing FFT through it:

1. `fft.F`: add `fft_serial_global` — a serial 3-D complex FFT over a *global* array, reusing the existing `gpfa` 1-D routine and the trig tables that `fft_init` already builds.
2. `rhofft.F`: when `Nodes > 1`, **gather** each rank's local mesh slice to one rank (in global natural order, using the same `getMeshBox(UNIFORM)` mapping that `order_rhog` uses), run `fft_serial_global`, and **scatter** the result back into the same local layout. For `Nodes == 1`, keep the original (correct) serial path. Keep `rhofft`'s existing unit scaling.
3. `m_fft_gpfa.F`: revert `gpfa` dummies to assumed-size `A(*)` and fix the two slice typos (correctness hardening; cheap and right regardless).

`m_rhog.F90` and `dhscf.F` are **unchanged** — the bypass deliberately reproduces the exact local layout (`J = 1 + J1 + N1·J2L + N1·N2·J3L`, full X, local contiguous (y,z) block) that `order_rhog` indexes, so `g2`↔`rhog` stay consistent.

**Cost:** the bypass gathers `product(mesh)` complex numbers to one rank per FFT call (~0.2 GB for this mesh) and does one serial 3-D FFT there. The FFT is cheap relative to the diagonalization, and the gather/scatter is a small fraction of an SCF step, so the overall hit is negligible — the diagonalization stays fully parallel.

## Result

rc210 (mesh 32×180×2430), warm DM, **np16**, `SCF.Mix charge` + Kerker:

```
scf: 1   -142776.6
scf: 2   +235,969     <- now SANE (was +82,000,000); identical to the serial run
scf: 3    +28,873
scf: 4    -38,160
 ...      converges toward the ground state, exactly like serial
```

rc070 regression (mesh 32×180×720, np16): converges to the Pulay baseline −120534.5.

The parallel charge-mixing path now produces the serial answer.
