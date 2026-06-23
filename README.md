# SIESTA 5.4.2 — parallel charge-density mixing (`SCF.Mix charge`) FFT fix

**TL;DR** — In SIESTA 5.4.2, reciprocal-space charge-density mixing (`SCF.Mix charge`, a.k.a. Kerker / ρ(G) mixing) **silently corrupts the density in parallel MPI runs on large meshes**, making the SCF energy explode to ~10⁷–10⁸ eV within a couple of iterations. The same calculation is **correct in serial** and correct for small meshes. The root cause is **not** the mixing math — it is SIESTA's distributed 3-D FFT (the `redistributePencil` pencil-transpose in `m_fft`), which does not round-trip correctly for large/elongated meshes in parallel (cf. upstream GitLab issues [#170](https://gitlab.com/siesta-project/siesta/-/issues/170) / [#169](https://gitlab.com/siesta-project/siesta/-/issues/169)).

This repo provides a **drop-in fix**: route the (cheap) charge-mixing FFT through a **gather → serial FFT → scatter** path. The serial FFT is bit-exact, so the parallel run reproduces the serial result. The expensive parts of the SCF loop (diagonalization) stay fully parallel.

> Bug found and fixed during a DFT/NEGF graphene-Cr contact study, where the large device cells (200+ atoms, elongated meshes) need Kerker charge mixing to tame metallic charge sloshing — and that path was unusable in parallel.

---

## Symptom

With `SCF.Mix charge` in an MPI run on a large mesh (here: a 210-atom cell, FFT mesh **32 × 180 × 2430**):

```
scf: 1   -142776.6      (correct, warm-started)
scf: 2  +81,972,021     <-- explodes
scf: 3  +81,972,021
 ... diverges / SIGSEGV
```

The *identical* run in serial (`-np 1`) is fine:

```
scf: 1   -142776.6
scf: 2   +235,969       (normal transition excursion)
scf: 3    +28,873
scf: 4    -38,160       <-- recovering, converges normally
```

## Root cause (short)

- The blow-up is **independent of the mixing scheme**: turning DIIS on/off and Kerker damping on/off (q0²=1 vs 0) gives a **bit-identical** scf:2 energy. So the mixer is a bystander.
- It is the **FFT density reconstruction**. A direct round-trip test `FFT⁻¹(FFT(ρ)) − ρ` gives **~1e-14 in serial and for small meshes, but a finite error in parallel for the large mesh** — and the error grows with the number of MPI pencil splits.
- It lives in the **shared parallel transpose** `redistributePencil` (used by *both* the internal gpfa FFT and a FFTW build — FFTW only replaces the 1-D kernels, not the transpose), so switching to FFTW does **not** help.
- This is the failure class of upstream **GitLab #170** ("Grid initialization … for certain numbers of MPI tasks … occurs in the fft routines … setting `ProcessorY` sometimes fixes it"). For this mesh, no `ProcessorY` choice avoids it.

Full story and evidence: **[DIAGNOSIS.md](DIAGNOSIS.md)**.

## The fix

Because the **serial** 3-D FFT is exact, the charge-mixing FFT is rerouted through it when running in parallel:

- **`rhofft.F`** — when `Nodes > 1`, gather the mesh to one rank, run the exact serial 3-D FFT there, scatter back (instead of the buggy distributed `fft`). Layout is kept identical to what `order_rhog` (in `m_rhog.F90`) already indexes, so no other file needs to change. The FFT is cheap relative to diagonalization, so the gather/scatter overhead is negligible.
- **`fft.F`** — adds `fft_serial_global`, a serial 3-D complex FFT over a global array (reuses the existing `gpfa` 1-D routine and trig tables). Also fixes two genuine `gpfa` slice typos (`f(IOffSet+1:n)` → `:2*n`).
- **`m_fft_gpfa.F`** — reverts the `gpfa` array dummies from assumed-shape `A(:)` back to assumed-size `A(*)` (the strided index arithmetic assumes assumed-size).

See `patches/` for the diffs against pristine SIESTA 5.4.2, and `src/` for the full patched files.

## Validation

Built with foss-2023a-style toolchain (gfortran 13, OpenMPI, ScaLAPACK/ELPA). All runs warm-started from the converged DM; `MeshCutoff 400 Ry`, mesh 32×180×2430 (rc210) / 32×180×720 (rc070).

| case  | mesh           | np  | result |
|-------|----------------|-----|--------|
| rc210 | 32×180×2430    | 1   | converges (reference) |
| rc210 | 32×180×2430    | 16  | **converges** — scf:2 = +235,969 (was +8.2e7), trajectory matches serial bit-for-bit |
| rc070 | 32×180×720     | 16  | converges to −120534.5 (Pulay baseline) ✓ |

## How to apply

```bash
# from the root of a pristine SIESTA 5.4.2 source tree
patch -p1 < patches/rhofft.F.patch
patch -p1 < patches/fft.F.patch
patch -p1 < patches/m_fft_gpfa.F.patch
# then rebuild (cmake/ninja or your usual flow)
```
or just copy the three files from `src/` over `Src/`.

## How to reproduce / test

See `tests/` — an `fdf` with `SCF.Mix charge` + Kerker, a run script, and the expected serial energy trajectory. Run it at `-np 1` and `-np 16`; before the fix the parallel run explodes at scf:2, after the fix it matches serial.

## Status / scope

- Tested for spin-unpolarized (`nspin = 1`) Γ-and-k diagonalization runs.
- The fix is a **workaround at the application layer** (it sidesteps the distributed-FFT bug rather than fixing the transpose). The underlying `redistributePencil` issue is upstream's to fix; this makes `SCF.Mix charge` usable in parallel today.

## References

- SIESTA GitLab #170 — parallel FFT grid-distribution bug: https://gitlab.com/siesta-project/siesta/-/issues/170
- SIESTA GitLab #169 — Y/Z process-grid ordering: https://gitlab.com/siesta-project/siesta/-/issues/169
- SIESTA GitLab MR !115 — FFT distribution rework + FFTW for Poisson
- SIESTA GitLab #82 — `SCF.RhoG.DIIS.UseSVD` (a separate, already-known charge-mixing bug; use `SCF.RhoG.DIIS.UseSVD False`)

## License

The patched SIESTA source files are GPL (© the SIESTA group), same as SIESTA. The patches, docs, and test scripts in this repo are released under the same terms.
