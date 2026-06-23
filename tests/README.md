# Reproducing the bug / testing the fix

The bug only shows up for a **large/elongated FFT mesh in parallel**. A convenient trigger is a long device cell whose mesh is something like `32 x 180 x 2430` at `MeshCutoff 400 Ry`. Use your own structure (any cell whose mesh has a large dimension); warm-start from a converged DM so iteration 1 is sane.

## fdf (reciprocal-space charge mixing)

```fdf
SystemName cm_test
SystemLabel CM
%include STRUCT_DEVICE.fdf        # <-- your structure + species/psml

SolutionMethod diagon
XC.functional GGA
XC.authors    PBE
PAO.BasisSize DZP
PAO.EnergyShift 0.02 Ry
MeshCutoff 400.0 Ry
ElectronicTemperature 1000 K
OccupationFunction MP
OccupationMPOrder 2

# --- reciprocal-space (Kerker) charge mixing: the path this fix repairs ---
SCF.Mix charge
DM.MixingWeight 0.10
SCF.Kerker.q0sq 1.0 Ry
SCF.RhoG.DIIS.Depth 8
SCF.RhoG.DIIS.UseSVD False         # GitLab #82 workaround (separate known bug)

SCF.DM.Tolerance 1.d-4
MaxSCFIterations 50
DM.UseSaveDM T                     # warm start
MD.NumCGSteps 0
```

## Run serial vs parallel and compare

```bash
SIESTA=~/siesta/build/Src/siesta       # your patched binary
for np in 1 16; do
  echo "=== np=$np ==="
  OMP_NUM_THREADS=1 mpirun -np $np $SIESTA cm_test.fdf > cm_np$np.log 2>&1
  grep 'scf:' cm_np$np.log | head -4
done
```

## Expected behaviour

**Stock SIESTA 5.4.2, parallel (np>1):**
```
scf: 1   -142776.6
scf: 2  +81,972,021      <-- explodes (garbage density from the broken parallel FFT)
```

**Patched, parallel (any np) — reproduces serial:**
```
scf: 1   -142776.6
scf: 2     +235,969      <-- sane; bit-for-bit equal to the np=1 run
scf: 3      +28,873
scf: 4      -38,160      <-- recovering; converges normally
```

The serial (`np=1`) run is correct with and without the patch. **The test passes if the
`np=16` scf:2 value matches the `np=1` value** (to ~6 sig figs) instead of exploding to ~1e8.

(Energies above are from the development system — a 210-atom graphene–Cr device. Your
absolute numbers will differ; what matters is parallel == serial, no blow-up.)
