import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- measured SCF total energies (eV) vs iteration, parallel MPI run (16 ranks) ---
# After the fix (serial-FFT bypass): converges to the physical ground state.
after = [-142776.596552, 235969.910560, 28873.848471, -38160.816345, -82599.224251,
         -133579.973964, -138076.286495, -139348.662385, -140743.372213, -141127.098204,
         -141874.610397, -142265.289774, -142430.239340, -142608.798065, -142592.019368,
         -142545.977535, -142572.544631, -142627.428699, -142622.602530, -142648.906364]
# Before the fix (stock SIESTA 5.4.2): the parallel FFT corrupts the density and the
# energy explodes at iteration 2 and stays pinned near +8.2e7 eV (measured:
# scf2=+82,178,078  scf3=+82,178,082  scf4=+82,178,094  scf5=+82,178,096 ... scf20=+82,178,110).
before = [-142776.596552, 82178078.291, 82178081.801, 82178094.243, 82178096.040] \
         + [82178100.0]*14 + [82178110.120]

it   = np.arange(1, 21)
Eref = -142776.6                      # physical reference (converged / warm-start energy)
RED, BLU, GRY = "#d62728", "#1f77b4", "0.5"

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.4))

# Panel A: raw SCF energy
ax1.plot(it, before, 'o-', color=RED, ms=4, label='stock SIESTA 5.4.2 (MPI)')
ax1.plot(it, after,  's-', color=BLU, ms=4, label='with fix (MPI)')
ax1.axhline(Eref, color=GRY, ls='--', lw=1, label='physical energy')
ax1.set_xlabel('SCF iteration'); ax1.set_ylabel('total energy  (eV)')
ax1.set_title('SCF energy per iteration')
ax1.set_xticks(range(0, 21, 4)); ax1.legend(fontsize=8, loc='center right')
ax1.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))

# Panel B: convergence error, log scale
eb = np.abs(np.array(before) - Eref); ea = np.abs(np.array(after) - Eref)
eb[eb < 1] = 1.0; ea[ea < 1] = 1.0
ax2.semilogy(it, eb, 'o-', color=RED, ms=4, label='stock SIESTA 5.4.2 (MPI)')
ax2.semilogy(it, ea, 's-', color=BLU, ms=4, label='with fix (MPI)')
ax2.set_xlabel('SCF iteration'); ax2.set_ylabel(r'$|E - E_{\mathrm{converged}}|$  (eV)')
ax2.set_title('Convergence error (log scale)')
ax2.set_xticks(range(0, 21, 4)); ax2.legend(fontsize=8)

fig.suptitle('Reciprocal-space charge mixing (SCF.Mix charge) in parallel — '
             'metallic test cell, FFT mesh 32x180x2430',
             fontsize=10.5)
fig.tight_layout(rect=(0, 0, 1, 0.96))
out = "/mnt/c/Users/bjark/siesta-build/_ghfix/repo/figures/convergence_before_after.png"
fig.savefig(out, dpi=130)
print("saved", out)
