import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

def load(path):
    E = []
    for l in open(path):
        s = l.split()
        if len(s) > 3 and s[0] == "scf:":
            try: E.append(float(s[3]))   # E_KS
            except: pass
    return np.array(E)

chg = load("/home/bjarke/rc210/COLDCHG.log")
Eref = -142776.6
it = np.arange(1, len(chg) + 1)
BLU, GRY = "#1f77b4", "0.5"

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.5))

# Panel A: energy descent
ax1.plot(it, chg, "s-", color=BLU, ms=4, label="cold charge mixing (Kerker + ρ(G) DIIS)")
ax1.axhline(Eref, color=GRY, ls="--", lw=1, label="converged ground state")
ax1.set_xlabel("SCF iteration"); ax1.set_ylabel("total energy  (eV)")
ax1.set_title("Cold start: SCF energy descent")
ax1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
ax1.annotate(f"atomic-density guess\n+{chg[0]/1e3:.0f}k eV", (1, chg[0]),
             xytext=(4, chg[0]*0.78), fontsize=8,
             arrowprops=dict(arrowstyle="->", color="0.4"))
ax1.legend(fontsize=8, loc="center right")

# Panel B: convergence error (log)
err = np.abs(chg - Eref); err[err < 1] = 1.0
ax2.semilogy(it, err, "s-", color=BLU, ms=4)
ax2.axhline(1.0, color=GRY, ls=":", lw=1)
ax2.set_xlabel("SCF iteration"); ax2.set_ylabel(r"$|E - E_{\mathrm{converged}}|$  (eV)")
ax2.set_title("Cold start: convergence error (log scale)")
ax2.annotate("energy converged\n(< 1 eV) by ~iter 20", (20, 1.0),
             xytext=(8, 50), fontsize=8,
             arrowprops=dict(arrowstyle="->", color="0.4"))

fig.suptitle("Cold start from the atomic-density guess (no prior density matrix) — "
             "fixed parallel charge mixing, 16 MPI ranks\n"
             "209-atom graphene–Cr cell, FFT mesh 32 × 180 × 2430", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.92))
out = "/mnt/c/Users/bjark/siesta-build/_ghfix/repo/figures/cold_start_convergence.png"
fig.savefig(out, dpi=135)
print("saved", out, "| npts:", len(chg), "| start:", chg[0], "| end:", chg[-1])
