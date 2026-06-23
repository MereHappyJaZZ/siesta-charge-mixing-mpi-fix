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

fix = load("/home/bjarke/rc210/COLDCHG.log")     # fixed parallel charge mixing (converges)
stk = load("/home/bjarke/rc210/COLDSTOCK.log")   # stock parallel, identical input (explodes)
Eref = -142776.6
RED, BLU, GRY = "#d62728", "#1f77b4", "0.5"

print("stock pts:", len(stk), "vals:", stk[:4])
print("fix   pts:", len(fix), "start/end:", fix[0], fix[-1])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.5))

# Panel A: SCF energy (symlog to span the explosion and the descent)
ax1.plot(np.arange(1, len(stk)+1), stk, "o-", color=RED, ms=4, label="stock SIESTA 5.4.2 (parallel)")
ax1.plot(np.arange(1, len(fix)+1), fix, "s-", color=BLU, ms=4, label="with fix (parallel)")
ax1.axhline(Eref, color=GRY, ls="--", lw=1, label="ground state (−142,777 eV)")
lt = max(3e5, np.abs(stk).max()/1e3 if len(stk) else 3e5)
ax1.set_yscale("symlog", linthresh=3e5)
ax1.set_xlabel("SCF iteration"); ax1.set_ylabel("total energy  (eV)")
ax1.set_title("Cold start: SCF energy per iteration")
ax1.legend(fontsize=8, loc="center right")

# Panel B: convergence error, log scale (the clean contrast)
ef = np.abs(fix - Eref); ef[ef < 1] = 1.0
es = np.abs(stk - Eref); es[es < 1] = 1.0
ax2.semilogy(np.arange(1, len(es)+1), es, "o-", color=RED, ms=4, label="stock (diverges)")
ax2.semilogy(np.arange(1, len(ef)+1), ef, "s-", color=BLU, ms=4, label="with fix (converges)")
ax2.set_xlabel("SCF iteration"); ax2.set_ylabel(r"$|E - E_{\mathrm{converged}}|$  (eV)")
ax2.set_title("Cold start: convergence error (log scale)")
ax2.legend(fontsize=8)

fig.suptitle("Reciprocal-space charge mixing (SCF.Mix charge) in parallel — COLD start from the atomic-density guess\n"
             "stock SIESTA 5.4.2 vs. with this fix · 209-atom graphene–Cr cell · 16 MPI ranks · FFT mesh 32 × 180 × 2430",
             fontsize=9.5)
fig.tight_layout(rect=(0, 0, 1, 0.90))
out = "/mnt/c/Users/bjark/siesta-build/_ghfix/repo/figures/convergence_before_after.png"
fig.savefig(out, dpi=135)
print("saved", out)
