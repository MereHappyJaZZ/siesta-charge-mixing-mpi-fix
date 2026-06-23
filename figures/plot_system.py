import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

FDF = "/home/bjarke/rc210/STRUCT_DEVICE.fdf"
lines = open(FDF).read().splitlines()

cell = []
for i, l in enumerate(lines):
    if "LatticeVectors" in l and "%block" in l:
        for j in range(1, 4):
            cell.append([float(x) for x in lines[i+j].split()[:3]])
        break
cell = np.array(cell); AX = cell[0, 0]   # periodic repeat along x

X, Y, Z, E = [], [], [], []
inblk = False
for l in lines:
    if "AtomicCoordinatesAndAtomicSpecies" in l and "%block" in l:
        inblk = True; continue
    if inblk and "%endblock" in l:
        break
    if inblk:
        p = l.split()
        if len(p) < 4: continue
        el = "?"
        if "#" in l and ":" in l.split("#")[1]:
            el = l.split("#")[1].split(":")[1].strip().split()[0]
        X.append(float(p[0])); Y.append(float(p[1])); Z.append(float(p[2])); E.append(el)
X, Y, Z = map(np.array, (X, Y, Z)); E = np.array(E)

COL = {"Cr": "#7f9bb3", "C": "#2b2b2b", "O": "#d23b22", "F": "#46b04a", "H": "#dddddd"}
SZ  = {"Cr": 120, "C": 46, "O": 70}

def cc_bonds(ax, h, v, x3, y3, z3, els, cut=1.85):
    idx = [k for k in range(len(els)) if els[k] in ("C", "O")]
    for a in range(len(idx)):
        for b in range(a+1, len(idx)):
            i, j = idx[a], idx[b]
            d = ((x3[i]-x3[j])**2 + (y3[i]-y3[j])**2 + (z3[i]-z3[j])**2) ** 0.5
            if d < cut:
                ax.plot([h[i], h[j]], [v[i], v[j]], "-", color="#8a8a8a", lw=1.2, zorder=1)

def spheres(ax, h, v, els, scale=1.0):
    for e in sorted(set(els), key=lambda x: -SZ.get(x, 40)):
        m = els == e
        ax.scatter(np.array(h)[m], np.array(v)[m], s=SZ.get(e, 40)*scale,
                   c=COL.get(e, "#888"), edgecolors="k", linewidths=0.4,
                   label=f"{e} ({(E==e).sum()})", zorder=3)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11.5, 6.4))

# --- Panel A: full device, side view (z transport vs y) ---
cc_bonds(ax1, Z, Y, X, Y, Z, E)
spheres(ax1, Z, Y, E)
ax1.set_xlabel("z   (Å)  —  transport direction  →"); ax1.set_ylabel("y  (Å)")
ax1.set_title("Graphene–Cr edge-contact device cell  "
              f"({len(E)} atoms;  cell {cell[0,0]:.2f} × {cell[1,1]:.2f} × {cell[2,2]:.1f} Å)")
ax1.annotate("Cr contact", (7, 13.5), fontsize=9, ha="center", color="#3a5a78")
ax1.annotate("O", (16, 4.2), fontsize=9, ha="center", color="#d23b22")
ax1.annotate("graphene channel (edge-on, 160 C)", (105, 4.0), fontsize=9, ha="center", color="#2b2b2b")
ax1.legend(loc="lower right", fontsize=8, ncol=3, framealpha=0.95)
ax1.set_aspect("equal", adjustable="datalim"); ax1.set_ylim(-3, 17)

# --- Panel B: top view of the contact (z vs x), x-periodic images replicated to show honeycomb ---
reg = (Z >= 8) & (Z <= 58)
zr, xr, yr, er = [], [], [], []
for s in (-AX, 0.0, AX):
    zr += list(Z[reg]); xr += list(X[reg] + s); yr += list(Y[reg]); er += list(E[reg])
zr, xr, yr = map(np.array, (zr, xr, yr)); er = np.array(er)
cc_bonds(ax2, zr, xr, xr, yr, zr, er)
spheres(ax2, zr, xr, er, scale=1.3)
ax2.set_xlabel("z   (Å)"); ax2.set_ylabel("x  (Å, periodic — 1 image each side)")
ax2.set_title("Top view of the Cr / graphene contact (honeycomb channel meeting the metal)")
ax2.set_aspect("equal", adjustable="datalim")
ax2.axhline(0, color="0.85", lw=0.6); ax2.axhline(AX, color="0.85", lw=0.6)

fig.tight_layout()
out = "/mnt/c/Users/bjark/siesta-build/_ghfix/repo/figures/system_structure.png"
fig.savefig(out, dpi=135)
print("saved", out)
