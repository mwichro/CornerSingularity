# Plot the leading Williams singular exponent alpha(nu) for the right-angle
# clamped--free (Dirichlet--Neumann) wedge, and mark the alpha = 1/2 homogeneity
# threshold that decides whether the reduced cubic localises to the bare tip.
#
# Reuses the secular-determinant solver from williams_clampedfree.py.

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from williams_clampedfree import smallest_root

nus = np.linspace(0.001, 0.499, 250)
alphas = []
for nu in nus:
    roots = smallest_root(3.0 - 4.0 * nu)
    alphas.append(min(roots) if roots else np.nan)
alphas = np.array(alphas)

fig, ax = plt.subplots(figsize=(6.2, 4.0))

# alpha > 1/2 region: cubic integral tip-CONVERGENT but not tip-dominated (b = b(nu, Khat)).
ax.axhspan(0.5, 1.0, color="#d9e8f5", alpha=0.7, zorder=0)
# alpha < 1/2 region: tip-DOMINATED, b would be a universal wedge functional.
ax.axhspan(0.0, 0.5, color="#f5e3d9", alpha=0.7, zorder=0)

ax.plot(nus, alphas, color="#103d6b", lw=2.2, zorder=3,
        label=r"$\alpha(\nu)$, clamped--free $90^\circ$")
ax.axhline(0.5, color="k", ls="--", lw=1.2, zorder=2)

# annotate the two regimes
ax.text(0.05, 0.93, r"$\alpha>\frac{1}{2}$: tip-convergent, not tip-dominated"
        "\n" r"$\Rightarrow\ b=b(\nu,\hat K)$ (phase diagram)",
        fontsize=8.5, va="top", color="#103d6b")
ax.text(0.05, 0.16, r"$\alpha<\frac{1}{2}$: tip-dominated"
        "\n" r"$\Rightarrow\ b$ universal wedge functional (not realised here)",
        fontsize=8.5, va="top", color="#7a3b16")

# reference markers
for nu0 in [0.3, 0.45, 0.499]:
    a0 = min(smallest_root(3.0 - 4.0 * nu0))
    ax.plot([nu0], [a0], "o", color="#103d6b", ms=4, zorder=4)
    off = (-78, -4) if nu0 > 0.48 else (4, 6)
    ax.annotate(rf"$({nu0:g},\,{a0:.3f})$", (nu0, a0),
                textcoords="offset points", xytext=off, fontsize=7.5)

ax.set_xlabel(r"Poisson ratio $\nu$ (plane strain, $\kappa=3-4\nu$)")
ax.set_ylabel(r"leading singular exponent $\alpha$")
ax.set_xlim(0.0, 0.5)
ax.set_ylim(0.0, 1.0)
ax.set_title(r"Williams exponent at the clamped--free corner and the $\alpha=\frac{1}{2}$ threshold",
             fontsize=10)
ax.legend(loc="upper right", fontsize=8.5, framealpha=0.9)
ax.grid(True, color="white", lw=0.6, zorder=1)

fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"figures/williams_exponent.{ext}", dpi=200)
print("wrote figures/williams_exponent.pdf and .png")
print(f"alpha range over nu in (0,1/2): [{np.nanmin(alphas):.4f}, {np.nanmax(alphas):.4f}]")
