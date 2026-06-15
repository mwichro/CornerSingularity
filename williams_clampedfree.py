# Williams singular exponent alpha for a right-angle (gamma = pi/2) elastic
# wedge with one edge CLAMPED (Dirichlet) and the opposite edge TRACTION-FREE
# (Neumann) -- the transverse pencil of the mixed corner of the paper.
#
# Method: Kolosov-Muskhelishvili eigenfunctions phi = C z^a, psi = E z^a
# (z = r e^{i theta}, a = exponent so that displacement ~ r^a). With C = c1+i c2,
# E = e1+i e2 this is 4 real DOF for the 4 homogeneous boundary conditions
#   theta = 0    (clamped):     u_r = u_theta = 0,
#   theta = pi/2 (traction-free): sigma_{theta theta} = sigma_{r theta} = 0.
# The secular equation is det M(a; kappa) = 0, kappa = 3 - 4 nu (plane strain).
# We report the smallest real root a in (0,1) as a function of nu and flag
# whether it sits above or below 1/2 (the convergence threshold for the
# coefficient integral 2 \int f'' d_phi^2 in 2D, where d_phi ~ r^{2(a-1)}).

import numpy as np

GAMMA = np.pi / 2.0  # right-angle opening; clamped at theta=0, free at theta=GAMMA


def bc_matrix(a, kappa):
    """Real 4x4 boundary-condition matrix in the unknowns (c1, c2, e1, e2)."""
    rows = []

    def disp_row(theta):
        # 2 mu (u_r + i u_theta) / r^a
        #   = kappa C e^{i(a-1)t} - a conj(C) e^{i(1-a)t} - conj(E) e^{-i(a+1)t}
        # Coefficients of C, conj(C), conj(E):
        cC = kappa * np.exp(1j * (a - 1) * theta)
        cCbar = -a * np.exp(1j * (1 - a) * theta)
        cEbar = -np.exp(-1j * (a + 1) * theta)
        # In terms of (c1,c2,e1,e2): C=c1+i c2, conj(C)=c1-i c2,
        # conj(E)=e1-i e2.
        coef = np.array([
            cC + cCbar,            # c1
            1j * (cC - cCbar),     # c2
            cEbar,                 # e1
            -1j * cEbar,           # e2
        ], dtype=complex)
        return coef

    def trac_row(theta):
        # (sigma_{tt} + i sigma_{rt}) / (a r^{a-1})
        #   = a C e^{i(a-1)t} + conj(C) e^{-i(a-1)t} + E e^{i(a+1)t}
        cC = a * np.exp(1j * (a - 1) * theta)
        cCbar = np.exp(-1j * (a - 1) * theta)
        cE = np.exp(1j * (a + 1) * theta)
        coef = np.array([
            cC + cCbar,            # c1
            1j * (cC - cCbar),     # c2
            cE,                    # e1
            1j * cE,               # e2
        ], dtype=complex)
        return coef

    # clamped edge theta = 0: real & imag parts of displacement vanish
    d0 = disp_row(0.0)
    rows.append(d0.real)
    rows.append(d0.imag)
    # free edge theta = GAMMA: real & imag parts of traction vanish
    tG = trac_row(GAMMA)
    rows.append(tG.real)
    rows.append(tG.imag)

    return np.array(rows)


def secular(a, kappa):
    return np.linalg.det(bc_matrix(a, kappa))


def smallest_root(kappa, lo=1e-4, hi=0.999, n=4000):
    xs = np.linspace(lo, hi, n)
    vals = np.array([secular(x, kappa) for x in xs])
    roots = []
    for i in range(len(xs) - 1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i + 1]) and vals[i] * vals[i + 1] < 0:
            a, b = xs[i], xs[i + 1]
            fa, fb = vals[i], vals[i + 1]
            for _ in range(80):
                m = 0.5 * (a + b)
                fm = secular(m, kappa)
                if fa * fm <= 0:
                    b, fb = m, fm
                else:
                    a, fa = m, fm
            roots.append(0.5 * (a + b))
    return roots


if __name__ == "__main__":
    print("Right-angle clamped--free wedge (gamma = pi/2), plane strain.")
    print(f"{'nu':>6} {'kappa':>8} {'alpha_min':>11}   {'side of 1/2':>12}")
    for nu in [0.0, 0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.49, 0.499]:
        kappa = 3.0 - 4.0 * nu
        roots = smallest_root(kappa)
        if roots:
            a = min(roots)
            side = "alpha > 1/2" if a > 0.5 else "alpha < 1/2"
            extra = "" if len(roots) == 1 else f"  (all in (0,1): {[round(x,4) for x in roots]})"
            print(f"{nu:6.3f} {kappa:8.3f} {a:11.5f}   {side:>12}{extra}")
        else:
            print(f"{nu:6.3f} {kappa:8.3f}   no real root in (0,1)")

    print()
    print("Cross-check: nu=0.3 known 90-deg clamped--free value ~ 0.7112 (lit.).")
    print()
    print("Homogeneity verdict (2D). The single-mode coefficient")
    print("  b_dir = 2 int_{B_rho} f''(J0) d_phi^2 dV,  d_phi = det grad(phi) ~ r^{2(alpha-1)},")
    print("has integrand d_phi^2 ~ r^{4(alpha-1)}, so int ... r dr ~ int r^{4 alpha - 3} dr,")
    print("which CONVERGES at the tip iff 4 alpha - 3 > -1, i.e. alpha > 1/2.")
    print("Since alpha(nu) in (0.60, 0.87) for nu in (0,1/2), we have alpha > 1/2 throughout:")
    print("the coefficient is tip-CONVERGENT but NOT tip-dominated -> b carries an O(1)")
    print("contribution from the matching annulus (the SIF / far field). Localization of b")
    print("to the bare wedge therefore FAILS; b must be computed on the graded-substrate")
    print("creasing model, which carries the extra dimensionless parameter K-hat.")
    print("(By contrast the threshold comparison t_c < t_coll, being scale-free, does localize.)")
