# Surface (Biot / creasing) instability threshold of a compressed compressible
# neo-Hookean half-space -- the buckling building block to which the corner
# instability reduces (Lemma biot). This furnishes lambda_star(nu) feeding (BC),
# and validates the incremental-buckling machinery against the classical
# incompressible value lambda_star ~ 0.544.
#
# 2D plane-strain neo-Hookean of the paper:
#   Psi(F) = (mu/2)(tr F^T F - 2) + f(J),  f(J) = -mu log J + Lam log^2 J,
#   P = mu F + f'(J) J F^{-T},
#   L_iJkL = mu d_ik d_JL + (f''J^2+f'J) Finv_iJ Finv_kL - f'J Finv_iL Finv_kJ,
# with f''J^2+f'J = 2 Lam,  f'J = -mu + 2 Lam log J.
#
# Base state: homogeneous F0 = diag(l1, l2), surface normal e2 traction-free
#   => P0_22 = 0 => mu l2^2 + f'(J)J = 0 => 2 Lam log(l1 l2) = mu(1 - l2^2).
# Incremental surface mode u = U exp(i k x1 + s k x2), decaying x2->-inf (Re s>0).
#   Q_ik(s) = -L_i1k1 + i s (L_i1k2 + L_i2k1) + s^2 L_i2k2,  det Q = 0 -> quartic.
# Combine the two decaying roots, impose incremental traction t_i = 0 at x2=0:
#   t_i = sum_a c_a [ i L_i2k1 + s_a L_i2k2 ] U^a_k = 0  -> 2x2 secular det.

import numpy as np
from scipy.optimize import brentq


def base_l2(l1, mu, Lam):
    """Solve 2 Lam log(l1 l2) = mu (1 - l2^2) for l2 in (0, ~2]."""
    def g(l2):
        return 2.0 * Lam * np.log(l1 * l2) - mu * (1.0 - l2 ** 2)
    lo, hi = 1e-6, 50.0
    return brentq(g, lo, hi, xtol=1e-14)


def tangent(l1, l2, mu, Lam):
    """4-index incremental modulus L[i,J,k,L] at F0 = diag(l1,l2)."""
    J = l1 * l2
    A2 = 2.0 * Lam               # f''J^2 + f'J
    A1 = -mu + 2.0 * Lam * np.log(J)   # f'J
    Fi = np.array([1.0 / l1, 1.0 / l2])
    Finv = np.diag(Fi)
    L = np.zeros((2, 2, 2, 2))
    for i in range(2):
        for Jx in range(2):
            for k in range(2):
                for Lx in range(2):
                    L[i, Jx, k, Lx] = (
                        mu * (i == k) * (Jx == Lx)
                        + A2 * Finv[i, Jx] * Finv[k, Lx]
                        - A1 * Finv[i, Lx] * Finv[k, Jx]
                    )
    return L


def secular(l1, mu, Lam):
    """Real surface-instability secular determinant at stretch l1."""
    l2 = base_l2(l1, mu, Lam)
    L = tangent(l1, l2, mu, Lam)

    # Q(s) = Q0 + Q1 s + Q2 s^2 (2x2, complex);  det Q(s) is a quartic in s.
    Q0 = -L[:, 0, :, 0]
    Q1 = 1j * (L[:, 0, :, 1] + L[:, 1, :, 0])
    Q2 = L[:, 1, :, 1]

    # det Q(s) = det(Q2) s^4 + ... ; assemble coefficients via polynomial of 2x2.
    # det(Q0+Q1 s+Q2 s^2) for 2x2 = a(s)d(s) - b(s)c(s), each entry quadratic in s.
    def entrypoly(M0, M1, M2, i, j):
        return np.array([M2[i, j], M1[i, j], M0[i, j]])  # coeff of s^2,s^1,s^0
    a = entrypoly(Q0, Q1, Q2, 0, 0)
    b = entrypoly(Q0, Q1, Q2, 0, 1)
    c = entrypoly(Q0, Q1, Q2, 1, 0)
    d = entrypoly(Q0, Q1, Q2, 1, 1)
    detpoly = np.convolve(a, d) - np.convolve(b, c)   # quartic, coeff s^4..s^0
    roots = np.roots(detpoly)
    dec = roots[roots.real > 1e-9]                  # decaying modes (Re s>0)
    if len(dec) != 2:
        # pick the two with largest positive real part as a fallback
        dec = roots[np.argsort(-roots.real)][:2]

    # null vectors and surface-traction operator for each decaying root
    cols = []
    for s in dec:
        Q = Q0 + Q1 * s + Q2 * s * s
        # null vector of 2x2 Q
        u = np.array([Q[0, 1], -Q[0, 0]])
        if np.linalg.norm(u) < 1e-12:
            u = np.array([Q[1, 1], -Q[1, 0]])
        u = u / np.linalg.norm(u)
        # traction row operator T_ik = i L_i2k1 + s L_i2k2
        T = 1j * L[:, 1, :, 0] + s * L[:, 1, :, 1]
        cols.append(T @ u)
    Smat = np.column_stack(cols)
    det = np.linalg.det(Smat)
    return det


def secular_real(l1, mu, Lam):
    """The secular determinant is purely imaginary; return its real coefficient."""
    return secular(l1, mu, Lam).imag


def lambda_star(mu, Lam):
    """Largest l1<1 at which the secular determinant changes sign."""
    l1s = np.linspace(0.3, 0.98, 400)
    vals = np.array([secular_real(l, mu, Lam) for l in l1s])
    for i in range(len(l1s) - 1, 0, -1):
        if np.isfinite(vals[i]) and np.isfinite(vals[i - 1]) and vals[i] * vals[i - 1] < 0:
            return brentq(lambda l: secular_real(l, mu, Lam), l1s[i - 1], l1s[i], xtol=1e-10)
    return np.nan


if __name__ == "__main__":
    mu = 1.0
    print("Compressible neo-Hookean surface-instability threshold lambda_star")
    print("(2 Lam = f''J^2+f'J; large Lam/mu -> incompressible Biot ~ 0.544)\n")
    print(f"{'Lam/mu':>8} {'nu(eff)':>9} {'lambda_star':>12}")
    for Lam in [0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 200.0, 1000.0]:
        # plane-strain effective Poisson nu = Lam/(2(Lam+mu))
        nu = Lam / (2.0 * (Lam + mu))
        ls = lambda_star(mu, Lam)
        print(f"{Lam/mu:8.1f} {nu:9.4f} {ls:12.6f}")
    print("\nReference: incompressible neo-Hookean Biot surface mode l_star = 0.5437.")
