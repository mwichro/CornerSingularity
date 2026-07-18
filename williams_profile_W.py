# Numerical verification of the profile-nondegeneracy hypothesis (W) for the
# right-angle clamped--free Williams mode (the tensor ansatz behind the vertex
# screening lemma):
#
#   G(theta) = alpha * Phi (x) e_r + Phi' (x) e_theta,
#   mu(theta) = det G = alpha * (Phi x Phi'),   trG = alpha Phi.e_r + Phi'.e_theta.
#
# (W) requires: mu >= 0 on [0, pi/2] with finitely many finite-order zeros and
# trG != 0 at each zero. mu's sign is intrinsic (quadratic in the mode), and
# J ~ K^2 r^{2(alpha-1)} mu at depth for EITHER load sign, so mu < 0 on a
# sector means the linear ansatz is orientation-reversing at depth there.
#
# Mode construction: null vector of the Kolosov-Muskhelishvili boundary-condition
# matrix of williams_clampedfree.py at the smallest root alpha(nu); the angular
# displacement profile is
#   u_x + i u_y = r^a e^{i th} W(th) / (2 mu_e),
#   W = kappa C e^{i(a-1)th} - a conj(C) e^{i(1-a)th} - conj(E) e^{-i(a+1)th},
# with analytic theta-derivative (checked against finite differences to 1e-9).
#
# VERDICT (recorded 2026-07-18, from the run below):
#   nu_W = 0.4042.  For nu in (nu_W, 1/2): mu >= 0 throughout, single order-2
#   zero at the clamped face (forced by Phi(0)=0), trG(0) != 0 -> (W) HOLDS.
#   For nu < nu_W: one interior zero theta_z(nu); mu < 0 on the clamped-side
#   sector (0, theta_z), with theta_z sweeping monotonically from the clamped
#   face (nu -> nu_W) to the free face as nu decreases (negative on the whole
#   opening for nu <~ 0.08). trG != 0 at every zero (layers fine), but the sign
#   clause of (W) FAILS on an O(1) sector -> the screened-vertex conclusion is
#   established only for nu > nu_W; below, the inner state on the negative
#   sector must leave the linear ansatz at depth (self-consistency caveat).

import numpy as np

from williams_clampedfree import bc_matrix, smallest_root


def profile_fields(nu, th):
    """alpha, normalised det G, normalised trG, singular values of G on grid th."""
    kappa = 3.0 - 4.0 * nu
    roots = smallest_root(kappa)
    if not roots:
        raise ValueError(f"no Williams root in (0,1) at nu={nu}")
    a = min(roots)
    M = bc_matrix(a, kappa)
    v = np.linalg.svd(M)[2][-1]  # null vector (residual ~ 1e-16)
    C, E = v[0] + 1j * v[1], v[2] + 1j * v[3]
    W = (kappa * C * np.exp(1j * (a - 1) * th)
         - a * np.conj(C) * np.exp(1j * (1 - a) * th)
         - np.conj(E) * np.exp(-1j * (a + 1) * th))
    Wp = ((1j * (a - 1)) * kappa * C * np.exp(1j * (a - 1) * th)
          - (1j * (1 - a)) * a * np.conj(C) * np.exp(1j * (1 - a) * th)
          + (1j * (a + 1)) * np.conj(E) * np.exp(-1j * (a + 1) * th))
    Z, Zp = np.exp(1j * th) * W, np.exp(1j * th) * (1j * W + Wp)
    Phi = np.stack([Z.real, Z.imag], -1)
    Phip = np.stack([Zp.real, Zp.imag], -1)
    er = np.stack([np.cos(th), np.sin(th)], -1)
    eth = np.stack([-np.sin(th), np.cos(th)], -1)
    detG = a * (Phi[:, 0] * Phip[:, 1] - Phi[:, 1] * Phip[:, 0])
    trG = a * np.einsum('ij,ij->i', Phi, er) + np.einsum('ij,ij->i', Phip, eth)
    G = (a * np.einsum('ti,tj->tij', Phi, er)
         + np.einsum('ti,tj->tij', Phip, eth))
    sv = np.linalg.svd(G, compute_uv=False)
    nrm = sv[:, 0].max()
    return a, detG / nrm**2, trG / nrm, sv / nrm


def theta2_coeff_at_clamped_face(nu):
    """Coefficient c2 in det G ~ c2 * theta^2 at theta = 0 (order-2 forced zero)."""
    th = np.array([1e-4, 2e-4])
    _, d, _, _ = profile_fields(nu, th)
    return float(np.mean(d / th**2))


def nu_W(lo=0.35, hi=0.45, iters=40):
    """Bisection for the sign change of c2(0): (W) holds for nu > nu_W."""
    flo = theta2_coeff_at_clamped_face(lo)
    for _ in range(iters):
        m = 0.5 * (lo + hi)
        if flo * theta2_coeff_at_clamped_face(m) <= 0:
            hi = m
        else:
            lo, flo = m, theta2_coeff_at_clamped_face(m)
    return 0.5 * (lo + hi)


if __name__ == "__main__":
    th = np.linspace(0.0, np.pi / 2, 40001)
    print("Hypothesis (W) check, right-angle clamped--free Williams profile.")
    print(f"nu_W (sign change of the theta^2 coefficient at the clamped face): "
          f"{nu_W():.4f}")
    print()
    print(f"{'nu':>6} {'alpha':>8} {'min detG~':>10} {'theta_z/(pi/2)':>14} "
          f"{'|trG| at zeros~':>15}  verdict")
    for nu in [0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40,
               0.4042, 0.417, 0.45, 0.49]:
        a, d, t, sv = profile_fields(nu, th)
        sc = [i for i in range(200, len(th) - 1) if d[i] * d[i + 1] < 0]
        zs = ", ".join(f"{th[i]/(np.pi/2):.3f}" for i in sc) or "-"
        tz = [abs(t[0])] + [abs(t[i]) for i in sc]
        verdict = "(W) holds" if (d.min() > -1e-10 and not sc) else \
                  "(W) FAILS on clamped-side sector"
        print(f"{nu:6.4f} {a:8.5f} {d.min():10.2e} {zs:>14} "
              f"{', '.join(f'{x:.3f}' for x in tz):>15}  {verdict}")
    print()
    print("Normalisation: detG by max|G|^2, trG by max|G| over theta (mode is")
    print("defined up to a real scalar; detG and the vanishing of trG are")
    print("scale- and sign-invariant). First |trG| column entry is the clamped")
    print("face theta=0; subsequent entries are the interior zeros theta_z.")
