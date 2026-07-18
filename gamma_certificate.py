# The inversion-free Gamma certificate of cor:cert / prop:certdeg, evaluated on
# the existing graded-cell discretisation (GradedCell of graded_buckling.py).
# This is the computation the paper recommends (rem:certnum) but had never run.
#
# Certificate (cor:cert): for any test field v in V^perp with <v,Qv> > 0,
#     b <= b_dir (1 - Gamma(v)),   Gamma(v) := E3[phi,phi,v]^2 / (2 b_dir <v,Qv>),
# so Gamma(v) > 1 certifies b < 0 with NO bordered solve.  In GradedCell's
# variables: the assembled load vector R obeys <R,v> = (1/2) E3[phi,phi,v] and
# the paper's frozen quartic b_dir equals the script's 4*E4, hence
#     Gamma(v) = 2 (R.v)^2 / (b_dir * v^T Q v),
# invariant under rescaling of phi and of v.  The supremum over v is
# Gamma_max = -b_psi/b_dir = -r  (r the ratio of the b tables), attained at the
# Lyapunov-Schmidt corrector; the point of the certificate is that CHEAP explicit
# v already exceed 1:
#   * v = R itself (projected to V^perp): zero linear solves;
#   * v = k-step conjugate-gradient iterates on Q v = R in V^perp: each iterate
#     is an explicit polynomial in Q applied to R, so the certificate remains
#     inversion-free (no near-singular bordered matrix; cond(B) never enters).
# We also evaluate v = phi_2 (the first excited mode), which probes the two
# hypotheses of prop:certdeg directly:
#   (a) neutrality:  neut = <phi2,Q phi2>/<phi2,M phi2>  (-> 0 in the flat limit?)
#   (b) coupling:    coup = E3[phi,phi,phi2]^2/(b_dir <phi2,M phi2>) >= c0 > 0.
#
# Caveats inherited from the surrogate cell (recorded in numerical.tex sub:numb):
# this is the flat graded half-space cell, not the prescribed wedge solve, and
# the discrete gap decays under refinement; Gamma certifies the sign of the
# DISCRETE b of this cell, at each tested resolution, nothing more.

import numpy as np

from graded_buckling import GradedCell, nu_of


def certificate_data(cell, l1bar, Khat, grading='linear', kmax=1000,
                     checkpoints=(1, 3, 10, 40, 160, 640, 1000)):
    """Assemble (phi, R, b_dir, Q) at the critical state and evaluate Gamma on
    the no-solve field v=R, on CG iterates, and on the excited mode phi2."""
    mat = cell.material(l1bar, Khat, grading)
    Q = cell.stiffness(mat)
    keep = cell.keep
    Qk = Q[np.ix_(keep, keep)]

    from scipy.linalg import eigh
    ev, V = eigh(Qk, subset_by_index=[0, 1])
    phi_k, phi2_k = V[:, 0], V[:, 1]
    phi = np.zeros(2 * cell.m)
    phi[keep] = phi_k
    gphi = np.array([cell.gradapply(p, phi) for p in range(4)])
    scale = np.sqrt(np.max(np.sum(gphi ** 2, axis=0)))
    phi, gphi = phi / scale, gphi / scale
    phi_k = phi[keep]

    # source R and frozen quartic b_dir = 4*E4, exactly as in landau_b
    from graded_buckling import MMAT
    cv_nd = cell._tile(mat['cv']).T
    fpp = cell._tile(mat['fpp'])
    fppp = cell._tile(mat['fppp'])
    fpppp = cell._tile(mat['fpppp'])
    p_ = np.sum(cv_nd * gphi, axis=0)
    dphi = gphi[0] * gphi[3] - gphi[1] * gphi[2]
    scal = fpp * dphi + 0.5 * fppp * p_ ** 2
    Sfield = fpp * p_[None, :] * (MMAT @ gphi) + scal[None, :] * cv_nd
    R = np.zeros(2 * cell.m)
    for q in range(4):
        iq, Dq = cell.comp_of[q], cell.Dop[cell.der_of[q]]
        R[iq * cell.m:(iq + 1) * cell.m] += Dq.T @ (cell.omega * Sfield[q])
    Rk = R[keep]
    e4 = (0.5 * fpp * dphi ** 2 + 0.5 * fppp * p_ ** 2 * dphi
          + (1.0 / 24.0) * fpppp * p_ ** 4)
    bdir = 4.0 * np.sum(e4 * cell.omega)

    # mass-orthogonal projector onto V^perp (constraint of the LS decomposition)
    Mphi_k = (cell.omega2 * phi)[keep]
    mnorm = Mphi_k @ phi_k

    def proj(v):
        return v - phi_k * (Mphi_k @ v) / mnorm

    def gamma(v):
        qv = v @ (Qk @ v)
        return 2.0 * (Rk @ v) ** 2 / (bdir * qv) if qv > 0 else np.nan

    # deflated conjugate gradients on Qk v = Rk in V^perp.  By Krylov
    # optimality Gamma_k is nondecreasing in k, and EVERY iterate certifies the
    # rigorous bound b <= b_dir (1 - Gamma_k); k1 is the first certifying k.
    gammas = {}
    k1 = None
    b = proj(Rk)
    x = np.zeros_like(b)
    r = b.copy()
    p = r.copy()
    rs = r @ r
    for k in range(1, kmax + 1):
        Ap = proj(Qk @ p)
        alpha = rs / (p @ Ap)
        x = x + alpha * p
        r = r - alpha * Ap
        rs_new = r @ r
        p = r + (rs_new / rs) * p
        rs = rs_new
        if k1 is None and gamma(proj(x)) > 1.0:
            k1 = k
        if k in checkpoints:
            gammas[k] = gamma(proj(x))
    gamma_R = gamma(proj(Rk))

    # excited mode: neutrality and coupling of prop:certdeg
    phi2 = proj(phi2_k)
    M2 = (cell.omega2[keep] if cell.omega2.ndim == 1 else None)
    mass2 = (cell.omega2[keep] * phi2) @ phi2
    neut = phi2 @ (Qk @ phi2) / mass2
    coup = 2.0 * (Rk @ phi2) ** 2 / (bdir * mass2)   # = E3^2/(2 b_dir ||phi2||_M^2)
    gamma_2 = gamma(phi2)

    return dict(gap=ev[1] - ev[0], kernel_eig=ev[0], bdir=bdir,
                gamma_R=gamma_R, gammas=gammas, gamma_2=gamma_2,
                neut=neut, coup=coup, k1=k1)


if __name__ == "__main__":
    mu = 1.0
    print("=" * 72)
    print("Gamma certificate (cor:cert) on the graded cell: b<0 iff certified")
    print("Gamma>1; Gamma(v) inversion-free for v=R and CG iterates of R.")
    print("=" * 72)

    print("\n(1) Rectangle sweep (grid of section D: NX=10, N=24, linear grading)")
    print(f"{'nu':>7} {'Khat':>6} {'l1c':>8} {'gap':>10} {'G(R)':>8} "
          f"{'G_cg40':>8} {'G_cg160':>8} {'G_cg1000':>9} {'k1':>5}  verdict")
    for Lam in [1.0, 2.0, 5.0, 10.0]:
        cell = GradedCell(mu, Lam, NX=10, N=24)
        for Khat in [0.10, 0.20, 0.30]:
            lc = cell.critical_l1(Khat)
            if not np.isfinite(lc):
                print(f"{nu_of(Lam,mu):7.3f} {Khat:6.2f}   (no critical l1)")
                continue
            d = certificate_data(cell, lc, Khat)
            g = d['gammas']
            ok = "CERTIFIED b<0" if g[1000] > 1 else "not certified"
            print(f"{nu_of(Lam,mu):7.3f} {Khat:6.2f} {lc:8.4f} {d['gap']:10.2e} "
                  f"{d['gamma_R']:8.3f} {g[40]:8.3f} {g[160]:8.3f} {g[1000]:9.3f} "
                  f"{str(d['k1']):>5}  {ok}")

    print("\n(2) Resolution robustness of the certificate (nu=0.417, Khat=0.15)")
    print(f"{'NX':>4} {'N':>4} {'gap':>10} {'G(R)':>8} {'G_cg160':>8} "
          f"{'G_cg640':>8} {'G_cg2000':>9} {'k1':>5}")
    for NX, N in [(10, 24), (12, 28), (14, 32), (16, 36)]:
        cell = GradedCell(mu, 5.0, NX=NX, N=N)
        lc = cell.critical_l1(0.15)
        d = certificate_data(cell, lc, 0.15, kmax=2000,
                             checkpoints=(1, 40, 160, 640, 2000))
        g = d['gammas']
        print(f"{NX:4d} {N:4d} {d['gap']:10.2e} {d['gamma_R']:8.3f} "
              f"{g[160]:8.3f} {g[640]:8.3f} {g[2000]:9.3f} {str(d['k1']):>5}")

    print("\n(3) prop:certdeg hypotheses as Khat -> 0 (nu=0.417, NX=12, N=28):")
    print("    (a) neut -> 0? (excited-mode neutrality)   (b) coup >= c0 > 0?")
    print(f"{'Khat':>6} {'l1c':>8} {'gap':>10} {'neut':>10} {'coup':>10} "
          f"{'G(phi2)':>8}")
    cell = GradedCell(mu, 5.0, NX=12, N=28)
    for Khat in [0.02, 0.05, 0.10, 0.20, 0.30]:
        lc = cell.critical_l1(Khat)
        if not np.isfinite(lc):
            print(f"{Khat:6.2f}   (no critical l1 in window)")
            continue
        d = certificate_data(cell, lc, Khat)
        print(f"{Khat:6.2f} {lc:8.4f} {d['gap']:10.2e} {d['neut']:10.3e} "
              f"{d['coup']:10.3e} {d['gamma_2']:8.3f}")
