# First computation of the reduced QUINTIC coefficient d of ass:subcrit on the
# graded cell -- the least-supported clause of the paper (d>0), for which no
# numerical value existed anywhere (rem:fbstatus).
#
# Sixth-order Lyapunov-Schmidt at (t_c, eps=0), symmetric case (c2=0).  With
# psi(w) = w^2 psi2 + w^3 psi3 (psi4 does not enter at order w^6, by
# stationarity of the reduced energy in psi), the reduced energy coefficient is
#
#   V6 = (1/720) E6[phi^6] + (1/24) E5[phi^4 psi2] + (1/4) E4[phi^2 psi2^2]
#      + (1/6) E3[psi2^3] + <S3, psi3> + (1/2) <psi3, L psi3>,       d = 6 V6,
#
# where psi2 = -(QLQ)^{-1} R (the chi of landau_b),
#       S3   = Q( E3[phi, psi2, .] + (1/6) E4[phi,phi,phi, .] ),
#       psi3 = -(QLQ)^{-1} S3   (so the last two terms = -1/2 <S3,(QLQ)^{-1}S3> <= 0:
#       the third-order feedback is ALWAYS destabilising for d, the exact mirror
#       of lem:fbsign at fourth order; the E5/E4/E3 middle terms are unsigned).
#
# In 2D the determinant is quadratic in the gradient, so all E_k (k>=3) are
# purely volumetric and closed-form in
#   q_h  = cof F0 : grad h          (the script's cv . g_h),
#   m_hk = polarised det            (= 1/2 g_h^T MMAT g_k, m_hh = det grad h),
# via Faa di Bruno.  Explicitly (all integrals against the cell weight):
#   E3[a,b,c]  = f''' qa qb qc + 2 f'' (qa m_bc + qb m_ac + qc m_ab)
#   E4[a,b,c,e]= f'''' q^4 + 2 f''' (6 pairings q q m) + 4 f'' (3 pairings m m)
#   E5, E6     = analogous, with multiplicities 2^k (# pairings); for equal
#   arguments: E6[phi^6] = f^(6) q^6 + 30 f^(5) q^4 d + 180 f'''' q^2 d^2
#              + 120 f''' d^3   (d = det grad phi), whose D1->0 limit recovers
#              d_dir = int f''' d_phi^3 of lem:volcoeff(ii)  [validation 1].
# Validation 2: V4 = (1/24)E4[phi^4] + (1/2)E3[phi,phi,psi2] + (1/2)<psi2,L psi2>
#              must reproduce b/4 with b from landau_b's independent assembly.
#
# f_log derivatives: f^(n) = (a_n + b_n log J)/J^n with a_{n+1} = b_n - n a_n,
# b_{n+1} = -n b_n, a1 = -mu, b1 = 2 Lam:
#   f^(5) = (-24 mu - 100 Lam + 48 Lam log J)/J^5
#   f^(6) = (120 mu + 548 Lam - 240 Lam log J)/J^6   (FD-checked).
#
# Caveats: same surrogate cell as the b tables (linear monotone grading, not the
# wedge; decaying gap under refinement); d, like b, scales as the 6th power of
# the mode amplitude, so we also report the invariant fold amplitude
# w_L^2 = -b/(2d) (paper eq:fold), meaningful at max|grad phi| = 1.

import numpy as np

from graded_buckling import GradedCell, MMAT, nu_of


def flog_derivs_56(J, mu, Lam):
    L = np.log(J)
    f5 = (-24.0 * mu - 100.0 * Lam + 48.0 * Lam * L) / J ** 5
    f6 = (120.0 * mu + 548.0 * Lam - 240.0 * Lam * L) / J ** 6
    return f5, f6


def quintic_data(cell, l1bar, Khat, grading='linear'):
    mat = cell.material(l1bar, Khat, grading)
    Q = cell.stiffness(mat)
    keep = cell.keep
    Qk = Q[np.ix_(keep, keep)]
    w = cell.omega

    from scipy.linalg import eigh
    ev, V = eigh(Qk, subset_by_index=[0, 1])
    phi = np.zeros(2 * cell.m)
    phi[keep] = V[:, 0]
    gphi = np.array([cell.gradapply(p, phi) for p in range(4)])
    scale = np.sqrt(np.max(np.sum(gphi ** 2, axis=0)))
    phi, gphi = phi / scale, gphi / scale

    cv = cell._tile(mat['cv']).T                       # (4, nodes)
    f2 = cell._tile(mat['fpp'])
    f3 = cell._tile(mat['fppp'])
    f4 = cell._tile(mat['fpppp'])
    J0 = cell._tile(mat['J0'])
    f5, f6 = flog_derivs_56(J0, cell.mu, cell.Lam)

    def qof(g):   return np.sum(cv * g, axis=0)
    def mof(ga, gb): return 0.5 * np.sum(ga * (MMAT @ gb), axis=0)

    qP = qof(gphi)
    dP = mof(gphi, gphi)                               # det grad phi

    # ---- load-vector assembly:  <ell, v> = int [A q_v + sum_h B_h m_{h,v}] w
    def load(A, Bs):
        """A: scalar field (pairs with q_v); Bs: list of (field, grad) pairing
        with m_{grad,v}."""
        Ltot = np.zeros(2 * cell.m)
        S = A[None, :] * cv
        for B, gh in Bs:
            S = S + B[None, :] * 0.5 * (MMAT @ gh)
        for qq in range(4):
            iq, Dq = cell.comp_of[qq], cell.Dop[cell.der_of[qq]]
            Ltot[iq * cell.m:(iq + 1) * cell.m] += Dq.T @ (w * S[qq])
        return Ltot

    # ---- order-2 source R  (as in landau_b) and psi2 via the bordered solve
    R = load(f2 * dP + 0.5 * f3 * qP ** 2, [(2.0 * f2 * qP, gphi)])
    Mphi = cell.omega2 * phi
    Mk = Mphi[keep]
    nk = len(keep)
    B = np.zeros((nk + 1, nk + 1))
    B[:nk, :nk] = Qk
    B[:nk, nk] = Mk
    B[nk, :nk] = Mk

    def bordered(rhs_full):
        rhs = np.zeros(nk + 1)
        rhs[:nk] = -rhs_full[keep]
        sol = np.linalg.solve(B, rhs)
        out = np.zeros(2 * cell.m)
        out[keep] = sol[:nk]
        return out

    psi2 = bordered(R)
    g2 = np.array([cell.gradapply(p, psi2) for p in range(4)])
    q2 = qof(g2)
    d2 = mof(g2, g2)                                   # det grad psi2
    mP2 = mof(gphi, g2)                                # m_{phi, psi2}

    # ---- b two ways (validation 2)
    E4_p4 = np.sum((f4 * qP ** 4 + 12 * f3 * qP ** 2 * dP + 12 * f2 * dP ** 2) * w)
    E3_pp2 = np.sum((f3 * qP ** 2 * q2
                     + 2 * f2 * (2 * qP * mP2 + q2 * dP)) * w)
    psi2Lpsi2 = psi2[keep] @ (Qk @ psi2[keep])
    V4 = E4_p4 / 24.0 + 0.5 * E3_pp2 + 0.5 * psi2Lpsi2
    b_ls = 4.0 * V4
    # landau_b-style independent assembly: b = (1/6)E4[phi^4] + 2 <R, psi2>
    b_direct = E4_p4 / 6.0 + 2.0 * (R @ psi2)

    # ---- S3 and psi3
    A_S3 = (f3 * qP * q2 + 2 * f2 * mP2
            + f4 * qP ** 3 / 6.0 + f3 * qP * dP)
    B_phi = 2 * f2 * q2 + f3 * qP ** 2 + 2 * f2 * dP
    B_psi = 2 * f2 * qP
    S3 = load(A_S3, [(B_phi, gphi), (B_psi, g2)])
    psi3 = bordered(S3)
    g3 = np.array([cell.gradapply(p, psi3) for p in range(4)])
    psi3Lpsi3 = psi3[keep] @ (Qk @ psi3[keep])
    T_psi3 = (S3 @ psi3) + 0.5 * psi3Lpsi3            # = -1/2 <S3,(QLQ)^{-1}S3> <= 0

    # ---- the four explicit V6 pieces
    E6_p6 = np.sum((f6 * qP ** 6 + 30 * f5 * qP ** 4 * dP
                    + 180 * f4 * qP ** 2 * dP ** 2 + 120 * f3 * dP ** 3) * w)
    E5_p4s = np.sum((f5 * qP ** 4 * q2
                     + 2 * f4 * (6 * qP ** 2 * q2 * dP + 4 * qP ** 3 * mP2)
                     + 4 * f3 * (3 * q2 * dP ** 2 + 12 * qP * dP * mP2)) * w)
    E4_p2s2 = np.sum((f4 * qP ** 2 * q2 ** 2
                      + 2 * f3 * (q2 ** 2 * dP + qP ** 2 * d2 + 4 * qP * q2 * mP2)
                      + 4 * f2 * (dP * d2 + 2 * mP2 ** 2)) * w)
    E3_s3 = np.sum((f3 * q2 ** 3 + 6 * f2 * q2 * d2) * w)

    d_frozen = 6.0 * E6_p6 / 720.0                    # = E6[phi^6]/120
    d_dir_pure = np.sum(f3 * dP ** 3 * w)             # lem:volcoeff limit
    d5 = 6.0 * E5_p4s / 24.0
    d4 = 6.0 * E4_p2s2 / 4.0
    d3 = 6.0 * E3_s3 / 6.0
    dpsi3 = 6.0 * T_psi3
    d = d_frozen + d5 + d4 + d3 + dpsi3

    wL2 = -b_ls / (2.0 * d) if d != 0 else np.nan
    return dict(gap=ev[1] - ev[0], b_ls=b_ls, b_direct=b_direct,
                d=d, d_frozen=d_frozen, d_dir_pure=d_dir_pure,
                d5=d5, d4=d4, d3=d3, dpsi3=dpsi3, wL2=wL2)


if __name__ == "__main__":
    mu = 1.0
    print("=" * 78)
    print("Reduced quintic d (6th-order LS) on the graded cell; ass:subcrit needs d>0")
    print("d = d_frozen + d5 + d4 + d3 + dpsi3;  dpsi3 <= 0 always (mirror of fbsign)")
    print("=" * 78)

    print("\n(1) Rectangle sweep (grid of section D: NX=10, N=24, linear grading)")
    print(f"{'nu':>7} {'Khat':>6} {'b':>9} {'|b_chk|':>8} {'d':>10} {'d_frozen':>10} "
          f"{'d5':>9} {'d4':>9} {'d3':>9} {'dpsi3':>9} {'w_L^2':>8}  sign")
    for Lam in [1.0, 2.0, 5.0, 10.0]:
        cell = GradedCell(mu, Lam, NX=10, N=24)
        for Khat in [0.10, 0.20, 0.30]:
            lc = cell.critical_l1(Khat)
            if not np.isfinite(lc):
                print(f"{nu_of(Lam,mu):7.3f} {Khat:6.2f}   (no critical l1)")
                continue
            r = quintic_data(cell, lc, Khat)
            chk = abs(r['b_ls'] - r['b_direct']) / abs(r['b_direct'])
            print(f"{nu_of(Lam,mu):7.3f} {Khat:6.2f} {r['b_ls']:9.4f} {chk:8.1e} "
                  f"{r['d']:10.4f} {r['d_frozen']:10.4f} {r['d5']:9.3f} "
                  f"{r['d4']:9.3f} {r['d3']:9.3f} {r['dpsi3']:9.3f} "
                  f"{r['wL2']:8.4f}  {'d>0' if r['d']>0 else 'd<0'}")

    print("\n(2) Resolution study (nu=0.417, Khat=0.15)")
    print(f"{'NX':>4} {'N':>4} {'gap':>10} {'b':>9} {'d':>10} {'d_frozen':>10} "
          f"{'dpsi3':>9} {'w_L^2':>8}")
    for NX, N in [(10, 24), (12, 28), (14, 32), (16, 36)]:
        cell = GradedCell(mu, 5.0, NX=NX, N=N)
        lc = cell.critical_l1(0.15)
        r = quintic_data(cell, lc, 0.15)
        print(f"{NX:4d} {N:4d} {r['gap']:10.2e} {r['b_ls']:9.4f} {r['d']:10.4f} "
              f"{r['d_frozen']:10.4f} {r['dpsi3']:9.3f} {r['wL2']:8.4f}")
