# Stage 3 -- the open computation of numerical.tex sub:numb:
# the subcritical sign b(nu, Khat) on the GRADED substrate.
#
# This replaces the unconverged Stage-2 of creasing_buckling.py. The Stage-2 bug
# (documented in handoff.md and at the foot of creasing_buckling.py) was mixing a
# STRONG-form (collocation) critical null vector into a WEAK-form (Galerkin)
# nonlinear assembly: the two formulations carry inconsistent spurious-mode
# content, so the second-order source grew with N and b_psi never converged.
#
# The fix, exactly as prescribed in the handoff / TODO block of numerical.tex:
#   (1) ONE consistent weak (Galerkin) formulation throughout.  The critical mode
#       phi is the kernel of the SAME symmetric Galerkin operator Q that is then
#       used to solve the Lyapunov-Schmidt second-order problem.
#   (2) Real nodal DOFs on a genuine 2D grid (X along the free face x Y in depth),
#       so the crease harmonic and its second harmonic are represented by real
#       cos/sin content with NO complex-conjugation ambiguity in the load.
#   (3) A GRADED base state lambda1 = lambda1(X) of amplitude Khat.  On the flat
#       (Khat=0) half-space the Biot threshold is scale-free: Q has a kernel for
#       EVERY wavenumber, the second harmonic sits on the same neutral curve, and
#       the LS solve Q chi = -source is resonant -> b is undefined
#       (Remark rem:flatdegenerate).  The grading lifts the degeneracy: the
#       critical mode becomes a discrete kernel and Q is invertible on V^perp, so
#       chi -- and hence b -- is well defined.  This script verifies that and
#       returns the sign of b over the (nu, Khat) rectangle.
#
# Model.  2D plane-strain compressible neo-Hookean of the paper,
#   Psi(F) = (mu/2)(F:F - 2) + f(J),  f(J) = -mu log J + Lam log^2 J.
# Base state graded slowly along the free face,
#   lambda1(X) = lambda1bar * (1 + Khat cos(kX)),   lambda2(X) from self-equilibrium
#   mu lambda2^2 + f'(J)J = 0  (traction-free surface, base_l2).
# The crease cell is one wavelength X in [0, 2pi/k) (k=1), periodic; depth Y<0 is
# mapped Chebyshev (free surface Y=0 natural, decay Y->-inf essential).
#
# Weak operator (single, consistent):
#   Q[u,v] = INT_cell  Mat4_pq(X)  G_p(u) G_q(v),
#   gradient 4-vector ordering p=(i,J): 0=dX uX, 1=dY uX, 2=dX uY, 3=dY uY,
#   Mat4(X) = mu I4 + f''(J0) cv(X) (x) cv(X) + f'(J0)J0/J0 ... = mu I4
#             + fpp cv cv^T + fp Mmat,   cv = (lambda2,0,0,lambda1) = cof F0 flattened,
#   Mmat the polarised-determinant matrix.  (Same density as creasing_buckling.py,
#   now X-dependent through the grading.)
#
# Lyapunov-Schmidt (energy/volumetric form; deviatoric is quadratic, drops out):
#   source   Sfield = f'' p (Mmat gphi) + (f'' dphi + 1/2 f''' p^2) cv,
#            p = cv:gphi,  dphi = det gphi,
#   bordered solve   [[Q, M phi],[ (M phi)^T, 0]] [chi; xi] = [-R; 0],   R_q = G_q^T(omega Sfield),
#   feedback Rchi = INT g(chi).Sfield,
#   E4    = INT [ 1/2 f'' dphi^2 + 1/2 f''' p^2 dphi + 1/24 f'''' p^4 ],
#   b     = 4 E4 + 2 Rchi.
# Sign of b is independent of the (arbitrary) normalisation of phi: E4 ~ phi^4 and
# Rchi ~ phi^2 * chi ~ phi^4, so rescaling phi rescales b by a positive factor.

import numpy as np
from surface_instability import base_l2, lambda_star
from creasing_buckling import cheb, depth_ops, clencurt, fderivs

# polarised-determinant matrix:  A^T Mmat B = M(A,B) = A00 B11 + B00 A11 - A01 B10 - B01 A10
MMAT = np.array([[0, 0, 0, 1.0],
                 [0, 0, -1.0, 0],
                 [0, -1.0, 0, 0],
                 [1.0, 0, 0, 0]])


class GradedCell:
    """Consistent weak operator + LS coefficient b on the monotonically graded cell.

    The free face is the coordinate X in [-Lx, Lx] (arc length along the corner
    face, rescaled).  The base compression is graded MONOTONICALLY through the Biot
    threshold, lambda1(X) = lambda1bar (1 + Khat * X/Lx), so the local stretch
    crosses lambda_star at a single turning point X=0.  Unlike a periodic
    modulation -- which has two crossings per period and re-imposes wavenumber
    structure, leaving the scale-free Biot degeneracy only weakly split -- the
    monotone grading produces a single LOCALISED neutral mode: a genuine discrete,
    isolated eigenvalue with an O(1) spectral gap, exactly the regularisation of
    Remark rem:flatdegenerate.  X is discretised by Chebyshev (mode decays away
    from the turning point; homogeneous Dirichlet at +-Lx), depth Y by mapped
    Chebyshev (free surface natural, decay at -inf essential)."""

    def __init__(self, mu, Lam, NX=20, N=40, Lmap=4.0, Lx=6.0):
        self.mu, self.Lam = mu, Lam
        self.NX, self.N, self.Lmap, self.Lx = NX, N, Lmap, Lx
        DY, DY2, Y, xi = depth_ops(N, Lmap)
        self.DY, self.Y = DY, Y
        n = N + 1
        self.n = n

        # Chebyshev in X on [-Lx, Lx]
        Dc, xc = cheb(NX)
        self.Xg = Lx * xc                       # (NX+1,) X nodes, descending
        DX = Dc / Lx
        nx = NX + 1
        self.nx = nx
        wXc = clencurt(NX) * Lx                  # Clenshaw-Curtis on [-Lx,Lx]
        self.wX = wXc

        wcc = clencurt(N)
        dYdxi = -2.0 * Lmap / (1.0 - xi) ** 2
        wY = wcc * np.abs(dYdxi)
        wY[0] = 0.0                       # far node (eta=inf): pinned, zero weight
        self.wY = wY
        # node-wise weight, ordering (iX, iY) flattened, length nx*n
        self.omega = np.outer(wXc, wY).ravel()

        # gradient block operators on a single scalar component (nx*n square)
        I_n = np.eye(n)
        I_NX = np.eye(nx)
        self.DXf = np.kron(DX, I_n)
        self.DYf = np.kron(I_NX, DY)

        # global DOF layout: [uX nodal | uY nodal], each nx*n
        m = nx * n
        self.m = m
        # gradient 4-vector ordering p=(i,J): comp i = p//2, deriv J = p%2 (0=X,1=Y)
        self.Dop = [self.DXf, self.DYf]          # derivative operators (m x m)
        self.comp_of = [0, 0, 1, 1]              # which component each g_p reads
        self.der_of = [0, 1, 0, 1]               # which derivative each g_p applies

        # essential BC: drop far-Y node (iY==0) and X-end nodes (iX in {0,NX})
        drop = set()
        for iX in range(nx):
            drop.add(iX * n + 0)                                # Y far
        for iX in (0, NX):                                      # X ends (decay)
            for iY in range(n):
                drop.add(iX * n + iY)
        drop = drop | {m + d for d in drop}
        self.keep = np.array([i for i in range(2 * m) if i not in drop])

        # full node weight for both components (mass), length 2m
        self.omega2 = np.concatenate([self.omega, self.omega])

    # -- per-node base-state material data for a given load & grading -------------
    def material(self, l1bar, Khat):
        mu, Lam = self.mu, self.Lam
        l1 = l1bar * (1.0 + Khat * (self.Xg / self.Lx))        # monotone grading
        l2 = np.array([base_l2(v, mu, Lam) for v in l1])
        J0 = l1 * l2
        fp, fpp, fppp, fpppp = fderivs(J0, mu, Lam)            # each (nx,)
        cv = np.zeros((self.nx, 4))
        cv[:, 0] = l2          # cof F0 = diag(l2, l1) -> flat (l2,0,0,l1)
        cv[:, 3] = l1
        return dict(l1=l1, l2=l2, J0=J0, fp=fp, fpp=fpp, fppp=fppp,
                    fpppp=fpppp, cv=cv)

    def _tile(self, perX):
        """expand a per-X array (nx,...) to per-node (nx*n,...) by repeating in Y."""
        return np.repeat(perX, self.n, axis=0)

    def mat4_nodes(self, mat):
        """4x4 Mat4 at every node -> array (4,4,nx*n)."""
        n = self.n
        fpp_nd = self._tile(mat['fpp'])
        fp_nd = self._tile(mat['fp'])
        cv_nd = self._tile(mat['cv'])                          # (NX*n,4)
        M4 = (self.mu * np.eye(4)[:, :, None]
              + fpp_nd[None, None, :] * cv_nd.T[:, None, :] * cv_nd.T[None, :, :]
              + MMAT[:, :, None] * fp_nd[None, None, :])
        return M4                                              # (4,4,NX*n)

    def gradapply(self, p, vec):
        """apply gradient operator G_p to a global DOF vector (length 2m)."""
        m = self.m
        comp = vec[self.comp_of[p] * m:(self.comp_of[p] + 1) * m]
        return self.Dop[self.der_of[p]] @ comp

    def stiffness(self, mat):
        """symmetric weak operator Q (2m x 2m), assembled in m-sized blocks."""
        M4 = self.mat4_nodes(mat)
        m = self.m
        Q = np.zeros((2 * m, 2 * m))
        for p in range(4):
            ip, Dp = self.comp_of[p], self.Dop[self.der_of[p]]
            for q in range(4):
                iq, Dq = self.comp_of[q], self.Dop[self.der_of[q]]
                w = self.omega * M4[p, q]                      # (m_node,) = (NX*n,)
                blk = Dp.T @ (w[:, None] * Dq)                 # (m x m)
                Q[ip * m:(ip + 1) * m, iq * m:(iq + 1) * m] += blk
        return 0.5 * (Q + Q.T)

    # -- critical load: smallest eigenvalue of Q crosses zero --------------------
    def min_eig(self, l1bar, Khat):
        from scipy.linalg import eigvalsh
        Q = self.stiffness(self.material(l1bar, Khat))
        Qk = Q[np.ix_(self.keep, self.keep)]
        return float(eigvalsh(Qk, subset_by_index=[0, 0])[0])

    def critical_l1(self, Khat, lo=0.40, hi=0.88, npts=13):
        from scipy.optimize import brentq
        f = lambda l: self.min_eig(l, Khat)
        # scan for a sign change (Q pos-def above threshold, indefinite below)
        grid = np.linspace(hi, lo, npts)
        vals = [f(g) for g in grid]
        for i in range(len(grid) - 1):
            if np.isfinite(vals[i]) and vals[i] * vals[i + 1] < 0:
                return brentq(f, grid[i], grid[i + 1], xtol=1e-5)
        return np.nan

    # -- the Landau coefficient b at a given (graded) critical state -------------
    def landau_b(self, l1bar, Khat):
        mat = self.material(l1bar, Khat)
        Q = self.stiffness(mat)
        keep = self.keep
        Qk = Q[np.ix_(keep, keep)]

        # critical mode phi = kernel of the SAME weak operator (lowest 2 only)
        from scipy.linalg import eigh
        ev, V = eigh(Qk, subset_by_index=[0, 1])
        gap = ev[1] - ev[0]                                    # spectral gap above kernel
        phi_k = V[:, 0]
        phi = np.zeros(2 * self.m)
        phi[keep] = phi_k

        # physical gradient of phi (real), 4 x (NX*n)
        gphi = np.array([self.gradapply(p, phi) for p in range(4)])
        scale = np.sqrt(np.max(np.sum(gphi ** 2, axis=0)))
        phi = phi / scale
        gphi = gphi / scale

        cv_nd = self._tile(mat['cv']).T                       # (4, NX*n)
        fpp_nd = self._tile(mat['fpp'])
        fppp_nd = self._tile(mat['fppp'])
        fpppp_nd = self._tile(mat['fpppp'])

        p_ = np.sum(cv_nd * gphi, axis=0)                     # cv:gphi
        dphi = gphi[0] * gphi[3] - gphi[1] * gphi[2]          # det gphi
        Mgphi = MMAT @ gphi                                   # (4, NX*n)

        # second-order source Sfield (4, NX*n), all real
        scal = fpp_nd * dphi + 0.5 * fppp_nd * p_ ** 2
        Sfield = fpp_nd * p_[None, :] * Mgphi + scal[None, :] * cv_nd

        # load vector R_q = G_q^T (omega Sfield_q)
        R = np.zeros(2 * self.m)
        for q in range(4):
            iq, Dq = self.comp_of[q], self.Dop[self.der_of[q]]
            R[iq * self.m:(iq + 1) * self.m] += Dq.T @ (self.omega * Sfield[q])

        # bordered LS solve on V^perp: [[Q, Mphi],[Mphi^T,0]][chi;xi] = [-R;0]
        Mphi = self.omega2 * phi                              # mass-weighted kernel
        Mk = Mphi[keep]
        nk = len(keep)
        B = np.zeros((nk + 1, nk + 1))
        B[:nk, :nk] = Qk
        B[:nk, nk] = Mk
        B[nk, :nk] = Mk
        rhs = np.zeros(nk + 1)
        rhs[:nk] = -R[keep]
        sol = np.linalg.solve(B, rhs)
        chi = np.zeros(2 * self.m)
        chi[keep] = sol[:nk]

        gchi = np.array([self.gradapply(p, chi) for p in range(4)])
        Rchi = np.sum(np.sum(gchi * Sfield, axis=0) * self.omega)

        e4 = (0.5 * fpp_nd * dphi ** 2
              + 0.5 * fppp_nd * p_ ** 2 * dphi
              + (1.0 / 24.0) * fpppp_nd * p_ ** 4)
        E4 = np.sum(e4 * self.omega)
        bdir = 4.0 * np.sum(0.5 * fpp_nd * dphi ** 2 * self.omega)

        b = 4.0 * E4 + 2.0 * Rchi
        # normalisation-invariant sign diagnostic: both 4E4 and 2Rchi scale as
        # phi^4, so their ratio is independent of the (arbitrary) mode amplitude;
        # b < 0  <=>  ratio < -1  (feedback overwhelms the adverse frozen term).
        ratio = (2.0 * Rchi) / (4.0 * E4) if E4 != 0 else np.nan

        # --- rotation-condition diagnostics (Lemma lem:rot / finding #4) -------
        # On the SAME critical mode, the second variation is
        #   <phi, L phi> = mu |grad phi|^2 + f''(J0) (dJ)^2 + 2 f'(J0) det(grad phi),
        # with the first-order volume change dJ = cof F0 : grad phi  (= p_ here)
        # and det(grad phi) = dphi.  The reviewer's worry is that the LAST term
        # (prefactor 2 f' ~ |log J0|/J0, more singular than f'' ~ |log J0|) is
        # neither bounded nor sign-controlled.  We measure all three pieces on the
        # actual mode, plus the dimensionless ROTATION DEFECT
        #   rho_rot = ||dJ||_2 / ||grad phi||_2   (amplitude-invariant),
        # small <=> the critical mode is asymptotically volume-preserving.
        fp_nd = self._tile(mat['fp'])
        w = self.omega
        grad_l2 = np.sqrt(np.sum(np.sum(gphi ** 2, axis=0) * w))
        dJ_l2 = np.sqrt(np.sum(p_ ** 2 * w))
        rho_rot = dJ_l2 / grad_l2                          # rotation defect (invariant)
        T_shear = np.sum(self.mu * np.sum(gphi ** 2, axis=0) * w)
        T_vol2 = np.sum(fpp_nd * p_ ** 2 * w)              # f''(dJ)^2   term
        T_det = np.sum(2.0 * fp_nd * dphi * w)             # 2 f' det    term (the omitted one)
        coen = T_shear + T_vol2 + T_det                    # = <phi, L phi>
        J0min = float(np.min(mat['J0']))
        # lem:rot predicts (at unit gradient) rho_rot = O(J0 / sqrt|log J0|):
        rate = rho_rot * np.sqrt(abs(np.log(J0min))) / J0min
        # how much the "omitted" 2 f' det term dominates the kept f''(dJ)^2 term:
        det_ratio = T_det / T_vol2 if T_vol2 != 0 else np.nan
        return dict(b=b, bdir=bdir, b_E4=4 * E4, b_feedback=2 * Rchi,
                    ratio=ratio, kernel_eig=ev[0], gap=gap,
                    rho_rot=rho_rot, coen=coen, T_shear=T_shear, T_vol2=T_vol2,
                    T_det=T_det, det_ratio=det_ratio, J0min=J0min, rate=rate)


def nu_of(Lam, mu):
    """plane-strain effective Poisson ratio."""
    return Lam / (2.0 * (Lam + mu))


if __name__ == "__main__":
    mu = 1.0
    print("=" * 72)
    print("Stage 3: consistent weak graded-substrate Landau coefficient b(nu,Khat)")
    print("=" * 72)

    # (A) Mechanism check: the homogeneous (Khat=0) cell is already indefinite at
    # the Biot stretch lambda_star (min eig < 0) -- the flat half-space is unstable
    # at/below onset, with NO discrete critical mode (the scale-free degeneracy of
    # Remark rem:flatdegenerate).  The grading (B) is what selects a discrete mode.
    print("\n(A) Homogeneous (Khat=0) cell at l1=lambda_star: operator already")
    print("    indefinite (min eig<0), no discrete critical mode -- degenerate")
    print(f"{'Lam':>6} {'nu':>7} {'lambda_star':>12} {'min_eig(Q)':>14}")
    for Lam in [1.0, 5.0]:
        cell = GradedCell(mu, Lam, NX=12, N=28)
        ls = lambda_star(mu, Lam)
        me = cell.min_eig(ls, 0.0)
        print(f"{Lam:6.1f} {nu_of(Lam,mu):7.4f} {ls:12.6f} {me:14.3e}")

    # (B) the degeneracy lift: spectral gap above the localised kernel vs Khat
    print("\n(B) Degeneracy lift -- spectral gap above the critical (localised) mode")
    print("    monotone grading -> discrete kernel (eig~0) isolated by an O(1) gap")
    Lam = 5.0
    cell = GradedCell(mu, Lam, NX=12, N=28)
    for Khat in [0.05, 0.10, 0.20, 0.30]:
        lc = cell.critical_l1(Khat)
        if not np.isfinite(lc):
            print(f"   Khat={Khat:4.2f}  (threshold outside scan window)")
            continue
        r = cell.landau_b(lc, Khat)
        print(f"   Khat={Khat:4.2f}  l1bar={lc:7.4f}  kernel_eig={r['kernel_eig']:+.2e}"
              f"  gap={r['gap']:.3e}")

    # (C) convergence of the SIGN diagnostic (the Stage-2 quantity that DIVERGED).
    # The absolute b carries the arbitrary mode-amplitude normalisation (here
    # max|grad phi|=1); as the mesh resolves the increasingly localised mode the
    # L4 integrals concentrate and |b| drifts.  The NORMALISATION-INVARIANT ratio
    # r = b_feedback/b_E4 (both ~ phi^4) is what converges, and b<0 <=> r<-1.
    print("\n(C) Convergence of the sign diagnostic ratio r = b_feedback/b_E4")
    print("    (the Stage-2 failure mode is now consistent; b<0 iff r<-1)")
    Lam, Khat = 5.0, 0.15
    print(f"    Lam={Lam}, nu={nu_of(Lam,mu):.3f}, Khat={Khat}")
    print(f"{'NX':>4} {'N':>4} {'ratio r':>10} {'sign(b)':>9} {'gap':>9}")
    for NX, N in [(10, 24), (12, 28), (14, 32), (16, 36)]:
        cell = GradedCell(mu, Lam, NX=NX, N=N)
        lc = cell.critical_l1(Khat)
        r = cell.landau_b(lc, Khat)
        print(f"{NX:4d} {N:4d} {r['ratio']:10.4f} {'b<0' if r['b']<0 else 'b>0':>9}"
              f" {r['gap']:9.2e}")

    # (D) phase diagram: invariant sign ratio r=b_feedback/b_E4 over (nu, Khat);
    #     b<0 (subcritical) wherever r<-1.
    print("\n(D) Sign ratio r = b_feedback/b_E4 over the rectangle (b<0 iff r<-1)")
    Lams = [1.0, 2.0, 5.0, 10.0]
    Khats = [0.10, 0.20, 0.30]
    header = "  nu \\ Khat " + "".join(f"{kh:>10.2f}" for kh in Khats)
    print(header)
    table = {}
    for Lam in Lams:
        row = f"  {nu_of(Lam,mu):8.3f}  "
        cell = GradedCell(mu, Lam, NX=10, N=24)
        for Khat in Khats:
            lc = cell.critical_l1(Khat)
            if not np.isfinite(lc):
                row += f"{'--':>10}"
                continue
            r = cell.landau_b(lc, Khat)
            table[(Lam, Khat)] = r['b']
            row += f"{r['ratio']:10.3f}"
        print(row)
    signs = [np.sign(v) for v in table.values()]
    print("\n  sign(b):", "uniform" if len(set(signs)) == 1 else "mixed",
          "->", "subcritical (b<0) everywhere" if all(s < 0 for s in signs) else
          ("supercritical (b>0)" if all(s > 0 for s in signs) else "geometry-dependent"))

    # (E) Rotation-condition diagnostic (finding #4): on the SAME critical mode,
    #     report the dimensionless rotation defect rho_rot = ||dJ||/||grad phi||
    #     and the second-variation decomposition.  The claim under scrutiny is that
    #     the "omitted" 2 f' det term does NOT dominate and that the mode is nearly
    #     volume-preserving (rho_rot small).  At the bifurcation J0>0 is finite, so
    #     the second variation 'coen' is finite -- the rigorous core of lem:rot.
    print("\n(E) Rotation condition on the critical mode (finding #4):")
    print("    rho_rot=||dJ||/||grad phi||  (small => mode ~ volume-preserving);")
    print("    coen=<phi,L phi> finite at the bifurcation; |T_det/T_vol2|<~1 means")
    print("    the 'omitted' 2f'det term does NOT dominate the kept f''(dJ)^2 term.")
    print(f"{'nu':>7} {'Khat':>6} {'J0min':>8} {'rho_rot':>9} {'coen':>10}"
          f" {'T_det/T_vol2':>13}")
    for Lam in [1.0, 5.0]:
        cell = GradedCell(mu, Lam, NX=12, N=28)
        for Khat in [0.10, 0.20, 0.30]:
            lc = cell.critical_l1(Khat)
            if not np.isfinite(lc):
                continue
            r = cell.landau_b(lc, Khat)
            print(f"{nu_of(Lam,mu):7.3f} {Khat:6.2f} {r['J0min']:8.4f}"
                  f" {r['rho_rot']:9.4f} {r['coen']:10.3e} {r['det_ratio']:13.4f}")

    # (F) Asymptotic-rate probe of lem:rot: drive the base compression DEEPER
    #     (fundamental J0 -> 0) at fixed grading, tracking the softest localised
    #     mode (lowest eigenvector of the symmetric weak operator).  lem:rot
    #     predicts rho_rot = O(J0/sqrt|log J0|), i.e. the invariant
    #     rate = rho_rot*sqrt|log J0|/J0 stays BOUNDED as J0 -> 0.
    print("\n(F) Asymptotic-rate probe (lem:rot): deeper compression, J0->0.")
    print("    rate = rho_rot*sqrt|log J0|/J0 bounded <=> rotation rate holds.")
    Lam, Khat = 5.0, 0.20
    cell = GradedCell(mu, Lam, NX=12, N=28)
    lc = cell.critical_l1(Khat)
    print(f"    Lam={Lam}, nu={nu_of(Lam,mu):.3f}, Khat={Khat}, l1c={lc:.4f}")
    print(f"{'l1bar':>8} {'J0min':>8} {'rho_rot':>9} {'rate':>9} {'T_det/T_vol2':>13}")
    for frac in [1.00, 0.85, 0.70, 0.55, 0.40]:
        l1b = lc * frac
        r = cell.landau_b(l1b, Khat)
        print(f"{l1b:8.4f} {r['J0min']:8.4f} {r['rho_rot']:9.4f}"
              f" {r['rate']:9.4f} {r['det_ratio']:13.4f}")
