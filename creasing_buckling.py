# Spectral (Chebyshev-in-depth) buckling solver for the compressed compressible
# neo-Hookean half-space, and the weakly-nonlinear cubic (Landau) coefficient b
# of the surface mode. This is the buckling + Lyapunov-Schmidt building block for
# the corner coefficient b(nu, Khat) of Proposition prop:bcount.
#
# Stage 1 (this file, validated first): the LINEAR buckling operator, whose
# critical stretch must match the independent Stroh solver surface_instability.py.
# Stage 2: the second-order field and the cubic coefficient via the exact 2D
# volumetric structure  J = J0 + D1[u] + D2[u]  (D1 linear, D2 = det grad u).
#
# Geometry: body Y<0, traction-free surface Y=0, compression along X (stretch l1).
# Incremental harmonic u = Re[ U(Y) e^{i k X} ];  depth Y mapped to Chebyshev.

import numpy as np
from surface_instability import base_l2, tangent


def cheb(N):
    if N == 0:
        return np.zeros((1, 1)), np.array([1.0])
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2.0, np.ones(N - 1), 2.0]) * (-1.0) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1.0 / c) / (dX + np.eye(N + 1))
    D = D - np.diag(D.sum(axis=1))
    return D, x


def depth_ops(N, Lmap):
    """Chebyshev operators for eta = -Y in [0, inf): eta = Lmap(1+xi)/(1-xi)."""
    D, xi = cheb(N)
    # node xi[0]=1 -> eta=inf (far);  xi[N]=-1 -> eta=0 (surface)
    a = (1.0 - xi) ** 2 / (2.0 * Lmap)     # dxi/deta, so d/deta = a * D
    Deta = np.diag(a) @ D
    DY = -Deta                              # d/dY
    DY2 = Deta @ Deta                       # d^2/dY^2
    eta = Lmap * (1.0 + xi) / (1.0 - xi)
    eta[0] = np.inf
    Y = -eta
    return DY, DY2, Y, xi


def buckling_matrix(l1, mu, Lam, N, k=1.0, Lmap=4.0):
    """Discrete incremental-equilibrium operator with surface/decay BCs."""
    l2 = base_l2(l1, mu, Lam)
    L = tangent(l1, l2, mu, Lam)
    DY, DY2, Y, xi = depth_ops(N, Lmap)
    n = N + 1
    I = np.eye(n)
    ik = 1j * k

    def blk(i, kk):
        return (L[i, 0, kk, 0] * (ik * ik) * I
                + (L[i, 0, kk, 1] + L[i, 1, kk, 0]) * ik * DY
                + L[i, 1, kk, 1] * DY2)

    A = np.block([[blk(0, 0), blk(0, 1)], [blk(1, 0), blk(1, 1)]]).astype(complex)

    surf, far = N, 0
    ra, rb = surf, n + surf          # surface rows (traction), per component
    fa, fb = far, n + far            # far rows (Dirichlet decay)

    # surface traction t_i = L_i2k1 (ik) u_k + L_i2k2 DY u_k = 0
    for i, row in ((0, ra), (1, rb)):
        A[row, :] = 0.0
        for kk, off in ((0, 0), (1, n)):
            A[row, off:off + n] = L[i, 1, kk, 0] * ik * I[surf, :] + L[i, 1, kk, 1] * DY[surf, :]
    # far-field decay: u = 0
    for row, off in ((fa, 0), (fb, n)):
        A[row, :] = 0.0
        A[row, off + far] = 1.0
    return A, L, l2, Y


def smallest_sv(l1, mu, Lam, N=60, k=1.0, Lmap=4.0):
    A, *_ = buckling_matrix(l1, mu, Lam, N, k, Lmap)
    return np.linalg.svd(A, compute_uv=False)[-1]


def critical_stretch(mu, Lam, N=60, Lmap=4.0):
    """Locate l1 minimising the smallest singular value (the buckling stretch)."""
    grid = np.linspace(0.45, 0.70, 60)
    sv = np.array([smallest_sv(l, mu, Lam, N, 1.0, Lmap) for l in grid])
    j = int(np.argmin(sv))
    lo, hi = grid[max(j - 1, 0)], grid[min(j + 1, len(grid) - 1)]
    # golden-section refine
    gr = (np.sqrt(5) - 1) / 2
    c, d = hi - gr * (hi - lo), lo + gr * (hi - lo)
    for _ in range(60):
        if smallest_sv(c, mu, Lam, N, 1.0, Lmap) < smallest_sv(d, mu, Lam, N, 1.0, Lmap):
            hi = d
        else:
            lo = c
        c, d = hi - gr * (hi - lo), lo + gr * (hi - lo)
    return 0.5 * (lo + hi)


if __name__ == "__main__":
    from surface_instability import lambda_star
    mu = 1.0
    print("Linear buckling: Chebyshev-in-depth  vs  Stroh (surface_instability.py)")
    print(f"{'Lam/mu':>8} {'l* (Chebyshev)':>15} {'l* (Stroh)':>12} {'|diff|':>10}")
    for Lam in [1.0, 2.0, 5.0, 10.0, 50.0]:
        lc = critical_stretch(mu, Lam, N=70, Lmap=4.0)
        ls = lambda_star(mu, Lam)
        print(f"{Lam/mu:8.1f} {lc:15.6f} {ls:12.6f} {abs(lc-ls):10.2e}")


# ---------------------------------------------------------------------------
# Stage 2: weakly-nonlinear cubic (Landau) coefficient  b = b_dir + b_psi.
# Energy (weak) form: all cubic/quartic terms are volumetric (deviatoric is
# quadratic), with J = J0 + D1[u] + D2[u] exact in 2D (D1=cofF0:grad u linear,
# D2 = det grad u quadratic). Second variation density:
#   C(A,B) = mu (A:B) + f''(c:A)(c:B) + f' M(A,B),   c = cof F0,
#   M(A,B) = A00 B11 + B00 A11 - A01 B10 - B01 A10   (polarised det).
# LS:  u = w phi + chi,  chi = w^2 chi2,  Q(chi2,v) = -3 T(phi,phi,v),
#   3T(phi,phi,v) = INT[ f''( p M(gphi,gv) + q d_phi ) + 1/2 f''' p^2 q ],
#   p = c:gphi,  q = c:gv,  d_phi = det gphi.
#   E4[phi] = INT[ 1/2 f'' d_phi^2 + 1/2 f''' p^2 d_phi + 1/24 f'''' p^4 ].
#   b = 4 E4[phi] + 2 R[chi2],   R[chi2] = -Q(chi2,chi2) <= 0   (the feedback).

def clencurt(N):
    """Clenshaw-Curtis quadrature weights on Chebyshev-Lobatto nodes (Trefethen)."""
    th = np.pi * np.arange(N + 1) / N
    w = np.zeros(N + 1)
    v = np.ones(N - 1)
    for k in range(1, N // 2 + 1):
        coef = 1.0 if (2 * k == N) else 2.0
        v -= coef * np.cos(2 * k * th[1:N]) / (4 * k * k - 1)
    w[1:N] = 2.0 * v / N
    w[0] = 1.0 / (N * N) if N % 2 == 0 else 1.0 / (N * N - 1)
    w[N] = w[0]
    return w


def fderivs(J0, mu, Lam):
    """f'(J0)J0?? no: return f', f'', f''', f'''' for f_log at J0."""
    fp = -mu / J0 + 2 * Lam * np.log(J0) / J0
    fpp = (mu + 2 * Lam - 2 * Lam * np.log(J0)) / J0 ** 2
    fppp = (4 * Lam * np.log(J0) - 6 * Lam - 2 * mu) / J0 ** 3
    fpppp = (6 * mu + 22 * Lam - 12 * Lam * np.log(J0)) / J0 ** 4
    return fp, fpp, fppp, fpppp


def grad_harm(Uvec, DY, m, n):
    """gradient 4-vector [g00,g01,g10,g11]=[dXuX,dYuX,dXuY,dYuY] for harmonic m."""
    UX, UY = Uvec[:n], Uvec[n:]
    return np.array([1j * m * UX, DY @ UX, 1j * m * UY, DY @ UY])  # (4,n)


def stiffness(m, DY, n, omega, Mat4):
    """Galerkin Q_m (2n x 2n, Hermitian) for harmonic m, weight omega (n,)."""
    # gradient operator blocks: rows g_p (n x 2n)
    I = np.eye(n)
    Z = np.zeros((n, n))
    G = [np.block([[1j * m * I, Z]]),    # g00 = dX uX
         np.block([[DY, Z]]),            # g01 = dY uX
         np.block([[Z, 1j * m * I]]),    # g10 = dX uY
         np.block([[Z, DY]])]            # g11 = dY uY
    W = np.diag(omega)
    Q = np.zeros((2 * n, 2 * n), dtype=complex)
    for p in range(4):
        for q in range(4):
            if Mat4[p, q] != 0.0:
                Q += Mat4[p, q] * (G[p].conj().T @ W @ G[q])
    return Q, G, W


def landau_b(mu, Lam, l1, N=80, Lmap=4.0, NX=8):
    l2 = base_l2(l1, mu, Lam)
    J0 = l1 * l2
    c = np.array([[l2, 0.0], [0.0, l1]])              # cof F0 (2D) = [[l2,0],[0,l1]]
    cv = np.array([c[0, 0], c[0, 1], c[1, 0], c[1, 1]])
    fp, fpp, fppp, fpppp = fderivs(J0, mu, Lam)

    DY, DY2, Y, xi = depth_ops(N, Lmap)
    n = N + 1
    wcc = clencurt(N)
    dYdxi = -2.0 * Lmap / (1.0 - xi) ** 2
    omega = wcc * np.abs(dYdxi)
    omega[0] = 0.0                                    # far node (eta=inf): field pinned to 0

    Mmat = np.array([[0, 0, 0, 1.0], [0, 0, -1.0, 0], [0, -1.0, 0, 0], [1.0, 0, 0, 0]])
    Mat4 = mu * np.eye(4) + fpp * np.outer(cv, cv) + fp * Mmat

    # --- critical mode phi (harmonic 1) from the strong-form operator ---
    A, *_ = buckling_matrix(l1, mu, Lam, N, k=1.0, Lmap=Lmap)
    U, S, Vh = np.linalg.svd(A)
    phi = Vh[-1].conj()                               # null vector -> Uhat_1 (2n,)
    # normalise so max|grad| ~ 1 (sign/scale irrelevant for the SIGN of b)
    gphi = grad_harm(phi, DY, 1, n)
    phi = phi / np.sqrt(np.max(np.abs(gphi) ** 2))
    gphi = grad_harm(phi, DY, 1, n)                  # (4,n) complex, harmonic-1 amp

    # pin far node (essential decay BC) for the chi solves
    keep = np.array([i for i in range(2 * n) if i % n != 0])

    # --- assemble harmonic sources S_m (m=0,2) physically in X, then project ---
    Xg = 2 * np.pi * np.arange(NX) / NX
    # physical gradient of phi: g(X,Y)=phi_hat e^{iX}+c.c. = 2 Re[gphi e^{iX}]
    # build g_phys (4,NX,n) real
    g_phys = np.real(gphi[:, None, :] * np.exp(1j * Xg)[None, :, None]
                     + np.conj(gphi)[:, None, :] * np.exp(-1j * Xg)[None, :, None])
    g00, g01, g10, g11 = g_phys
    p_ph = cv[0] * g00 + cv[1] * g01 + cv[2] * g10 + cv[3] * g11          # c:gphi
    dphi_ph = g00 * g11 - g01 * g10                                       # det gphi
    Mg = np.array([Mmat @ g_phys[:, i, :].reshape(4, -1) for i in range(NX)])  # unused direct
    # Mmat @ gphi physical: (4,NX,n)
    Mg_phys = np.einsum('pq,qXn->pXn', Mmat, g_phys)
    # Sfield (4,NX,n) = f'' p (Mmat gphi) + (f'' dphi + 1/2 f''' p^2) cv
    scal = (fpp * dphi_ph + 0.5 * fppp * p_ph ** 2)                       # (NX,n)
    Sfield = fpp * p_ph[None, :, :] * Mg_phys + scal[None, :, :] * cv[:, None, None]

    # FFT in X to get harmonic-m components of Sfield: S_m[p,Y] = (1/NX) sum_X S e^{-imX}
    Shat = np.fft.fft(Sfield, axis=1) / NX                                # (4,NX,n); index=m
    I = np.eye(n); Z = np.zeros((n, n))

    def Gops(m):
        return [np.block([[1j * m * I, Z]]), np.block([[DY, Z]]),
                np.block([[Z, 1j * m * I]]), np.block([[Z, DY]])]

    def solve_chi(m):
        Qm, _, _ = stiffness(m, DY, n, omega, Mat4)
        # load R_m[v] = INT gv . S_m  ->  vector  sum_p G_p^H W S_m^p
        W = np.diag(omega)
        Sm = Shat[:, m, :]                                                # (4,n)
        Gm = Gops(m)
        Rm = np.zeros(2 * n, dtype=complex)
        for p in range(4):
            Rm += Gm[p].conj().T @ (omega * Sm[p])
        # solve Q_m chi = -R_m with far node pinned to 0
        Qk = Qm[np.ix_(keep, keep)]
        Rk = Rm[keep]
        chik = np.linalg.solve(Qk, -Rk)
        chi = np.zeros(2 * n, dtype=complex); chi[keep] = chik
        return chi, Qm, Rm

    chi0, Q0, R0 = solve_chi(0)
    chi2, Q2, R2 = solve_chi(2)

    # --- R[chi2] feedback (physical, real):  sum over harmonics of R_m . chi_m ---
    # chi_phys = chi0 (real, m=0) + 2 Re[chi2 e^{2iX}]
    # R[chi2] = INT g(chi_phys) . Sfield  (Sfield already built from phi,phi)
    gchi0 = grad_harm(chi0, DY, 0, n)
    gchi2 = grad_harm(chi2, DY, 2, n)
    gchi_phys = (np.real(gchi0)[:, None, :]
                 + np.real(gchi2[:, None, :] * np.exp(2j * Xg)[None, :, None]
                           + np.conj(gchi2)[:, None, :] * np.exp(-2j * Xg)[None, :, None]))
    Rchi = np.sum(np.einsum('pXn,pXn->Xn', gchi_phys, Sfield) * omega[None, :]) \
        * (2 * np.pi / NX)

    # --- E4[phi] (physical) ---
    e4dens = 0.5 * fpp * dphi_ph ** 2 + 0.5 * fppp * p_ph ** 2 * dphi_ph \
        + (1.0 / 24.0) * fpppp * p_ph ** 4
    E4 = np.sum(e4dens * omega[None, :]) * (2 * np.pi / NX)
    # b_dir alone (frozen: drop p and chi):  4 * INT 1/2 f'' dphi^2
    bdir = np.sum(0.5 * fpp * dphi_ph ** 2 * omega[None, :]) * (2 * np.pi / NX) * 4.0

    b = 4.0 * E4 + 2.0 * Rchi
    return dict(b=b.real, bdir=bdir, b_E4=4 * E4.real, b_feedback=2 * Rchi.real,
                J0=J0, l2=l2, fpp=fpp)


if __name__ == "__main__":
    # VALIDATED: the frozen single-mode part b_dir = 2 INT f''(J0) d_phi^2 > 0,
    # confirming the paper's "adverse" frozen sign (Lemma volcoeff(ii)).
    print("\nFrozen single-mode part b_dir = 2 INT f''(J0) d_phi^2  (expect > 0):")
    print(f"{'Lam/mu':>7} {'l*':>9} {'b_dir':>10}")
    for Lam in [1.0, 2.0, 5.0, 10.0]:
        ls = critical_stretch(mu := 1.0, Lam, N=70, Lmap=4.0)
        r = landau_b(1.0, Lam, ls, N=70, Lmap=4.0, NX=8)
        print(f"{Lam:7.1f} {ls:9.5f} {r['bdir']:10.4f}")
    # NOT YET CONVERGED: the Lyapunov-Schmidt feedback term r['b_feedback'].
    # It diverges with N because phi is taken from the STRONG-form (collocation)
    # null space and fed into a WEAK-form (Galerkin) nonlinear assembly; the two
    # formulations' spurious-mode content is inconsistent and pollutes the
    # second-order source. The fix is a single consistent weak formulation:
    # take phi from the Galerkin Q_1 generalised eigenproblem and represent the
    # m=2 harmonic by real (cos/sin) DOFs to remove the conjugation ambiguity.
    print("\n[Stage-2 feedback term b_psi: implemented but NOT converged -- see note]")
