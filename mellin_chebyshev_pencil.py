# Mellin x Chebyshev spectral solver for the linear-elasticity corner pencil
# (Kondratiev pencil) at the right-angle clamped--free wedge.
#
# Separable ansatz u = r^lambda (a(theta), b(theta)) in the plane-strain Navier
# equations on the wedge theta in (0, pi/2). This turns the corner problem into
# an angular ODE eigenvalue problem -- a QUADRATIC eigenvalue problem in lambda
# -- with Chebyshev collocation in theta. The singular exponents alpha are the
# eigenvalues with Re(lambda) in (0,1).
#
# This is the "Mellin for the corner" half of the spectral programme: scale
# covariance (r -> c r) is diagonalised by r^lambda (= Fourier in log r), and
# the smooth angular profile is captured spectrally. Validated against the
# secular-determinant roots of williams_clampedfree.py.
#
# Navier (plane strain), divided by mu, with kappa1 = 1/(1-2nu):
#   lap u_r - u_r/r^2 - (2/r^2) d_theta u_theta + kappa1 d_r e = 0
#   lap u_theta - u_theta/r^2 + (2/r^2) d_theta u_r + kappa1 (1/r) d_theta e = 0
# e = (1/r) d_r(r u_r) + (1/r) d_theta u_theta.  Substituting u = r^l (a,b):
#   Eq1:  a'' + c(l^2-1) a + [kappa1(l-1) - 2] b' = 0
#   Eq2:  c b'' + (l^2-1) b + [kappa1(l+1) + 2] a' = 0,   c = 1 + kappa1.
# BC theta=0 (clamped): a=b=0.   BC theta=pi/2 (free): sigma_tt = sigma_rt = 0:
#   [g(l+1)+2] a + (g+2) b' = 0,   a' + (l-1) b = 0,   g = kappa1 - 1 = Lame/mu.

import numpy as np
from scipy.linalg import eig


def cheb(N):
    """Chebyshev-Gauss-Lobatto nodes x in [1,-1] and diff matrix (Trefethen)."""
    if N == 0:
        return np.zeros((1, 1)), np.array([1.0])
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.hstack([2.0, np.ones(N - 1), 2.0]) * (-1.0) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1.0 / c) / (dX + np.eye(N + 1))
    D = D - np.diag(D.sum(axis=1))
    return D, x


def pencil_eigs(nu, N=40):
    """Return eigenvalues lambda of the wedge pencil (plane strain)."""
    kappa1 = 1.0 / (1.0 - 2.0 * nu)
    c = 1.0 + kappa1
    g = kappa1 - 1.0

    Dx, x = cheb(N)
    # map x in [1,-1] -> theta = (pi/4)(x+1) in [pi/2, 0]; node 0 is theta=pi/2 (free),
    # node N is theta=0 (clamped).
    s = 4.0 / np.pi          # dx/dtheta
    D = s * Dx
    D2 = D @ D
    n = N + 1
    I = np.eye(n)
    Z = np.zeros((n, n))

    def blk(TLa, TLb, BLa, BLb):
        return np.block([[TLa, TLb], [BLa, BLb]])

    # interior coefficient matrices of  (l^2 M2 + l M1 + M0) phi = 0,  phi=[a;b]
    M2 = blk(c * I, Z, Z, I)
    M1 = blk(Z, kappa1 * D, kappa1 * D, Z)
    M0 = blk(D2 - c * I, -(kappa1 + 2.0) * D,
             (kappa1 + 2.0) * D, c * D2 - I)

    free = 0      # theta = pi/2
    clamp = N     # theta = 0
    ra, rb = free, n + free          # free-end rows (a-eq, b-eq)
    ca, cb = clamp, n + clamp        # clamped-end rows

    for M in (M2, M1, M0):
        for row in (ra, rb, ca, cb):
            M[row, :] = 0.0

    # clamped theta=0:  a_N = 0,  b_N = 0   (in M0)
    M0[ca, clamp] = 1.0
    M0[cb, n + clamp] = 1.0
    # free theta=pi/2, BC1: l*(g a0) + [(g+2)a0 + (g+2)(D b)_0] = 0
    M1[ra, free] = g
    M0[ra, free] = g + 2.0
    M0[ra, n:] = (g + 2.0) * D[free, :]
    # free theta=pi/2, BC2: l*(b0) + [(D a)_0 - b0] = 0
    M1[rb, n + free] = 1.0
    M0[rb, :n] = D[free, :]
    M0[rb, n + free] += -1.0

    # companion linearisation:  A y = l B y,  y = [phi; l phi]
    A = np.block([[Z2 := np.zeros((2 * n, 2 * n)), np.eye(2 * n)],
                  [-M0, -M1]])
    B = np.block([[np.eye(2 * n), np.zeros((2 * n, 2 * n))],
                  [np.zeros((2 * n, 2 * n)), M2]])
    w = eig(A, B, right=False)
    w = w[np.isfinite(w)]
    return w


def smallest_alpha(nu, N=40, tol=1e-6):
    w = pencil_eigs(nu, N)
    cand = [z.real for z in w if abs(z.imag) < tol and tol < z.real < 1.0 - 1e-9]
    return min(cand) if cand else np.nan


if __name__ == "__main__":
    from williams_clampedfree import smallest_root

    print("Mellin x Chebyshev pencil  vs.  secular determinant")
    print(f"{'nu':>6} {'alpha (spectral)':>17} {'alpha (secular)':>16} {'|diff|':>10}")
    for nu in [0.1, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.49]:
        a_spec = smallest_alpha(nu, N=40)
        roots = smallest_root(3.0 - 4.0 * nu)
        a_sec = min(roots) if roots else np.nan
        diff = abs(a_spec - a_sec)
        print(f"{nu:6.3f} {a_spec:17.6f} {a_sec:16.6f} {diff:10.2e}")

    # Spectral (exponential) convergence in the number of angular modes N,
    # the concrete advantage over a graded-mesh FEM (algebraic in DOFs).
    print("\nSpectral convergence at nu=0.3 (alpha_ref from N=64):")
    a_ref = smallest_alpha(0.3, N=64)
    print(f"{'N':>4} {'alpha':>14} {'|alpha-alpha_ref|':>20}")
    for N in [6, 8, 12, 16, 20, 28, 36]:
        a = smallest_alpha(0.3, N=N)
        print(f"{N:4d} {a:14.10f} {abs(a - a_ref):20.2e}")
