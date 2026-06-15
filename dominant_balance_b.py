## confirm the dominant balance of b numerically before writing up Direction 2.


import numpy as np
mu, lam = 1.0, 1.0
def fpp(J):   return (mu+2*lam-2*lam*np.log(J))/J**2
def fppp(J):  return (4*lam*np.log(J)-6*lam-2*mu)/J**3
def fpppp(J): return (6*mu+22*lam-12*lam*np.log(J))/J**4
print(f'{\"J0\":>10} {\"A=fpppp*D1^4\":>14} {\"B=12fppp D1^2 D2\":>16} {\"C=12 fpp D2^2\":>14}  6b')
for J0 in [1e-1,1e-2,1e-3,1e-4,1e-5]:
    D1 = J0/np.sqrt(abs(np.log(J0)))  # barrier-saturating
    D2 = 1.0
    A = fpppp(J0)*D1**4
    B = 12*fppp(J0)*D1**2*D2
    C = 12*fpp(J0)*D2**2
    print(f'{J0:10.0e} {A:14.4f} {B:16.4f} {C:14.2f}  {A+B+C:.2f}')
print('=> C (positive, ~|logJ0|/J0^2) dominates; b ~ 2*fpp*d_phi^2 > 0')

