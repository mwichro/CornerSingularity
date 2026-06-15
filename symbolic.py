import sympy as sp
w = sp.symbols('w')
J0,D1,D2,D3,D4 = sp.symbols('J0 D1 D2 D3 D4')
mu,lam = sp.symbols('mu_e lambda_e', positive=True)
Jsym = sp.Function('f')
# J(w) up to w^4 (general d); for d=2 D3=D4=0
Jw = J0 + D1*w + D2*w**2 + D3*w**3 + D4*w**4
f = lambda J: -mu*sp.log(J) + lam*sp.log(J)**2
g = f(Jw)
E3 = sp.diff(g,w,3).subs(w,0)
E4 = sp.diff(g,w,4).subs(w,0)
print('E3 (=2*c2):'); sp.pprint(sp.simplify(E3))
print('E4 (=6*b):'); sp.pprint(sp.simplify(E4))
# Faa di Bruno generic check (replace f-derivs by symbols)
fp,fpp,fppp,fpppp = sp.symbols('fp fpp fppp fpppp')
fg = sp.Function('F')
gg = fg(Jw)
import sympy
E3g = sp.diff(gg,w,3).subs(w,0)
E4g = sp.diff(gg,w,4).subs(w,0)
sub = {sp.Derivative(fg(J0),(J0,1)):fp, sp.Derivative(fg(J0),(J0,2)):fpp, sp.Derivative(fg(J0),(J0,3)):fppp, sp.Derivative(fg(J0),(J0,4)):fpppp}
print('E3 generic:', sp.expand(E3g.subs(fg(J0),0).doit()) if False else sp.simplify(E3g))
print('E4 generic:', sp.simplify(E4g))

# --- quintic d = E6/120, plus frozen-psi (D1->0) direct coefficients, d=2 ---
# c2 = E3/2, b = E4/6, d = E6/120  ; for d=2 set D3=D4=0
Jw2 = J0 + D1*w + D2*w**2
gg2 = fg(Jw2)
E5g = sp.simplify(sp.diff(gg2,w,5).subs(w,0))
E6g = sp.simplify(sp.diff(gg2,w,6).subs(w,0))
print('E5 generic (d=2):', E5g)
print('E6 generic (d=2):', E6g)

# direct (frozen-psi) leading coefficients: D1 small -> keep D1^0 terms
b_dir = sp.Rational(1,6)*E6g.subs(D1,0)*0 + sp.Rational(1,6)*(12*D2**2*sp.Derivative(fg(J0),(J0,2)))  # E4 at D1=0
d_dir = sp.Rational(1,120)*(120*D2**3*sp.Derivative(fg(J0),(J0,3)))                                    # E6 at D1=0
print('b_direct (D1=0):', sp.simplify(b_dir), '   => 2*D2^2*f\'\'')
print('d_direct (D1=0):', sp.simplify(d_dir), '   => D2^3*f\'\'\'')

# f_log signs in compression J0<1:  f''>0, f'''<0
fp2  = sp.simplify(sp.diff(f(J0),J0,2)*J0**2)   # = mu + 2lam - 2lam log J0  > 0  for J0<1
fp3  = sp.simplify(sp.diff(f(J0),J0,3)*J0**3)   # = 4lam log J0 - 6lam - 2mu  < 0 for J0<1
print("f''*J0^2 =", fp2, "  (>0 for J0<1 -> b_direct>0)")
print("f'''*J0^3 =", fp3, "  (<0 for J0<1 -> d_direct<0)")

