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

