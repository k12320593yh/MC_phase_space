import sympy as sp
from sympy.abc import x,y,a,b,c

aa = sp.solve([x**2+a**2],[x])
print(aa)