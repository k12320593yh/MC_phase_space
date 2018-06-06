import sympy as sym
import numpy as np
import LorentzGenerators as lg
import Lorentz_Transformation as lt

x = sym.Symbol('x')
y = sym.Symbol('y')
a = 3.14
b = 2.71


aa = sym.solve([x**2,y**2],[x,y])
print(aa)