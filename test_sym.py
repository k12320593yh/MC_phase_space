import numpy as np
import sympy as sp
from sympy.abc import x,y,a,b,c

# aa = sp.solve([x**2+a**2],[x])
# print(aa)
# print(np.random.randint(0,2))
x = 10
a = np.random.rand(3,)
b = np.random.rand(3,)
c = np.random.rand(3,)
print(np.vstack((a,b,c)))
# np.save()