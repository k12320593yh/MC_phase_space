import contraction as lc
import Lorentz_Transformation as lt
import LorentzGenerators as lg
import numpy as np

# a = np.linspace(0,26,27)
# b = a.reshape(3,3,3)
# print(b)
# print(b.ravel())
k1 = p2 = k2 = p1 = np.array([1,0,0,0])

print(-16*lt.fcc((k1), (p2))*lt.fcc((k2), (p1)))