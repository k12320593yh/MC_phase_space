import contraction as lc
import Lorentz_Transformation as lt
import LorentzGenerators as lg
import numpy as np

a = np.linspace(0,26,27)
b = a.reshape(3,3,3)
print(b)
print(b.ravel())