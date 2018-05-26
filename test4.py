import LorentzGenerators as lg
import Lorentz_Transformation as lt
import numpy as np
np.random.seed(42)

a = lt.Lorentz4vector(components=[1,0,0,0])
b = lt.Lorentz4vector(components=[3,0,0,0])
print(lt.four_vector_contraction(a,b))