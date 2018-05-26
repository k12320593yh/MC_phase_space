import LorentzGenerators as lg
import Lorentz_Transformation as lt
# import MonteCarlo as mc
import numpy as np
import scipy.linalg as sp
import inspect as insp

test_4_momentum = lt.Lorentz4vector(name='test4momentum',components=np.array([5,4,0,0]),mass=3)
cnt = 0
for i in range(10):
    test_transformation = lt.LorentzTransFormation(name='test', parameter=np.random.rand(6, ))
    test_4_momentum.lorentz_transformation(test_transformation)
    if not lt.check_on_shell_ness(test_4_momentum.get_4_vector(),test_4_momentum.get_mass()):
        cnt += 1
print(cnt)
# True
# After a certain lorentz transformation the new four vector is still on shell
# The numerical error accumulates really fast. Do not perform consectutive Lorentz Transformatin before eliminating the error.