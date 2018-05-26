import LorentzGenerators as lg
import Lorentz_Transformation as lt
import numpy as np

p2 = lt.Lorentz4vector(components=[85,84,12,4])
p2.lorentz_rot_toz()
p2.zrapidity()
p2.z_boost_to_rest_frame()
print(p2.get_4_vector())
print(p2.rotmat)
print(p2.boostmat)

print(np.matmul(p2.boostmat,p2.rotmat))