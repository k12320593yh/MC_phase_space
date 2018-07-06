import vegas as vegas
import numpy as np
import math
import scatteringamplitude
import MonteCarlo as mc
import Lorentz_Transformation as lt



def f(x):
	return (x[0]**2/x[1]**2)*np.log(x[2]+20)*np.exp(x[3]**2)

def g(raw_dot):
    s_sqrt = 500
    p_1 = [s_sqrt / 2, 0, 0, s_sqrt / 2]
    p_2 = [s_sqrt / 2, 0, 0, -s_sqrt / 2]
    mv =  np.sqrt((s_sqrt**2-20*s_sqrt))*raw_dot[0]
    # energy of final state photon. Must be > 10 GeV to avoid soft photon singularity and be seen by detectors.
    e_1 = (s_sqrt ** 2 - mv ** 2) / (2 * s_sqrt)
    e_2 = e_1 * raw_dot[1] + s_sqrt / 2 - e_1
    e_3 = s_sqrt-e_1-e_2
    # return [e_1,e_2,s_sqrt-e_1-e_2]
    phi_1 = 2 * np.pi * raw_dot[2]
    phi_12 = 2 * np.pi * raw_dot[3]
    cos_theta1 = -0.984 + 2 * raw_dot[4]*0.984
    #photon angle cut.
    sin_theta1 = np.sin(np.arccos(cos_theta1))
    cos_theta12 = ((s_sqrt - e_1 - e_2) ** 2 - e_1 ** 2 - e_2 ** 2) / (2 * e_1 * e_2)
    cos_theta2 = cos_theta1 * cos_theta12 + sin_theta1 * np.sqrt(1 - cos_theta12 ** 2) * np.cos(phi_12)
    sin_theta2 = np.sin(np.arccos(cos_theta2))
    phi_2 = phi_1 - np.arccos((cos_theta12 - cos_theta2 * cos_theta1) / (sin_theta1 * sin_theta2))
    k1 = e_1 * np.array([1, sin_theta1 * np.cos(phi_1), sin_theta1 * np.sin(phi_1), cos_theta1])
    k2 = e_2 * np.array([1, sin_theta1 * np.cos(phi_2), sin_theta2 * np.sin(phi_2), cos_theta2])
    k3 = np.array([s_sqrt, 0, 0, 0]) - k1 - k2
    # if np.random.randint(1,):
    #     k2,k3 = k3,k2

    weight = 1 / (8*(2*np.pi**5)*s_sqrt)
    return weight*scatteringamplitude.eenunugamma1(vectors=np.array([k1,k2,k2]))

# print(phase_space_23_dot_nt2_massless_cut(s_sqrt=500,raw_dot=np.random.rand(5,)))

integ = vegas.Integrator([[0, 1], [0, 1], [0, 1], [0, 1],[0, 1]])
result = integ(g, nitn=10, neval=1000)
print(result.summary())
print(result)
print(result.Q)
# print('result = %s Q = %.+2f' % (result, result.Q))