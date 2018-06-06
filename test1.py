import sympy as sp
import numpy as np
import LorentzGenerators as lg
import Lorentz_Transformation as lt
# np.random.seed(42)

def phase_space_dot_23_nt(s_sqrt,masses):
    raw_dot = np.random.rand(5,)
    mv = raw_dot[0]*(s_sqrt-np.sum(masses))+np.sum(masses[1:])
    E_gamma = (s_sqrt**2-mv**2)/(2*s_sqrt)
    E_4 = raw_dot[3]*(s_sqrt-E_gamma-np.sum(masses[1:]))+np.sum(masses[2:])
    p4norm = np.sqrt(E_4**2-masses[1]**2)
    E_5 = s_sqrt-E_4-E_gamma
    theta_g = np.pi*raw_dot[1]
    phi_g = 2*np.pi*raw_dot[2]
    theta_4 = np.pi*raw_dot[4]
    pgx = np.sin(theta_g)*np.cos(phi_g)*E_gamma
    pgy = np.sin(theta_g)*np.sin(phi_g)*E_gamma
    pgz = np.cos(theta_g)*E_gamma
    p4z = np.cos(theta_4)*p4norm
# It's a nice approach to generate photon energy using chain decay idea.
# However in this case all we need is E_gamma.