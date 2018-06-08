import sympy as sym
import numpy as np
import LorentzGenerators as lg
import Lorentz_Transformation as lt
# np.random.seed(42)

# def rest_component_finder(E1,E2,E3,p1x,p1y,p1z,p2z):

def phase_space_dot_23_nt(s_sqrt,masses):
    raw_dot = np.random.rand(5,)
    mv = raw_dot[0]*(s_sqrt-np.sum(masses))+np.sum(masses[1:])
    E_gamma = (s_sqrt**2-mv**2)/(2*s_sqrt)
    # E_4 = raw_dot[3]*(s_sqrt-E_gamma-np.sum(masses[1:]))+np.sum(masses[2:])
    # p4norm = np.sqrt(E_4**2-masses[1]**2)
    # E_5 = s_sqrt-E_4-E_gamma
    # p5norm = np.sqrt(E_5**2-masses[2]**2)
    theta_g = np.pi*raw_dot[1]
    phi_g = 2*np.pi*raw_dot[2]
    theta_4 = np.pi*raw_dot[4]
    pgx = np.sin(theta_g)*np.cos(phi_g)*E_gamma
    pgy = np.sin(theta_g)*np.sin(phi_g)*E_gamma
    pgz = np.cos(theta_g)*E_gamma
    # pg = lt.Lorentz4vector(components=[E_gamma,pgx,pgy,pgz],mass=0)
    # print(pg.get_on_shell_ness())
    # print(E_gamma)
    # print(E_4)
    # print(E_5)
    #All four components of pg is now generated.

    # p4z = np.cos(theta_4)*p4norm
    # p5z = -p4z-pgz

    # For p4 and p5, two of all four components is known.
    # All we have to do is to solve the 4 remaining components using four momentum conservation.

    p4x = sym.Symbol('p4x')
    p4y = sym.Symbol('p4y')
    p5x = sym.Symbol('p5x')
    p5y = sym.Symbol('p5y')
    p5theta = sym.Symbol('p5theta')
    # E_5 = sym.Symbol('E_5')
    #
    # aa = sym.solve([p4x+p5x+pgx,p4y+p5y+pgy,p4x**2+p4y**2+p4z**2+masses[1]**2-E_4**2,p5x**2+p5y**2+p5z**2+masses[2]**2-E_5**2],[p4x,p4y,p5x,p5y])
    # aa = sym.solve([p5norm*sym.cos(p5theta)+p4z+pgz],[p5theta])
    aa = sym.solve([])
    print(aa)
    # print(aa[0][1]+aa[0][3]+pgy)
    # return [E_gamma,E_4,E_5]
# E1 = 0
# E2 = 0
# E3 = 0
# for i in range(10000):
#     event = phase_space_dot_23_nt(s_sqrt=500,masses=[0,0,0])
#     E1 += event[0]
#     E2 += event[1]
#     E3 += event[2]
# print(E1,E2,E3)
phase_space_dot_23_nt(s_sqrt=100,masses=[0,0,0])
#####
# It's a nice approach to generate photon energy using chain decay idea.
# However in this case all we need is E_gamma.