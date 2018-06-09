import sympy as sym
import numpy as np
import LorentzGenerators as lg
import Lorentz_Transformation as lt
import time
# np.random.seed(42)

# def rest_component_finder(E1,E2,E3,p1x,p1y,p1z,p2z):
def phase_space_dot_23_nt_massless(s_sqrt):
    while 1:
        raw_dot = np.random.rand(5,)
        E_gamma = ((s_sqrt/2)-10)*raw_dot[0]+10
        E_1 = raw_dot[3]*E_gamma+((s_sqrt/2)-E_gamma)
        E_2 = s_sqrt-E_1-E_gamma
        # 0.984 is cos10.
        theta_g = np.pi*raw_dot[1]
        phi_g = 2*np.pi*raw_dot[2]
        theta_1 = np.pi*raw_dot[4]
        cos_theta_1 = np.cos(theta_1)
        p3_g = lt.three_vector_gen_off_t_p(E_gamma,theta_g,phi_g)
        p_g = lt.Lorentz4vector(components=[E_gamma]+p3_g)
        pg = p_g.get_4_vector()
        p1z = E_1*np.cos(theta_1)
        p2z = -p3_g[-1]-p1z
        cos_theta_2 = np.arccos(p2z/E_2)
        p1_visible = E_1<10 or np.abs(cos_theta_1)>0.984
        p2_visible = E_2<10 or np.abs(cos_theta_2)>0.984
        if not (p1_visible and p2_visible):
            continue
        ####################
        # Energy generating algorithm validation. Positive result confirmed.
        # flag = True
        # if E_1+E_2<=E_gamma or abs(E_1)-E_2>=E_gamma or E_1+E_gamma<=E_2 or abs(E_1-E_gamma) >= E_2 or E_2 + E_gamma <= E_1 or abs(E_2 - E_gamma) >= E_1:
        #     flag = False
        # return flag
        #######################
        p1x = sym.Symbol('p1x')
        p1y = sym.Symbol('p1y')
        p2x = sym.Symbol('p2x')
        p2y = sym.Symbol('p2y')
        # p5theta = sym.Symbol('p5theta')
        # E_5 = sym.Symbol('E_5')
        #
        aa = sym.solve([p1x+p2x+pg[1],p1y+p2y+pg[2],p1x**2+p1y**2+p1z**2-E_1**2,p2x**2+p2y**2+p2z**2-E_2**2],[p1x,p1y,p2x,p2y])
        # print(E_gamma+E_1+E_2-s_sqrt)
        # print(aa)
        try:
            aa = np.float64(aa)
        except TypeError:
            continue
        p1 = lt.Lorentz4vector(components=np.array([E_1, aa[0][0], aa[0][1], p1z]), mass=0)
        p2 = lt.Lorentz4vector(components=np.array([E_2, aa[0][2], aa[0][3], p2z]), mass=0)
        # print(p1 + p2 + pg)
        # print(p1.get_on_shell_ness())
        # print(p2.get_on_shell_ness())
        return aa[0]
        # print(type(aa[0][0]))
        # print([E_1,aa[0][0],aa[0][1],p1z])

    # aa = sym.solve([p5norm*sym.cos(p5theta)+p4z+pgz],[p5theta])
    # aa = sym.solve([])

    # pg = lt.Lorentz4vector(components=[E_gamma,E_gamma*(np.sin())])
start = time.clock()
for i in range(100):
    phase_space_dot_23_nt_massless(500)
    # print('bing!')
end = time.clock()
print(end-start)
######
# Summary on massless case:
# 1. It works sometimes. But not always.
# 2. Sometimes the equation solver just cannot find a proper Real root and hence returns a complex root set.
# 3. The chance is 50/50.
# 4. When the collider energy is 500Gev, I believe the massless generator is good enough
#    since m_e is only ~500keV, one in a millionth of the total energy.
#####
# print(phase_space_dot_23_nt_massless(1))
# F = 0
# for i in range(10000):
#     if not phase_space_dot_23_nt_massless(1):
#         F += 1
# print(F)
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
# phase_space_dot_23_nt(s_sqrt=100,masses=[0,0,0])
#####
# It's a nice approach to generate photon energy using chain decay idea.
# However in this case all we need is E_gamma.