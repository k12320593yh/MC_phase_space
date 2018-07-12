import numpy as np
import inspect as insp
import Lorentz_Transformation as lt
import scipy.linalg as sp
import sympy as sym
import time
# np.random.seed(7)

def test_function_0(r,t,p):
    return (r**2)*np.sin(t)

test_bound_0 = np.array([[0,1],[0,np.pi],[0,2*np.pi]])

def test_function_1(x,y,z):
    return 1

#simple 1D MC integrator
def MCI_1D(function,startpoint,endpoint,number_of_points=100000):
    summ = 0
    length = endpoint - startpoint
    dots = np.random.rand(number_of_points,)
    for i in range(dots.shape[0]):
        tmp = dots[i]*length+startpoint
        if i <10:
            print(tmp)
        summ += dots[i]
    return (summ/number_of_points)*length

# print(MCI_1D(f,0,np.pi))

#string washer
def string_washer(f):
    st = insp.getargspec(f)
    tmp = []

#multi-dimensional MC integration
#Requires explicit integrand and bound to work
#Input shape:
#function:
#   multi-variable function as integrand.
#bound:
#   (dim,2) ndarray object. Each line of the arrary contains 2 elements,
#   first being the start point while second being the end point.
#number_of_points:
#   how many dots we cast to do the MC integraions. This affects the accuracy.
#   Casting less dots speeds up the programme at a cost of accuracy.

def MCI_MD(function,bound,number_of_points=100000):
    summ = 0
    args = insp.getargspec(function).args
    dim = len(args)
    lengths = bound[:,1]-bound[:,0]
    #get the variable list and dimensionality of the integration
    #and the integraion region length of each variable
    dots_0 = np.random.rand(number_of_points,dim)
    dots = np.zeros(number_of_points)
    for i in range(number_of_points):
        for j in range(dim):
            dots_0[i][j] = bound[j][0] + lengths[j]*dots_0[i][j]
        dots[i] = function(*dots_0[i])
        summ += dots[i]
    return (summ/number_of_points)*np.product(lengths)


# Primitive MC integrator. Acuraccy ~ 0.1% with 100000 dots.

#encodes all the information necessary to describe an 2->3 event.
#Generation of phase space dot
#Takes three arguments:
#momentum: p1 + p2 ,initial state particle momentum, a Lorentz 4 vector
#massess: 3, ndarry. Masses of final state particles.
#bound: phase space integration bound.
class Event(object):
    def __init__(self,vectors,weight,final_state_particles,masses,raw_dot):
        self.__vectors = vectors
        self.__weight = weight
        self.number = final_state_particles
        self.__masses = masses
        self.raw_dot = raw_dot

    def get_vector(self):
        return self.__vectors

    def get_weight(self):
        return self.__weight

    def get_mass(self):
        return self.__masses


#Naive Case
def two_two_phase_space_dot(s_sqrt,masses):
    raw_dot = np.random.rand(2,)
    momentum = lt.Lorentz4vector(components=[s_sqrt,0,0,0])
    p_1 = lt.Lorentz4vector(components=[s_sqrt/2,0,0,s_sqrt/2])
    # momentum.lorentz_rot_toz()
    momentum.zrapidity()
    # momentum.z_boost_to_rest_frame()
    mass_0 = momentum.get_4_vector()[0]
    s = lt.four_vector_contraction(momentum,momentum)
    E_1 = (mass_0**2+masses[0]**2-masses[1]**2)/(2*mass_0)
    p_3_mag = np.sqrt(E_1 ** 2 - masses[0] ** 2)

    p_3 = lt.Lorentz4vector(components=[E_1,0,0,p_3_mag],mass=masses[0])
    p_4 = lt.Lorentz4vector(components=momentum.get_4_vector()-p_3.get_4_vector(),mass=masses[1])
    random_theta_rot_0 = lt.general_matrix_form(parameter=[0,np.arccos(raw_dot[0]),0,0,0,0])
    # print(random_theta_rot_0)
    random_phi_theta_rot_0 = np.matmul(lt.general_matrix_form(parameter=[0,0,raw_dot[1]*(2*np.pi),0,0,0]),random_theta_rot_0)
    # print(random_phi_theta_rot_0)
    # print(p_3.get_4_vector())
    p_3.lorentz_transformation_off_matrix(random_phi_theta_rot_0)
    p_4.lorentz_transformation_off_matrix(random_phi_theta_rot_0)
    # p_3.sh()
    # print(np.matmul(random_phi_theta_rot_0,p_3.get_4_vector()))
    weight= (1/(64*(np.pi**2)*s))*(p_3.get_mag()/p_1.get_mag())
    # print(weight)

    return Event(vectors=np.array([p_3.get_4_vector(), p_4.get_4_vector()]),raw_dot=raw_dot,weight=weight,final_state_particles=2,masses=masses)

# two_two_phase_space_dot(s_sqrt=500,masses=[0,0])

#Cut mechanism. Returns a boolean value.
def cut(event):
    return True
    p1 = lt.Lorentz4vector(event.get_vector()[0],mass=event.get_mass()[0])
    p2 = lt.Lorentz4vector(event.get_vector()[1],mass=event.get_mass()[1])
    p3 = lt.Lorentz4vector(event.get_vector()[2],mass=event.get_mass()[2])
    z_vector = lt.Lorentz4vector(components=[1,0,0,1],mass=0)
    # p3 will be taken as the momentum of the photon.
    # Reject photons that are too soft, which might lead to real emission singularity.
    if p3.get_4_vector()[0] <= 10:
        # print("Photon to soft")
        return False
    # Reject photons too close to the beam.
    if np.abs(lt.cos_theta(p3,z_vector))>0.9848:
        # print("Photon collinear")
        return False

    if (p1.get_4_vector()[0]<10 or np.abs(lt.cos_theta(p1,z_vector))>0.9848)&(p2.get_4_vector()[0]<10 or np.abs(lt.cos_theta(p2,z_vector))>0.9848):
        return True
    else:
        return False


# This works, but is too slow. The dot generating algorithm needs to be modified to
# cast dot into a particular region we're interested in.
# And this is what next function does.

def focus_phase_space_dot(bounds):
    pass
# planned to be finished within tomorrow.


def two_three_phase_space_dot(s_sqrt,masses):
#All final state particles must be on shell hence there're 3n-4 free variable to be integrated over. In our case,
#it's 5.
    # print(s_sqrt)
    raw_dot = np.random.rand(5,)
    momentum = lt.Lorentz4vector(components=[s_sqrt,0,0,0])
    momentum.zrapidity()
    #s in Mandelstam variables.
    mass_0 = momentum.get_4_vector()[0]
    s = lt.four_vector_contraction(momentum,momentum)
    if s_sqrt - np.sum(masses) < 0:
        print('Impossible process: energy conservation violation detected!')
        return None

    #determine how much energy is carried away by the first particle. Have to save enough for the next two!
    mediate_virtual_particle_mass = raw_dot[0]*(s_sqrt-np.sum(masses))+masses[1]+masses[2]
    # print(mediate_virtual_particle_mass)
    # print(mediate_virtual_particle_mass)
    #After the mediator mass is determined, |p1| is determined as well:

    E_1 = (mass_0**2+masses[0]**2-mediate_virtual_particle_mass**2)/(2*mass_0)
    p_1_mag = np.sqrt(E_1**2 - masses[0]**2)
    #Introducing angular randomness
    # random_theta_rot_0 = lt.general_matrix_form(parameter=[0,raw_dot[1]*(np.pi),0,0,0,0])
    # random_phi_theta_rot_0 = np.matmul(lt.general_matrix_form(parameter=[0,0,raw_dot[2]*(2*np.pi),0,0,0]),random_theta_rot_0)
    theta_0 = raw_dot[1]*np.pi
    phi_0 = raw_dot[2]*2*np.pi
    theta_1 = raw_dot[3]*np.pi
    phi_1 = raw_dot[4]*2*np.pi

    random_theta_rot_0 = lt.general_matrix_form(parameter=[0,raw_dot[1]*(np.pi),0,0,0,0])
    random_phi_theta_rot_0 = np.matmul(lt.general_matrix_form(parameter=[0,0,raw_dot[2]*(2*np.pi),0,0,0]),random_theta_rot_0)

    p_1_0 = lt.Lorentz4vector(components=[E_1,0,0,p_1_mag],mass=masses[0])

    #Generate final state particle 1
    p_1 = lt.Lorentz4vector(name='final state particle 1',components=[E_1,-p_1_mag*np.sin(theta_0)*np.cos(phi_0)
        ,p_1_mag*np.sin(theta_0)*np.sin(phi_0),p_1_mag*np.cos(theta_0)],mass=masses[0])
    p_mediator = lt.Lorentz4vector(name='mediator particle',
                                         components=momentum.get_4_vector() - p_1.get_4_vector(),
                                         mass=mediate_virtual_particle_mass)
    p_mediator_duplicate = lt.Lorentz4vector(name='mediator duplicate',
                                             components=[s_sqrt-E_1,
                                                         0,0,-p_1_mag])
    # p_mediator_duplicate =    lt.Lorentz4vector(name='mediator particle',
    #                                          components=momentum.get_4_vector()-p_1.get_4_vector(),
    #                                          mass=mediate_virtual_particle_mass)



    # p_mediator.lorentz_transformation_off_matrix(random_phi_theta_rot_0)
    # print(p_1.components[0]+p_mediator.components[0])

    #The energy is not distributed equally among p_2 and p_3
    #Further investigation planned.

    # print(p_mediator.components[0])
    # print(p_mediator.get_4_vector())
    # p_1.sh(message='test8p_1 before rotation 0')
    p_1_0.lorentz_transformation_off_matrix(random_phi_theta_rot_0)
    # print("p_1 after rotation: " + str(p_1.get_4_vector()))

    # print(p_1.get_4_vector()+p_mediator.get_4_vector())
    #For now, the decay of original particle in to p1 and virtual particle conserves four momentum.

    # p_mediator_duplicate.lorentz_transformation_off_matrix(random_phi_theta_rot_0)
    # print(p_mediator_duplicate.get_4_vector())
    # So far so good.
    p_mediator_duplicate.zrapidity()
    p_mediator_duplicate.z_boost_to_rest_frame()
    # p_mediator_duplicate.sh()
    # print(p_mediator_duplicate.get_4_vector())
    # Reproduce last step and generate the next two. p_mediator is the new momentum l4v.
    decay_1_axial_rapidity = p_mediator_duplicate.get_zrapidity()
    decay_1_global_boost = lt.boost_matrix(rapidities=[0,0,decay_1_axial_rapidity])
    # p_mediator_duplicate.lorentz_transformation_off_matrix(decay_1_global_boost)
    # print(p_mediator_duplicate.get_4_vector())
    # decay_1_global_boost works as expected.
    # For now, p_mediator is in its rest frame.

    ##############################

    #No Problem above this line.

    ##############################

    # And that is where the next pair of particles are generated.
    E_2 = (mediate_virtual_particle_mass**2+masses[1]**2-masses[2]**2)/(2*mediate_virtual_particle_mass)
    p_2_mag = np.sqrt(E_2**2 - masses[1]**2)

    #Generate random angle to rotate the second pair of particle.
    random_theta_rot_1 = lt.general_matrix_form(parameter=[0,raw_dot[3]*(np.pi),0,0,0,0])
    random_phi_theta_rot_1 = np.matmul(lt.general_matrix_form(parameter=[0,0,raw_dot[4]*(2*np.pi),0,0,0]),random_theta_rot_1)
    # random_phi_theta_rot_1 = lt.general_matrix_form(parameter=[0,raw_dot[3]*(2*np.pi),raw_dot[4]*(2*np.pi),0,0,0])
    # random_phi_theta_rot_1 = np.diag((1,np.sin(np.pi*raw_dot[3])*np.cos(2*np.pi*raw_dot[4])
    #                                   ,np.sin(np.pi*raw_dot[3])*np.sin(2*np.pi*raw_dot[4])
    #                                   ,np.cos(np.pi*raw_dot[3])))
    # After trying various rotation algorithm, the original one seems correct. Perfect!

    # p_2 = lt.Lorentz4vector(name='final state particle 2',components=[E_2,0,0,p_2_mag],mass=masses[1])
    # p_3 = lt.Lorentz4vector(components=p_mediator_duplicate.get_4_vector()-p_2.get_4_vector(),mass=masses[2])
    p_2 = lt.Lorentz4vector(name='final state particle 2',components=[E_2,-p_2_mag*np.sin(theta_1)*np.cos(phi_1),
                                                                      p_2_mag*np.sin(theta_1)*np.sin(phi_1),p_2_mag*np.cos(theta_1)],mass=masses[1])
    p_3 = lt.Lorentz4vector(name='final state particle 3',components=p_mediator_duplicate.get_4_vector() - p_2.get_4_vector(), mass=masses[2])

    # p_2.lorentz_transformation_off_matrix(random_phi_theta_rot_1)
    # p_3.lorentz_transformation_off_matrix(random_phi_theta_rot_1)


    p_2.lorentz_transformation_off_matrix(decay_1_global_boost)
    p_3.lorentz_transformation_off_matrix(decay_1_global_boost)


    p_2.lorentz_transformation_off_matrix(random_phi_theta_rot_0)
    p_3.lorentz_transformation_off_matrix(random_phi_theta_rot_0)


    # print(p_1+p_2+p_3)
    # p4x = sym.Symbol('p4x')
    # p4y = sym.Symbol('p4y')
    # p5x = sym.Symbol('p5x')
    # p5y = sym.Symbol('p5y')
    #
    # pgx = p_1.get_4_vector()[1]
    # pgy = p_1.get_4_vector()[2]
    # p4z = p_2.get_4_vector()[3]
    # E_4 = p_2.get_4_vector()[0]
    # E_5 = p_3.get_4_vector()[0]
    # p5z = p_3.get_4_vector()[3]
    # p1 = p_1.get_4_vector()
    # p2 = p_2.get_4_vector()
    # p3 = p_3.get_4_vector()
    # # print(p_2.get_on_shell_ness())
    # # print(p_3.get_on_shell_ness())
    # # print(p1[1]+p2[1]+p3[1])
    # # print(p1[2]+p2[2]+p3[2])
    # # print(p2[1]**2+p2[2]**2+p2[3]**2+masses[1]**2-p2[0]**2)
    # # print(p3[1]**2+p3[2]**2+p3[3]**2+masses[2]**2-p3[0]**2)
    #
    # aa = sym.solve([p4x + p5x + pgx, p4y + p5y + pgy, p4x ** 2 + p4y ** 2 + p4z ** 2 + masses[1] ** 2 - E_4 ** 2,
    #                 p5x ** 2 + p5y ** 2 + p5z ** 2 + masses[2] ** 2 - E_5 ** 2], [p4x, p4y, p5x, p5y])
    # p_2p = lt.Lorentz4vector(components=[E_4,aa[0][0],aa[0][1],p4z],mass=masses[1])
    # p_3p = lt.Lorentz4vector(components=[E_5,aa[0][2],aa[0][3],p5z],mass=masses[2])
    # print(p_1+p_2+p_3)
    # print(p_1+p_2p+p_3p)
    # print(p_2p.get_on_shell_ness())
    # print(p_3p.get_on_shell_ness())
    # print(aa[0])
    # print(p2[1],p2[2],p3[1],p3[2])
    #
    # print(p2[1]+p3[1]+pgx)
    # print(aa[0][0]+aa[0][2]+pgx)
    # print(p2[2]+p3[2]+pgy)
    # print(aa[0][1]+aa[0][3]+pgy)
    # print(aa[0][0]**2+aa[0][1]**2+p4z**2-E_4**2+masses[1]**2)
    # print(p2[1]**2+p2[2]**2+p2[3]**2-E_4**2+masses[1]**2)
    # print(aa[0][2]**2+aa[0][3]**2+p5z**2-E_5**2+masses[2]**2)
    # print(p3[1]**2+p3[2]**2+p3[3]**2-E_5**2+masses[2]**2)
    # return (aa[0][0]-p2[1])
    # print(aa[0][1] + aa[0][3] + pgy)
    # print(p_1+p_2+p_3)
    # print(p_2+p_3)
    #Let's briefly recap this:
    #Generate a 1->2 decay event invovling a stationary particle with four momentum (sqrt(s),0,0,0), final state particle
    #and a virtual
    # After the Energy bug was fixed, the equation-solving based algorithm seems to work in some circumstaces.
    # However, not in every case.
    weight = 1/((4*np.pi)**5*(s_sqrt/2)**2*p_1.get_4_vector()[0]*p_2.get_4_vector()[0]*p_3.get_4_vector()[0])
    return Event(vectors=np.array([p_1.get_4_vector(),p_2.get_4_vector(),p_3.get_4_vector()]),
                 weight=weight,final_state_particles=3,masses=masses,raw_dot=raw_dot)

# # #
# p1 = np.zeros(4, )
# p2 = np.zeros(4, )
# p3 = np.zeros(4, )
# for i in range(1000):
#     a = two_three_phase_space_dot(s_sqrt=500,masses=[0,0,0]).get_vector()
#     p1 += a[0]
#     p2 += a[1]
#     p3 += a[2]
#     # print(a[0]+a[1]+a[2])
# print(p1)
# print(p2)
# print(p3)
#Not enough randomness generated. Further review scheduled.
# T = 0
# F = 0
# for i in range(100):
#     result = two_three_phase_space_dot(s_sqrt=500, masses=[0, 5, 5])
#     if result<1e-6:
#         T+=1
#     else:
#         F+=1
# print(T)
# print(F)
# print(event.get_vector())
#     cut(event)
#This function extracts what is necessary from the final state particle momentums.
#M_square is Lorentz invariant, hence it only depends on the contraction of five ps.
#As the initial state momentum can be viewed as fixed,


#M_square must be a function of p1,p2 and p3.
#If not, then it will not be Lorentz invariable.

class Amplitude(object):
    def __init__(self,function,dimension):
        self.function = function
        self.dimension = dimension

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
    #This solves a full set of 4-variable quadratic euqations and drags the entire program significantly.
    print(aa)


# An attempt to generate dots within
# def phase_space_dot_23_nt_cut_massless(s_sqrt):
#     raw_dot = np.random.rand(5,)
#     E_gamma = raw_dot[0]*((s_sqrt/2)-10)+10
#     cos_theta_gamma = -0.984+2*raw_dot[1]*0.984
#     phi = 2*np.pi*raw_dot[2]
    # Generate photons hard enough.
    # Fatal failure. After the energies are determined, the angles are no longer totally free, but heavily restricted.


def phase_space_23_dot_nt2_massless_cut(s_sqrt):
    while 1:
        raw_dot = np.random.rand(5, )
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
        if np.abs(cos_theta2)<0.984:
            continue
        sin_theta2 = np.sin(np.arccos(cos_theta2))
        phi_2 = phi_1 - np.arccos((cos_theta12 - cos_theta2 * cos_theta1) / (sin_theta1 * sin_theta2))
        if e_3 > 10 and np.abs((e_1*cos_theta1+e_2*cos_theta2)/e_3)<0.984:
            continue
        k1 = e_1 * np.array([1, sin_theta1 * np.cos(phi_1), sin_theta1 * np.sin(phi_1), cos_theta1])
        k2 = e_2 * np.array([1, sin_theta1 * np.cos(phi_2), sin_theta2 * np.sin(phi_2), cos_theta2])
        k3 = np.array([s_sqrt, 0, 0, 0]) - k1 - k2
        # if np.random.randint(1,):
        #     k2,k3 = k3,k2

        weight = 1 / (8*(2*np.pi**5)*s_sqrt)
        return Event(vectors=np.array([k1,k2,k3]),
                     weight=weight, final_state_particles=3, masses=[0,0,0], raw_dot=raw_dot)

# start = time.time()
# a = np.zeros(4,)
# b = np.zeros(4,)
# d = np.zeros(4,)
# for i in range(100):
#     c = phase_space_23_dot_nt2_massless_cut(500)
#     a += c[0]
#     b += c[1]
#     d += c[2]
# end = time.time()
# print(start-end)
# #roughly 200 events/s with cut. Far better than before!
# print(a/100)
# print(b/100)
# print(d/100)

def phase_space_integration(function,masses,number_of_dots=100,s_sqrt=500,dimension=2,mode=0):
    summ = 0
    total_space = 0
    M_square = Amplitude(function=function,dimension=dimension)

    if mode == 1:
        #load already generated phase space dots into python instead of generating new.
        vectors = np.load('1000000_eenunua.npy')[1:]
        for i in range(number_of_dots):
            j = np.random.randint(0,10000,1)[0]
            phase_space_dot = Event(vectors=vectors[3*j:3*j+3],
                                    weight=1/((4*np.pi)**5*(s_sqrt/2)**2*vectors[3*j][0]*vectors[3*j+1][0]*
                                              vectors[3*j+2][0]),
                                    masses=[0,0,0],raw_dot=[0,0,0,0,0],final_state_particles=3)
            a = M_square.function(vectors=phase_space_dot.get_vector(),
                                  masses=masses)
            summ += a * phase_space_dot.get_weight()
            total_space += phase_space_dot.get_weight()
        return summ/number_of_dots
            # i += 1
            # print(i)

    # count = 0
    momentum = lt.Lorentz4vector(components=[s_sqrt,0,0,0],name='decaying mediator',mass=s_sqrt)
    abnormal_dots = []
    if dimension == 3:
        i = 0
        while i <= number_of_dots:
            # print("#######")
            phase_space_dot = phase_space_23_dot_nt2_massless_cut(s_sqrt=s_sqrt)
            # print(phase_space_dot.get_vector()[2])
            if cut(phase_space_dot):
                a = M_square.function(vectors=phase_space_dot.get_vector(),
                                      masses=masses)
                summ += a*phase_space_dot.get_weight()
                total_space += phase_space_dot.get_weight()
                i += 1
                print(i)
            else:
                pass


    elif dimension == 2:
        for i in range(number_of_dots):
            phase_space_dot = two_two_phase_space_dot(s_sqrt=s_sqrt,masses=masses)
            a = M_square.function(vectors=phase_space_dot.get_vector(),
                                      masses=masses) * phase_space_dot.get_weight()
            summ += a
            total_space += phase_space_dot.get_weight()
            # print(a)
            # count+=1
    print(summ)
    return summ/number_of_dots
#
# x = np.load("1000000_eenunua.npy")[1:]
# print(x.shape)
# for i in range(100):
#     y = lt.Lorentz4vector(components=x[i],mass=0)
#     print(y.get_on_shell_ness())
# def unit_func():
#     return 1
# start = time.clock()
# cnt = 0
# array = np.zeros(4,)
# while cnt < 100000:
#     event = two_three_phase_space_dot(s_sqrt=500,masses=[0,0,0])
#     if cut(event):
#         cnt+=1
#         array = np.vstack((array,event.get_vector()))
#         # print(cnt)
# np.save("100000_eenunua",array)
# end = time.clock()
#
# print(end-start)

#see if the processes is allowed by 4-momen
