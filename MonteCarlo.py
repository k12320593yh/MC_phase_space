import numpy as np
import inspect as insp
import Lorentz_Transformation as lt
import scipy.linalg as sp

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


#encodes all the information necessary to describe an 2->3 event.
#Generation of phase space dot
#Takes three arguments:
#momentum: p1 + p2 ,initial state particle momentum, a Lorentz 4 vector
#massess: 3, ndarry. Masses of final state particles.
#bound: phase space integration bound.
class Event(object):
    def __init__(self,vectors,weight):
        self.__vectors = vectors
        self.__weight = weight

    def get_vector(self):
        return self.__vectors

    def get_weight(self):
        return self.__weight

def two_three_phase_space_dot(momentum_0,momentum_1,masses):
#All final state particles must be on shell hence there're 3n-4 free variable to be integrated over. In our case,
#it's 5.
    raw_dot = np.random.rand(5,)
    print(raw_dot)
    momentum = lt.four_vector_add(momentum_0,momentum_1)
    momentum.lorentz_rot_toz()
    momentum.zrapidity()
    momentum.z_boost_to_rest_frame()
    #s in Mandelstam variables.
    mass_0 = momentum.get_4_vector()[0]
    s = lt.four_vector_contraction(momentum,momentum)
    # print(s)
    s_sqrt = np.sqrt(s)
    # print(s_sqrt)
    if s_sqrt - np.sum(masses) < 0:
        print('Impossible process: energy conservation violation detected!')
        return None

    #determine how much energy is carried away by the first particle. Have to save enough for the next two!
    mediate_virtual_particle_mass = raw_dot[0]*(s_sqrt-np.sum(masses))+masses[1]+masses[2]
    # print(mediate_virtual_particle_mass)
    #After the mediator mass is determined, |p1| is determined as well:

    E_1 = (mass_0**2+masses[0]**2-mediate_virtual_particle_mass**2)/(2*mass_0)
    # print('E_1 = ',E_1)
    p_1_mag = np.sqrt(E_1**2 - masses[0]**2)
    #Introducing angular randomness
    random_theta_rot_0 = lt.general_matrix_form(parameter=[0,np.arccos(raw_dot[1]),0,0,0,0])
    random_phi_theta_rot_0 = np.matmul(lt.general_matrix_form(parameter=[0,0,raw_dot[2]*(2*np.pi),0,0,0]),random_theta_rot_0)

    #Generate final state particle 1
    p_1 = lt.Lorentz4vector(name='final state particle 1',components=[E_1,0,0,p_1_mag],mass=masses[0])
    p_1.lorentz_transformation_off_matrix(random_phi_theta_rot_0)

    weight = p_1_mag/s_sqrt
    # Reproduce last step and generate the next two. p_mediator is the new momentum l4v.
    p_mediator = lt.Lorentz4vector(name='mediator particle',components=momentum.get_4_vector()-p_1.get_4_vector(),mass=mediate_virtual_particle_mass)

    #lt the mediator 4momemtum into rest fram at convenience of calulate next two particles.
    E_2 = (mediate_virtual_particle_mass**2+masses[1]**2-masses[2]**2)/(2*mediate_virtual_particle_mass)

    # E_2 = p_mediator.get_4_vector()[0]
    p_2_mag = np.sqrt(E_2**2 - masses[1]**2)

    random_theta_rot_1 = lt.general_matrix_form(parameter=[0,np.arccos(raw_dot[3]),0,0,0,0])
    random_phi_theta_rot_1 = np.matmul(lt.general_matrix_form(parameter=[0,0,raw_dot[4]*(2*np.pi),0,0,0]),random_theta_rot_1)

    p_2 = lt.Lorentz4vector(name='final state particle 2',components=[E_2,0,0,p_2_mag],mass=masses[1])
    p_2.lorentz_transformation_off_matrix(random_phi_theta_rot_1)

    p_3 = lt.Lorentz4vector(name='final state particle 3',components=p_mediator.get_4_vector()-p_2.get_4_vector(),mass=masses[2])

    weight *= mediate_virtual_particle_mass/p_2_mag

    return Event(vectors=np.array([momentum_0.get_4_vector(),momentum_1.get_4_vector(),p_1.get_4_vector(),p_2.get_4_vector(),p_3.get_4_vector()]),
                 weight=weight)



#This function extracts what is necessary from the final state particle momentums.
#M_square is Lorentz invariant, hence it only depends on the contraction of five ps.
#As the initial state momentum can be viewed as fixed,


#M_square must be a function of p1,p2 and p3.
#If not, then it will not be Lorentz invariable.
def MC_23_phase_space_integration(M_square,number_of_dots,s_sqrt):
    summ = 0
    total_space = 0
    momentum = lt.Lorentz4vector(components=[s_sqrt,0,0,0],name='decaying mediator',mass=s_sqrt)
    for i in range(number_of_dots):
        phase_space_dot = two_three_phase_space_dot(momentum,s_sqrt)
        summ += M_square(*phase_space_dot.get_vector())*phase_space_dot[-1]
        total_space += phase_space_dot[-1]

    return summ/total_space


#see if the processes is allowed by 4-momen
