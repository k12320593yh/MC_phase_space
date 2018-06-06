import numpy as np
import inspect as insp
import Lorentz_Transformation as lt
import scipy.linalg as sp


# np.random.seed(7)
def test_function_0(r, t, p):
    return (r ** 2) * np.sin(t)


test_bound_0 = np.array([[0, 1], [0, np.pi], [0, 2 * np.pi]])


def test_function_1(x, y, z):
    return 1


# simple 1D MC integrator
def MCI_1D(function, startpoint, endpoint, number_of_points=100000):
    summ = 0
    length = endpoint - startpoint
    dots = np.random.rand(number_of_points, )
    for i in range(dots.shape[0]):
        tmp = dots[i] * length + startpoint
        if i < 10:
            print(tmp)
        summ += dots[i]
    return (summ / number_of_points) * length


# print(MCI_1D(f,0,np.pi))

# string washer
def string_washer(f):
    st = insp.getargspec(f)
    tmp = []


# multi-dimensional MC integration
# Requires explicit integrand and bound to work
# Input shape:
# function:
#   multi-variable function as integrand.
# bound:
#   (dim,2) ndarray object. Each line of the arrary contains 2 elements,
#   first being the start point while second being the end point.
# number_of_points:
#   how many dots we cast to do the MC integraions. This affects the accuracy.
#   Casting less dots speeds up the programme at a cost of accuracy.

def MCI_MD(function, bound, number_of_points=100000):
    summ = 0
    args = insp.getargspec(function).args
    dim = len(args)
    lengths = bound[:, 1] - bound[:, 0]
    # get the variable list and dimensionality of the integration
    # and the integraion region length of each variable
    dots_0 = np.random.rand(number_of_points, dim)
    dots = np.zeros(number_of_points)
    for i in range(number_of_points):
        for j in range(dim):
            dots_0[i][j] = bound[j][0] + lengths[j] * dots_0[i][j]
        dots[i] = function(*dots_0[i])
        summ += dots[i]
    return (summ / number_of_points) * np.product(lengths)


# encodes all the information necessary to describe an 2->3 event.
# Generation of phase space dot
# Takes three arguments:
# momentum: p1 + p2 ,initial state particle momentum, a Lorentz 4 vector
# massess: 3, ndarry. Masses of final state particles.
# bound: phase space integration bound.
class Event(object):
    def __init__(self, vectors, weight, final_state_particles, masses, raw_dot):
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


# Naive Case
def two_two_phase_space_dot(s_sqrt, masses):
    raw_dot = np.random.rand(2, )
    momentum = lt.Lorentz4vector(components=[s_sqrt, 0, 0, 0])
    p_1 = lt.Lorentz4vector(components=[s_sqrt / 2, 0, 0, s_sqrt / 2])
    # momentum.lorentz_rot_toz()
    momentum.zrapidity()
    # momentum.z_boost_to_rest_frame()
    mass_0 = momentum.get_4_vector()[0]
    s = lt.four_vector_contraction(momentum, momentum)
    E_1 = (mass_0 ** 2 + masses[0] ** 2 - masses[1] ** 2) / (2 * mass_0)
    p_3_mag = np.sqrt(E_1 ** 2 - masses[0] ** 2)

    p_3 = lt.Lorentz4vector(components=[E_1, 0, 0, p_3_mag], mass=masses[0])
    p_4 = lt.Lorentz4vector(components=momentum.get_4_vector() - p_3.get_4_vector(), mass=masses[1])
    random_theta_rot_0 = lt.general_matrix_form(parameter=[0, np.arccos(raw_dot[0]), 0, 0, 0, 0])
    # print(random_theta_rot_0)
    random_phi_theta_rot_0 = np.matmul(lt.general_matrix_form(parameter=[0, 0, raw_dot[1] * (2 * np.pi), 0, 0, 0]),
                                       random_theta_rot_0)
    # print(random_phi_theta_rot_0)
    # print(p_3.get_4_vector())
    p_3.lorentz_transformation_off_matrix(random_phi_theta_rot_0)
    p_4.lorentz_transformation_off_matrix(random_phi_theta_rot_0)
    # p_3.sh()
    # print(np.matmul(random_phi_theta_rot_0,p_3.get_4_vector()))
    weight = (1 / (64 * (np.pi ** 2) * s)) * (p_3.get_mag() / p_1.get_mag())
    # print(weight)

    return Event(vectors=np.array([p_3.get_4_vector(), p_4.get_4_vector()]), raw_dot=raw_dot, weight=weight,
                 final_state_particles=2, masses=masses)


# two_two_phase_space_dot(s_sqrt=500,masses=[0,0])

# Cut mechanism. Returns a boolean value.
def cut(event):
    p1 = lt.Lorentz4vector(event.get_vector()[0], mass=event.get_mass()[0])
    p2 = lt.Lorentz4vector(event.get_vector()[1], mass=event.get_mass()[1])
    p3 = lt.Lorentz4vector(event.get_vector()[2], mass=event.get_mass()[2])
    z_vector = lt.Lorentz4vector(components=[1, 0, 0, 1], mass=0)
    # p3 will be taken as the momentum of the photon.
    # Reject photons that are too soft, which might lead to real emission singularity.
    if p3.e <= 10:
        # print("Photon to soft")
        return False
    # Reject photons too close to the beam.
    if np.abs(lt.cos_theta(p3, z_vector)) > 0.9848:
        # print("Photon collinear")
        return False

    p1visible = True
    p2visible = True
    if p1.e < 10 or np.abs(lt.cos_theta(p1, z_vector)) > 0.9848:
        p1visible = False
    if p2.e < 10 or np.abs(lt.cos_theta(p2, z_vector)) > 0.9848:
        p2visible = False

    if (not p1visible) and (not p2visible):
        return True

    return False


# This works, but is too slow. The dot generating algorithm needs to be modified to
# cast dot into a particular region we're interested in.
# And this is what next function does.

def focus_phase_space_dot(bounds):
    pass


# planned to be finished within tomorrow.


def two_three_phase_space_dot(s_sqrt, masses):
    # All final state particles must be on shell hence there're 3n-4 free variable to be integrated over. In our case,
    # it's 5.
    # print(s_sqrt)
    raw_dot = np.random.rand(5, )
    momentum = lt.Lorentz4vector(components=[s_sqrt, 0, 0, 0])
    # momentum.lorentz_rot_toz()
    momentum.zrapidity()
    # momentum.z_boost_to_rest_frame()
    # s in Mandelstam variables.
    mass_0 = momentum.get_4_vector()[0]
    s = lt.four_vector_contraction(momentum, momentum)
    if s_sqrt - np.sum(masses) < 0:
        print('Impossible process: energy conservation violation detected!')
        return None

    # determine how much energy is carried away by the first particle. Have to save enough for the next two!
    mediate_virtual_particle_mass = raw_dot[0] * (s_sqrt - np.sum(masses)) + masses[1] + masses[2]
    # print(mediate_virtual_particle_mass)
    # print(mediate_virtual_particle_mass)
    # After the mediator mass is determined, |p1| is determined as well:

    E_1 = (mass_0 ** 2 + masses[0] ** 2 - mediate_virtual_particle_mass ** 2) / (2 * mass_0)
    # print('E_1 = ',E_1)
    p_1_mag = np.sqrt(E_1 ** 2 - masses[0] ** 2)
    # Introducing angular randomness
    random_theta_rot_0 = lt.general_matrix_form(parameter=[0, raw_dot[1] * (2 * np.pi), 0, 0, 0, 0])
    random_phi_theta_rot_0 = np.matmul(lt.general_matrix_form(parameter=[0, 0, raw_dot[2] * (2 * np.pi), 0, 0, 0]),
                                       random_theta_rot_0)

    # theta_0_inverse = lt.general_matrix_form(parameter=[0,-np.arccos(raw_dot[1]),0,0,0,0])
    # phi_0_inverse = lt.general_matrix_form(parameter=[0,0,-raw_dot[2]*(2*np.pi),0,0,0])
    # Generate final state particle 1

    p_1 = lt.Lorentz4vector(name='final state particle 1', components=[E_1, 0, 0, p_1_mag], mass=masses[0])

    p_mediator = lt.Lorentz4vector(name='mediator particle',
                                   components=momentum.get_4_vector() - p_1.get_4_vector(),
                                   mass=mediate_virtual_particle_mass)
    # p_1_z = lt.Lorentz4vector(components=p_1.get_4_vector(),mass=p_1.get_mass())
    p_mediator_duplicate = lt.Lorentz4vector(name='mediator particle',
                                             components=momentum.get_4_vector() - p_1.get_4_vector(),
                                             mass=mediate_virtual_particle_mass)
    p_mediator.lorentz_transformation_off_matrix(random_phi_theta_rot_0)
    # print(p_1.components[0]+p_mediator.components[0])

    # The energy is not distributed equally among p_2 and p_3
    # Further investigation planned.

    # print(p_mediator.components[0])
    # print(p_mediator.get_4_vector())
    # p_1.sh(message='p_1 before rotation 0')
    p_1.lorentz_transformation_off_matrix(random_phi_theta_rot_0)
    # print("p_1 after rotation: " + str(p_1.get_4_vector()))
    # print(p_1.get_4_vector()+p_mediator.get_4_vector())
    # For now, the decay of original particle in to p1 and virtual particle conserves four momentum.

    # p_mediator_duplicate.lorentz_transformation_off_matrix(random_phi_theta_rot_0)
    # print(p_mediator_duplicate.get_4_vector())
    # So far so good.
    p_mediator_duplicate.zrapidity()
    p_mediator_duplicate.z_boost_to_rest_frame()
    # p_mediator_duplicate.sh()
    # print(p_mediator_duplicate.get_4_vector())
    # Reproduce last step and generate the next two. p_mediator is the new momentum l4v.
    decay_1_axial_rapidity = p_mediator_duplicate.get_zrapidity()
    decay_1_global_boost = lt.boost_matrix(rapidities=[0, 0, decay_1_axial_rapidity])
    # p_mediator_duplicate.lorentz_transformation_off_matrix(decay_1_global_boost)
    # print(p_mediator_duplicate.get_4_vector())
    # decay_1_global_boost works as expected.
    # For now, p_mediator is in its rest frame.

    ##############################

    # No Problem above this line.

    ##############################

    # And that is where the next pair of particles are generated.
    E_2 = (mediate_virtual_particle_mass ** 2 + masses[1] ** 2 - masses[2] ** 2) / (2 * mediate_virtual_particle_mass)
    p_2_mag = np.sqrt(E_2 ** 2 - masses[1] ** 2)

    # Generate random angle to rotate the second pair of particle.
    random_theta_rot_1 = lt.general_matrix_form(parameter=[0, raw_dot[3] * (2 * np.pi), 0, 0, 0, 0])
    random_phi_theta_rot_1 = np.matmul(lt.general_matrix_form(parameter=[0, 0, raw_dot[4] * (2 * np.pi), 0, 0, 0]),
                                       random_theta_rot_1)

    p_2 = lt.Lorentz4vector(name='final state particle 2', components=[E_2, 0, 0, p_2_mag], mass=masses[1])
    p_3 = lt.Lorentz4vector(name='final state particle 3',
                            components=p_mediator_duplicate.get_4_vector() - p_2.get_4_vector(), mass=masses[2])

    p_2.lorentz_transformation_off_matrix(random_phi_theta_rot_1)
    p_3.lorentz_transformation_off_matrix(random_phi_theta_rot_1)

    p_2.lorentz_transformation_off_matrix(decay_1_global_boost)
    p_3.lorentz_transformation_off_matrix(decay_1_global_boost)

    p_2.lorentz_transformation_off_matrix(random_phi_theta_rot_0)
    p_3.lorentz_transformation_off_matrix(random_phi_theta_rot_0)

    # Let's briefly recap this:
    # Generate a 1->2 decay event invovling a stationary particle with four momentum (sqrt(s),0,0,0), final state particle
    # and a virtual
    weight = 1 / ((4 * np.pi) ** 5 * (s_sqrt / 2) ** 2 * p_1.e * p_2.e * p_3.e)
    return Event(vectors=np.array([p_1.get_4_vector(), p_2.get_4_vector(), p_3.get_4_vector()]),
                 weight=weight, final_state_particles=3, masses=masses, raw_dot=raw_dot)



p1 = np.zeros(4, )
p2 = np.zeros(4, )
p3 = np.zeros(4, )
for i in range(10000):
    a = two_three_phase_space_dot(s_sqrt=5,masses=[0,0,0]).get_vector()
    p1 += a[0]
    p2 += a[1]
    p3 += a[2]

print(p1)
print(p2)
print(p3)
# Not enough randomness generated. Further review scheduled.
# for i in range(10000):
#     event = two_three_phase_space_dot(s_sqrt=500, masses=[0, 0, 0])
#     cut(event)
# This function extracts what is necessary from the final state particle momentums.
# M_square is Lorentz invariant, hence it only depends on the contraction of five ps.
# As the initial state momentum can be viewed as fixed,


# M_square must be a function of p1,p2 and p3.
# If not, then it will not be Lorentz invariable.

class Amplitude(object):
    def __init__(self, function, dimension):
        self.function = function
        self.dimension = dimension


def phase_space_integration(function, masses, number_of_dots=100, s_sqrt=500, dimension=2):
    summ = 0
    total_space = 0
    # count = 0
    momentum = lt.Lorentz4vector(components=[s_sqrt, 0, 0, 0], name='decaying mediator', mass=s_sqrt)
    M_square = Amplitude(function=function, dimension=dimension)
    abnornal_dots = []
    if dimension == 3:
        i = 0
        while i <= number_of_dots:
            # print("#######")
            phase_space_dot = two_three_phase_space_dot(s_sqrt=s_sqrt, masses=masses)
            # print(phase_space_dot.get_vector()[2])
            if cut(phase_space_dot):
                a = M_square.function(vectors=phase_space_dot.get_vector(),
                                      masses=masses)
                summ += a[0] * phase_space_dot.get_weight()
                total_space += phase_space_dot.get_weight()
                i += 1
                print(i)
            else:
                pass


    elif dimension == 2:
        for i in range(number_of_dots):
            phase_space_dot = two_two_phase_space_dot(s_sqrt=s_sqrt, masses=masses)
            a = M_square.function(vectors=phase_space_dot.get_vector(),
                                  masses=masses) * phase_space_dot.get_weight()
            summ += a
            total_space += phase_space_dot.get_weight()
            # print(a)
            # count+=1
    print(summ)
    return summ / number_of_dots


