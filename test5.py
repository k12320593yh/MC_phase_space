import numpy as np
import Lorentz_Transformation as lt

print(lt.gamma_matrices_trace_generator(5))

# a = np.ones((4,4))
# b = lt.Lorentz_order_2_tensor(components=a)
# a = np.array([a,a,a,a])
# print(a.shape)
# In this way we can stack the tensor order up to any particular order we like.
#
# b = np.random.rand(4,4,4,4,4,)
# b = np.random.rand(4,4,4,4,4,)
# print(b[1,2,3,2,2])
# c = (1,2,3,2,3)
# print(type(c))
# print(b[c])
#it is possible to index a np array using a tuple. 