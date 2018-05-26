import numpy as np

a = np.ones((4,4))
a = np.array([a,a,a,a])
print(a.shape)
#In this way we can stack the tensor order up to any particular order we like.

b = np.random.rand(4,4,4,4,4,)
print(b[1,2,3,2,2])
c = (1,2,3,2,3)
print(type(c))
print(b[c])
#it is possible to index a np array using a tuple. 