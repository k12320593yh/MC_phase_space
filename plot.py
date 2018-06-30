import matplotlib.pyplot as plt
import numpy as np

def y(r):
    return (0.4*r)/(1-0.8*r)

x = np.linspace(0,1,10000)
# yy = np.zeros_like(x)
# for i in range(x.shape[0]):
#     yy[i] = y(x)
yy = y(x)
plt.plot(x,yy)
plt.xlabel('r')
plt.ylabel('y')
plt.show()
