import MC_integration as mci
import matplotlib.pyplot as plt
import numpy as np

#Area of a circle
def MC_demo():
    count = 0
    dots_in = []
    dots_out = []
    for i in range(10000):
        dot_coordinate = np.array([-1,-1])+2*np.random.rand(2,)
        if dot_coordinate[0]**2 + dot_coordinate[1]**2 < 1:
            dots_in.append(dot_coordinate)
            count += 1
        else:
            dots_out.append(dot_coordinate)
    dots_in = np.array(dots_in)
    dots_out = np.array(dots_out)
    print(dots_in.shape)
    print(dots_in[0])
    plt.figure(num=3, figsize=(10, 10))
    theta = np.linspace(0,2*np.pi,10000)
    x = np.cos(theta)
    y = np.sin(theta)

    plt.scatter(dots_in[:,0],dots_in[:,1])
    plt.scatter(dots_out[:,0],dots_out[:,1])
    plt.plot(x,y,color = 'b')
    plt.xlim((-1,1))
    plt.xlabel('x')
    plt.ylim((-1,1))
    plt.ylabel('y')
    plt.show()
    print(count/10000)

MC_demo()

def circle():
    theta = np.linspace(0,2*np.pi,10000)
    x = np.cos(theta)
    y = np.sin(theta)
    plt.figure(num=3, figsize=(10, 10))
    plt.plot(x,y,color = 'b')
    plt.xlim((-1,1))
    plt.xlabel('x')
    plt.ylim((-1,1))
    plt.ylabel('y')
    plt.show()
# circle()