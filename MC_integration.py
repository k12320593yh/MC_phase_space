import numpy as np

def integration():
    n = int(input('Please enter how many dots to cast:\n'))
    integral = 0
    for i in range(n):
        x = np.random.rand()
        integral += np.sin(x)*np.cos(x)
    return integral/n


if __name__ == '__main__':
    print(integration())