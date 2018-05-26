import numpy as np
import math

def measurement():
    n = int(input('Please enter how many dots to cast:\n'))
    in_counter = 0
    total_counter = 0
    for i in range(n):
        x = np.random.rand()
        y = np.random.rand()
        if x**2 + y**2 <= 1:
            in_counter += 1
        else:
            pass
        total_counter += 1
    return 4*(in_counter/total_counter)

if __name__ == '__main__':
    print(measurement())