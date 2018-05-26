import numpy as np
import math

def measurement():
    n = int(input('Please enter how many dots to cast:\n'))
    in_counter = 0
    total_counter = 0
    position = np.array([0,0])
    for i in range(n):
        step = 0.5*(2 * np.random.rand(2, ) - np.array([1, 1]))
        next_position = position + step
        if abs(next_position[0]) <= 1 and abs(next_position[1]) <= 1:
            position = next_position
            in_counter += 1
        else:
            pass
        total_counter += 1
    return 4*(in_counter/total_counter)

if __name__ == '__main__':
    print(measurement())