import numpy as np
import matplotlib.pyplot as plt

test1 = 'uv_init_height: '
test2 = 'u_init_data: '
test3 = 'v_init_data: '
with open('../InitialFields_Wind_Control.txt', 'r') as f:
    for line in f:
        if test1 in line:
            line = line.strip(test1)
            line = line.strip('\n')
            line = line.split(',')
            heights = np.array([float(z) for z in line])
        elif test2 in line:
            line = line.strip(test2)
            line = line.strip('\n')
            line = line.split(',')
            u = np.array([float(U) for U in line])
        elif test3 in line:
            line = line.strip(test3)
            line = line.strip('\n')
            line = line.split(',')
            v = np.array([float(V) for V in line])

wind_spd = np.sqrt(u**2 + v**2)

print 'uv_init_heights: '
print [round(z, 2) for z in heights]
print 'u_init_data: '
print u
print 'v_init_data: '
print v
print 'wind_spd: '
print wind_spd
print 'new_u: '
print [-round(w, 2) for w in wind_spd]
print 'new_v: '
print [0.00 for x in wind_spd]
