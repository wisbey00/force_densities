import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import math
from force_density import (construct_matrices,enforce_planarity,
                            solve_and_plot,plot,get_bars,get_forces)

# Define the starting coords of each node.
# Slightly unnecessary, only need to define the fixed nodes.
coords = np.array([
    [0,0,0],

    [1,0,0],
    [0,1,0],
    [-1,0,0],
    [0,-1,0],

    [1,1,0],
    [-1,1,0],
    [-1,-1,0],
    [1,-1,0],

    [2,2,0],
    [-2,2,0],
    [-2,-2,0],
    [2,-2,0],

    [2,0,0],
    [0,2,0],
    [-2,0,0],
    [0,-2,0],

    [2,1,0],
    [1,2,0],
    [-1,2,0],
    [-2,1,0],
    [-2,-1,0],
    [-1,-2,0],
    [1,-2,0],
    [2,-1,0],

    [2.5,2.5,0],
    [-2.5,2.5,0],
    [-2.5,-2.5,0],
    [2.5,-2.5,0],

    [2.5,1,0],
    [1,2.5,0],
    [-1,2.5,0],
    [-2.5,1,0],
    [-2.5,-1,0],
    [-1,-2.5,0],
    [1,-2.5,0],
    [2.5,-1,0]
])

# Due to the symmetry of the design, I grouped the nodes
# into 7 groups for which the nodes are identical.
groups = [
    [0],
    [1,2,3,4],
    [5,6,7,8],
    [9,10,11,12],
    [13,14,15,16],
    [17,18,19,20,21,22,23,24],
    [25,26,27,28],
    [29,30,31,32,33,34,35,36]
]
reverse_groups = np.zeros(37,dtype=np.int32)
for i in range(len(groups)):
    for entry in groups[i]:
        reverse_groups[entry] = i

# Set up the neighbour list array.
# The ith entry is a list of the neighbours of the ith node.
n_list = -1*np.ones((37,5))
n_list[0,:4] = groups[1]
for i,node in enumerate(groups[1]):  
    n_list[node,:4] = [0,groups[2][i],groups[2][i-1],groups[4][i]]
for i,node in enumerate(groups[2]):
    if i<3:     n_list[node,:4] = [groups[1][i],groups[1][i+1],groups[5][2*i],groups[5][2*i+1]]
    else:       n_list[node,:4] = [groups[1][-1],groups[1][i],groups[5][2*i],groups[5][2*i+1]]
for i,node in enumerate(groups[3]):
    n_list[node,:5] = [groups[5][2*i],groups[5][2*i+1],groups[6][i],groups[7][2*i],groups[7][2*i+1]]
for i,node in enumerate(groups[4]):
    if i==0:    n_list[node,:5] = [groups[1][i],groups[5][2*i],groups[5][2*i-1],groups[7][2*i],groups[7][2*i-1]]
    else:       n_list[node,:5] = [groups[1][i],groups[5][2*i-1],groups[5][2*i],groups[7][2*i-1],groups[7][2*i]]
for i,node in enumerate(groups[5]):
    if i<7:     n_list[node,:4] = [groups[2][i//2],groups[3][i//2],groups[4][(i+1)//2],groups[7][i]]
    else:       n_list[node,:4] = [groups[2][i//2],groups[3][i//2],groups[4][0],groups[7][i]]
for i,node in enumerate(groups[6]):
    n_list[node,:3] = [groups[3][i],groups[7][2*i],groups[7][2*i+1]]
for i,node in enumerate(groups[7]):
    if i<7 and i%2==0:     n_list[node,:4] = [groups[3][i//2],groups[4][(i+1)//2],groups[5][i],groups[7][i-1]]
    elif i<7 and i%2==1:   n_list[node,:4] = [groups[3][i//2],groups[4][(i+1)//2],groups[5][i],groups[7][i+1]]
    else:                  n_list[node,:4] = [groups[3][i//2],groups[4][0],groups[5][i],groups[7][0]]
n_list = np.array(n_list,dtype=np.int32)

coords[groups[0],2] = 0
coords[groups[1],2] = -1
coords[groups[2],2] = -2
coords[groups[3],2] = -8.88
coords[groups[4],2] = -4.42
coords[groups[5],2] = -5.42
coords[groups[6],2] = -12.63
coords[groups[7],2] = -12.63


plot(coords,get_bars(n_list),'airy','Airy stress function for equilibrium shape.')