import numpy as np
import matplotlib.pyplot as plt
import math
from force_density import construct_matrices, enforce_planarity,solve_and_plot,plot,get_bars

n_list = np.array([
        [1,2,-1,-1,-1,-1],
        [0,2,3,-1,-1,-1],
        [0,1,4,9,-1,-1],
        [1,4,5,7,-1,-1],
        [2,3,6,7,8,9],
        [3,6,7,-1,-1,-1],
        [4,5,7,8,10,-1],
        [3,4,5,6,-1,-1],
        [4,6,9,10,-1,-1],
        [2,4,8,10,11,-1],
        [6,8,9,11,-1,-1],
        [9,10,-1,-1,-1,-1]
        ],dtype=np.int32)

# fixed_nodes=[0,5,8,11]
# free_nodes=[1,2,3,4,6,7,9,10,12]

# force_densities=np.array([1,1,-2,2,-2,1,2,2,-2,1,2,-2,-2,-2,-2,-2,-2,1,-2,2,1,2,-2,1,2,1,2,-2])
# half_loaded_nodes = [0,5,11,8]
# unit_loaded_nodes = [1,2,6,10]
# double_loaded_nodes = [3,7,9,12]
# p = np.zeros((13,3))
# p[unit_loaded_nodes,2] = -1
# p[half_loaded_nodes,2] = -0.5
# p[double_loaded_nodes,2] = -1
# print(p)


coords = np.array([
    [0,0,9],
    [0,1,8],
    [-1,1,7],
    [0,2,6],
    [-1,2,0],
    [0,2.5,4.5],
    [-1,2.5,3.5],
    [-0.5,2.25,0],
    [-1.5,2.25,0],
    [-2,2,3],
    [-2,2.5,1.5],
    [-2.5,2.5,0]
])

coords[4,2] = enforce_planarity(coords[1],coords[2],coords[3],coords[4])
coords[8,2] = (coords[9,2]+coords[4,2]+coords[6,2]+coords[10,2])/4
coords[7,2] = (coords[9,2]+coords[4,2]+coords[6,2]+coords[10,2])/4

bars = get_bars(n_list)
print(bars)
plot(coords,bars,'quarter_des1')