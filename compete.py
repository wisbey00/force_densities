import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import math
from force_density import construct_matrices, enforce_planarity,solve_and_plot,plot,get_bars

n_list = np.array([
        [1,10,-1,-1,-1,-1,-1],
        [0,2,3,-1,-1,-1,-1],
        [1,4,10,12,-1,-1,-1],
        [1,4,5,6,-1,-1,-1],
        [2,3,6,7,-1,-1,-1],
        [3,6,-1,-1,-1,-1,-1],
        [3,4,5,7,8,-1,-1],
        [4,6,8,9,11,12,15],
        [6,7,9,-1,-1,-1,-1],
        [7,8,15,-1,-1,-1,-1],
        [0,2,14,-1,-1,-1,-1],
        [7,12,13,14,15,-1,-1],
        [2,7,11,14,-1,-1,-1],
        [11,14,-1,-1,-1,-1,-1],
        [10,11,12,13,-1,-1,-1],
        [7,9,11,-1,-1,-1,-1]
        ],dtype=np.int32)

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

Fi = 20
Fo = 15
force_densities_mat = np.zeros((8,8))
force_densities_mat[0,1] = -Fi
force_densities_mat[1,2] = -Fo
force_densities_mat[1,4] = -Fi
force_densities_mat[2,5] = -Fo
force_densities_mat[3,5] = 0.5*Fi + Fo
force_densities_mat[3,6] = Fi + 2*Fo
force_densities_mat[3,7] = Fi + 2*Fo
force_densities_mat[4,7] = -Fi
force_densities_mat[4,5] = 0.5*Fi + Fo
force_densities_mat[5,7] = -2*Fo
force_densities_mat[6,7] = -2*Fi/3-4*Fo/3
force_densities_mat[7,7] = 0.5*Fi

return_bars = get_bars(n_list)
bars = []
for i,bar in enumerate(return_bars):
    print(bar)
    if bar in bars: pass
    elif [bar[1],bar[0]] in bars: pass
    else: bars.append(bar)
force_densities = np.zeros(len(bars))
for i,bar in enumerate(bars):
    group1 = reverse_groups[bar[0]] 
    group2 = reverse_groups[bar[1]]
    force_densities[i] = force_densities_mat[group1,group2]
# print(reverse_groups[bars], '\n \n', force_densities)
print(len(bars))

free_nodes = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,29,30,31,32,33,34,35,36]
fixed_nodes = [25,26,27,28]

p = np.zeros((37,3))
p[groups[0],2] = -1
p[groups[1],2] = -1

plot(coords,get_bars(n_list),'complete_des')

# print(np.shape(coords))

D,Df,bars = construct_matrices(n_list,fixed_nodes,free_nodes,force_densities)
new_coords = solve_and_plot(D,Df,free_nodes,fixed_nodes,p,coords,bars,'bigger_des2')
plot(new_coords,bars,'bigger_des2')

