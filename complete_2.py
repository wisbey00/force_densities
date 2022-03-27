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

# Define the fixed nodes.
free_nodes = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,29,30,31,32,33,34,35,36]
fixed_nodes = [25,26,27,28]

# Define the loading on the nodes.
p = np.zeros((37,3))
p[groups[0],2] = -1
p[groups[1],2] = -1
p[groups[3],2] = -1/4
p[groups[4],2] = -1/2
p[groups[5],2] = -1/4
p[groups[6],2] = -1/2
p[groups[7],2] = -1/2
p[groups[2],2] = -3         # This is to model the weight of a lighting rig attached to nodes in group 1.

# plot(coords,get_bars(n_list),'complete_des')

for i in range(1):

    # Define the horizontal compression forces in each arch.
    # A value of 1 is the force equal to the gravity 
    # load applied to nodes in the central cross.
    Fo = 5
    Fi = 5

    # Define the force in the outer bars between nodes in groups 7. 
    # Force is factor 1 times force in central arch plus factor 2 
    # times force in outer arch.
    factor1 = 2.5
    factor2 = 1

    # Define the force densities as a function of the horizontal forces in the arches, 
    # dependent upon the distribution of the horizontal forces between bars 7,7 and 6,7.
    force_densities_mat = np.zeros((8,8))
    force_densities_mat[7,7] = factor1*Fi + factor2*Fo
    force_densities_mat[0,1] = -Fi
    force_densities_mat[1,2] = -Fo
    force_densities_mat[1,4] = -Fi
    force_densities_mat[2,5] = -Fo
    force_densities_mat[3,5] = 1.5*Fi+Fo-2*force_densities_mat[7,7]
    force_densities_mat[3,6] = 4*Fi+4*Fo-4*force_densities_mat[7,7]
    force_densities_mat[3,7] = Fi + 2*Fo
    force_densities_mat[4,7] = -Fi
    force_densities_mat[4,5] = 1.5*Fi+Fo-2*force_densities_mat[7,7]
    force_densities_mat[5,7] = -2*Fo
    force_densities_mat[6,7] = -4*Fi/3-4*Fo/3+4*force_densities_mat[7,7]/3

    # Define the bars in the system and check that they are unique.
    return_bars = get_bars(n_list)
    bars = []
    for i,bar in enumerate(return_bars):
        if bar in bars: pass
        elif [bar[1],bar[0]] in bars: pass
        else: bars.append(bar)
    # Take the force densities from the groupwise matrix to define them for each bar.
    force_densities = np.zeros(len(bars))
    for i,bar in enumerate(bars):
        group1 = reverse_groups[bar[0]] 
        group2 = reverse_groups[bar[1]]
        force_densities[i] = force_densities_mat[group1,group2]

    # Get the system matrices, solve for the z coordinates in the system, use this to 
    # get the forces in each bar then plot the resulting shape.
    D,Df,bars = construct_matrices(n_list,fixed_nodes,free_nodes,force_densities)
    new_coords = solve_and_plot(D,Df,free_nodes,fixed_nodes,p,coords,bars,'bigger_des2')
    print('Max height of the structure = {}'.format(np.round(4*np.max(new_coords[:,2])),3))
    print(get_forces(new_coords,bars,force_densities))
    plot(new_coords,bars,'complete_des','Equilibrium shape for the structure \n'
            'with arch forces of Fo = {} and Fi = {} \n and factor 1 = {}, factor 2 = {}'.format(Fo,Fi,factor1,factor2))

