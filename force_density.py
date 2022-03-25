import numpy as np
import matplotlib.pyplot as plt
import math
from mpl_toolkits import mplot3d

def construct_matrices(n_list,fixed_nodes,free_nodes,force_densities):

    nnodes = np.shape(n_list)[0]
    max_neigh = np.shape(n_list)[1]

    nbars = 0
    for node in n_list:
        for entry in node:
            if entry !=-1:
                nbars +=1
    nbars = int(0.5*nbars)
    nbars = 80

    Ca = np.zeros((nbars,nnodes))
    bars = []
    k=0
    for i in range(nnodes):
        for j in range(max_neigh):
            if n_list[i,j]==-1:
                pass
            else:
                if n_list[i,j]>i:
                    Ca[k,i]=1
                    Ca[k,n_list[i,j]]=-1
                    bars.append([i,n_list[i,j]])
                    k+=1
                else:
                    pass
    
    Cf = Ca[:,fixed_nodes]
    C = Ca[:,free_nodes]

    Q = np.eye(nbars)
    for i in range(nbars):
        Q[i,i] *= force_densities[i]

    D = np.matmul(np.matmul(np.transpose(C),Q),C)
    Df = np.matmul(np.matmul(np.transpose(C),Q),Cf)


    # print('Cf = {} \n \n C = {} \n \n D = {} \n \n Df = {}'.format(Cf,C,D,Df))
    return D, Df, bars

def solve_and_plot(D,Df,free_nodes,fixed_nodes,p,coords,bars,name):
    coords[free_nodes,2] = list(np.matmul(np.linalg.inv(D),p[free_nodes,2]-np.matmul(Df,coords[fixed_nodes,2])))
    # coords[free_nodes,1] = list(np.matmul(np.linalg.inv(D),p[free_nodes,1]-np.matmul(Df,coords[fixed_nodes,1])))
    # coords[free_nodes,0] = list(np.matmul(np.linalg.inv(D),p[free_nodes,0]-np.matmul(Df,coords[fixed_nodes,0])))

    return coords

def plot(coords,bars,name):

    for bar in bars:
        plt.plot(coords[bar,0],coords[bar,2])
    plt.savefig('{}_xz'.format(name))
    plt.close()

    for bar in bars:
        plt.plot(coords[bar,0],coords[bar,1])
    plt.savefig('{}_xy'.format(name))
    plt.close()

    for bar in bars:
        plt.plot(coords[bar,1],coords[bar,2])
    plt.savefig('{}_yz'.format(name))
    plt.close()

    ax = plt.axes(projection='3d')
    for bar in bars:
        ax.plot3D(coords[bar,0],coords[bar,1],coords[bar,2])
    plt.show()
    plt.savefig('{}_3D'.format(name))
    plt.close()

def get_bars(n_list):

    nnodes = np.shape(n_list)[0]
    max_neigh = np.shape(n_list)[1]

    nbars = 0
    for node in n_list:
        for entry in node:
            if entry !=-1:
                nbars +=1
    nbars = int(0.5*nbars)

    bars = []
    k=0
    for i in range(nnodes):
        for j in range(max_neigh):
            if n_list[i,j]==-1:
                pass
            else:
                if n_list[i,j]>i:
                    bars.append([i,n_list[i,j]])
                    k+=1
                else:
                    pass

    return bars


def enforce_planarity(coord1,coord2,coord3,coord4):

    matrix = np.array([
        [coord1[0],coord1[1],1],
        [coord2[0],coord2[1],1],
        [coord3[0],coord3[1],1]
    ])
    vector = np.array([coord1[2],coord2[2],coord3[2]])

    (A,B,C) = np.linalg.solve(matrix,vector)
    z4 = A*coord4[0] + B*coord4[1] + C

    return z4
