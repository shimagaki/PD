#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cm
size =64 

X,Y = np.meshgrid(np.arange(size), np.arange(size))
fname = "all-spins-tMC0.dat"
fname_vortex = "vortex-tMC0.dat"

table = np.zeros((size,size)) 
i=0  
for line in open(fname,"r"):
    itemList = line.split(' ')
    #itemList = line.split('  ')
    del itemList[-1]
    if(i==0):
        len_x = len(itemList)
    for j in range(len_x):
        table[i][j] = itemList[j]
    i+=1
mat = np.copy(table) 
U = np.cos(mat)
V = np.sin(mat)
c = np.zeros(mat.shape) 

for line in open(fname_vortex):
    itemList = line[:-1].split(' ')
    i_x,i_y = int(itemList[0]), int(itemList[1])
    c[i_x][i_y] = 100 

Q = plt.quiver(X, Y, U, V, c,angles='xy',scale_units='xy')
name_image = "./test-output-t0-rapid-cooling.pdf"
plt.savefig(name_image)
plt.close()
