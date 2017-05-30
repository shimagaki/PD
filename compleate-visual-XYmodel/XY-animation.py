#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cm
size = 64 
fig, ax  = plt.subplots()
mats = []
n_sample = int(raw_input("t_mc:")) 
#X,Y = np.meshgrid(np.arange(0, 2 * np.pi, .2), np.arange(0, 2 * np.pi, .2))
X,Y = np.meshgrid(np.arange(size), np.arange(size))
list_sample = np.arange(0,50) 
np.append(list_sample,np.arange(n_sample,n_sample+50))
#for n in list_sample:
for n in list_sample:
    fname = "xy-tMC"+str(n)+".dat"
    table = np.zeros((size,size)) 
    i=0  
    for line in open(fname):
        itemList = line.split('  ')
        del itemList[-1]
        if(i==0):
            len_x = len(itemList)
        for j in range(len_x):
            table[i][j] = itemList[j]
        i+=1
    mat = np.copy(table) 
    U = np.cos(mat)
    V = np.sin(mat)
    c = np.cos(mat) 
    Q = plt.quiver(X, Y, U, V, c,angles='xy',scale_units='xy')
    #qk = plt.quiverkey(Q, 0.9, 0.9, 2, r'$2 \frac{m}{s}$', labelpos='E',coordinates='figure')
    #mat = ax.matshow(np.copy(table),animated=True)
    mats.append([Q]) 
ani = animation.ArtistAnimation(fig, mats, interval=800)
plt.show()
