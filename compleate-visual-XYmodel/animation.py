#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

line_max =101 
fig, ax  = plt.subplots()
mats = []
n_sample = int(raw_input("t_mc:")) 
for j in range(n_sample):
    fname = "xy-beta40.000000-tMC"+str(j)+".dat"
    table = np.zeros((line_max,line_max)) 
    i=0  
    for line in open(fname):
        itemList = line.split('  ')
        del itemList[-1]
        if(i==0):
            len_x = len(itemList)
        for j in range(len_x):
            table[i][j] = itemList[j]
        i+=1
     
    mat = ax.matshow(np.copy(table),animated=True)
    mats.append([mat]) 
ani = animation.ArtistAnimation(fig, mats, interval=200)
plt.colorbar(mat) 
plt.show()
