#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cm
size = 64 

n_sample = int(raw_input("t_mc:")) 
#X,Y = np.meshgrid(np.arange(0, 2 * np.pi, .2), np.arange(0, 2 * np.pi, .2))
X,Y = np.meshgrid(np.arange(size), np.arange(size))
list_sample = np.arange(0,3) 
list_sample = np.append(list_sample,np.arange(n_sample+1,n_sample+3))
print list_sample
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
    #fig, ax = plt.subplots(Q)
    name_image = "./output-"+str(n)+".pdf"
    qk = plt.quiverkey(Q, 0.9, 0.9, 2, r'$2 \frac{m}{s}$', labelpos='E',coordinates='figure')
    plt.savefig(name_image)
    plt.close()
    #mat = ax.matshow(np.copy(table),animated=True)
    #plt.show()
