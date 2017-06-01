#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cm
size = 30 

#n_sample = int(raw_input("t_mc:")) 
#X,Y = np.meshgrid(np.arange(0, 2 * np.pi, .2), np.arange(0, 2 * np.pi, .2))
X,Y = np.meshgrid(np.arange(size), np.arange(size))
t0 = 3.14*3
table = np.reshape(np.arange(0,t0,t0/(size*size))[:size*size],(size,size))
print "shape table = ", table
mat = np.copy(table) 
U = np.cos(mat)
V = np.sin(mat)
c =np.cos(mat) #0.5 * (np.copy(mat)+ np.ones((i,i))) 
print "c=", c

Q = plt.quiver(X, Y, U, V, c,angles='xy',scale_units='xy')
#fig, ax = plt.subplots(Q)
name_image = "./color-check3.pdf"
plt.savefig(name_image)
plt.close()
#mat = ax.matshow(np.copy(table),animated=True)
#plt.show()
