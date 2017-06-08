#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
fname = "correlation.dat"
count=0
fig = plt.figure()
ims = []
for line in open(fname,"r"):
    itemList = line.split(' ')
    del itemList[-1]
    del itemList[0] # Corr(r=0)
    del itemList[1] # Corr(r=0)
    del itemList[2] # Corr(r=0)
    len_low = len(itemList)
    if(count==0):
        data = np.zeros(len_low)
        x = np.arange(len_low)
    
    for j in range(len_low):
        data[j] = itemList[j]
    im = plt.plot(x,data)  
    ims.append(im) 
    count +=1
ani = animation.ArtistAnimation(fig,ims,interval=500)
plt.show()
