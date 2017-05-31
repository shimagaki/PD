import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
fig = plt.figure()
ims = []
for i in range(10):
        
    line_max =11 
    i=0
    table = np.zeros((line_max,line_max)) 
    fname = "xy-"+str(i)+".dat"
    for line in open(fname):
        itemList = line.split('  ')
        del itemList[-1]
        if(i==0):
            len_x = len(itemList)
        for j in range(len_x):
            table[i][j] = itemList[j]
        i+=1
    print "table =\n", table
    
    rand = np.random.randn(10)
    im = plt.plot(rand)        
    ims.append(im)             
ani = animation.ArtistAnimation(fig, ims, interval=100)
plt.show()
