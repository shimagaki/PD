import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
fig = plt.figure()
ims = []
for i in range(10):
        rand = np.random.randn(10)
        im = plt.plot(rand)        
        ims.append(im)             
ani = animation.ArtistAnimation(fig, ims, interval=100)
plt.show()
