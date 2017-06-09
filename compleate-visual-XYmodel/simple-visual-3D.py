from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import numpy as np
d = 16 
fig = plt.figure()
ax = fig.gca(projection='3d')
fname = "config-T0.100000-tMC0-random-field-3D.dat"
data = np.loadtxt(fname)
data=data.reshape((d,d,d))

x, y, z = np.meshgrid(np.arange(d),np.arange(d),np.arange(d))
u=np.cos(data) 
v=np.sin(data) 
w=np.ones(data.shape) 
#w=np.ones(data.shape) 
ax.quiver(x, y, z, u, v, w,length=0.1, norm=True)
plt.show()
