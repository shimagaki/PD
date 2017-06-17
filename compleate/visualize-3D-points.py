import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d

fname0 = "moc-pyramid-N4-d-0.3-l5.0.dat"
data0 = np.loadtxt(fname0)
size_data = len(data0) 
# Create plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, axisbg="1.0")
ax = fig.gca(projection='3d')
ax.scatter3D(data0.T[0],data0.T[1],data0.T[2])
plt.title("N="+str(size_data))
plt.legend(loc=2)
plt.show()
