import numpy as np
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d

pi = np.pi
N =  1000 
d = 0.3 # size of triangle. 
l = 5.0 # distance among triangles.
theta = pi / 3.0 # angular of isoscales triangle. 
#for triangler#
pi = np.pi
#p1, p2 are fixed#
fname="moc-pyramid-N"+str(N)+"-d-"+str(d)+"-l"+str(l)+".dat"
def rotM(p):
    px = p[0]
    py = p[1]
    pz = p[2]

    Rx = np.array([[1, 0, 0],
                   [0, np.cos(px), np.sin(px)],
                   [0, -np.sin(px), np.cos(px)]])
    Ry = np.array([[np.cos(py), 0, -np.sin(py)],
                   [0, 1, 0],
                   [np.sin(py), 0, np.cos(py)]])
    Rz = np.array([[np.cos(pz), np.sin(pz), 0],
                   [-np.sin(pz), np.cos(pz), 0],
                   [0, 0, 1]])
    R = Rz.dot(Ry).dot(Rx)
    return R

def gen_triangle_random(theta):
    p1 =  d*np.array([np.cos(-theta/2.0),np.sin(-theta/2.0),0])
    p2 =  d*np.array([np.cos(theta/2.0),np.sin(theta/2.0),0])
    p3 =  d*np.array([0,0,0])
    # rotation axis vector 
    f=open(fname, "w")
    for i in range(N):
        Pr = l*np.random.randn(3)
        ax = np.zeros(3)
        t1,t2 = 2*pi*np.random.rand(), 2*pi*np.random.rand()
        ax[0]= np.cos(t1)*np.cos(t2) 
        ax[1]= np.cos(t1)*np.sin(t2) 
        ax[2]= np.sin(t1)
        R_ax = rotM(ax)
        q1 = np.dot(R_ax,p1) + Pr 
        q2 = np.dot(R_ax,p2) + Pr
        q3 = np.dot(R_ax,p3) + Pr
        f.write(str(q1[0])+" "+str(q1[1])+" "+str(q1[2])+"\n") 
        f.write(str(q2[0])+" "+str(q2[1])+" "+str(q2[2])+"\n") 
        f.write(str(q3[0])+" "+str(q3[1])+" "+str(q3[2])+"\n") 
    f.close()


if __name__ == '__main__': 
    gen_triangle_random(theta)   
            
            
