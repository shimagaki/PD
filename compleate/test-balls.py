import numpy as np
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d

pi = np.pi
L = 30  
theta = 2 * pi / L  
l = 10.0
c = 1.0
a = (l/np.sqrt(L))
#fname="moc-triangle-N"+str(N)+"-d-"+str(d)+"-l"+str(l)+".dat"
fname="moc-spher-N"+str(L)+"-a-"+str(a)+".dat"
#for triangler#
#p1, p2 are fixed#
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

def gen_ball_random(theta):
    f=open(fname, "w")
    p1 =  np.array([np.cos(-theta/2.0),np.sin(-theta/2.0),0])
    p2 =  np.array([np.cos(theta/2.0),np.sin(theta/2.0),0])
    p3 =  np.array([0.5,0,np.sqrt(3)/2.0])
    # rotation axis vector 
    ax = np.zeros(3)
    t1,t2 = 2*pi*np.random.rand(), 2*pi*np.random.rand()
    ax[0]= np.cos(t1)*np.cos(t2) 
    ax[1]= np.cos(t1)*np.sin(t2) 
    ax[2]= np.sin(t1)
    R_ax = rotM(ax)
    for i in range(L):
        for j in range(L):
            px = a*np.cos(theta*i)*np.cos(theta*j)+c
            py = a*np.cos(theta*i)*np.sin(theta*j)+c
            pz = a*np.sin(theta*i)+c
            p = [px,py,pz]
            q = np.dot(R_ax,p)
            f.write(str(q[0])+" "+str(q[1])+" "+str(q[2])+"\n") 
    f.close()

if __name__ == '__main__': 
    gen_ball_random(theta)   
   
