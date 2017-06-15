import numpy as np

pi = np.pi
N = 10 
d = 0.3 # size of triangle. 
l = 10.0 # distance among triangles.
theta = pi / 3 # angular of isoscales triangle. 

#         z     p0 
#         ^     /\ 
#         |    /  \                 p1____p2    
#         |   /____\ <== theta      |      |
#         |   \    /<------------- -|      |
#         |    \  /                 |______|
#         |     \/p5              p3        p4
#        ----------------> xy
#        Theta is the angle of the middle and top(botom).

fname="moc-octahedral-N"+str(N)+"-d-"+str(d)+"-l"+str(l)+"-theta-"+str(theta)+".dat"

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
def gen_octahedral_random():
    p0 =  d*np.array([ 0, 0, np.tan(theta)])
    p1 =  d*np.array([ 0.5, 0.5, 0])
    p2 =  d*np.array([-0.5, 0.5, 0])
    p3 =  d*np.array([-0.5,-0.5, 0])
    p4 =  d*np.array([ 0.5,-0.5, 0])
    p5 =  d*np.array([0, 0, -np.tan(theta)])
    # rotation axis vector 
    f=open(fname, "w")
    for i in range(N):
        Pr = l * np.random.randn(3)
        ax = np.zeros(3)
        t1,t2 = 2*pi*np.random.rand(), 2*pi*np.random.rand()
        ax[0]= np.cos(t1)*np.cos(t2) 
        ax[1]= np.cos(t1)*np.sin(t2) 
        ax[2]= np.sin(t1)
        R_ax = rotM(ax)

        q0 =  np.dot(R_ax,p0) + Pr  #p0#
        q1 =  np.dot(R_ax,p1) + Pr  #p1#
        q2 =  np.dot(R_ax,p2) + Pr  #p2#
        q3 =  np.dot(R_ax,p3) + Pr  #p3#
        q4 =  np.dot(R_ax,p4) + Pr  #p4#
        q5 =  np.dot(R_ax,p5) + Pr  #p5#
        
        f.write(str(q0[0])+" "+str(q0[1])+" "+str(q0[2])+"\n") 
        f.write(str(q1[0])+" "+str(q1[1])+" "+str(q1[2])+"\n") 
        f.write(str(q2[0])+" "+str(q2[1])+" "+str(q2[2])+"\n") 
        f.write(str(q3[0])+" "+str(q3[1])+" "+str(q3[2])+"\n") 
        f.write(str(q4[0])+" "+str(q4[1])+" "+str(q4[2])+"\n") 
        f.write(str(q5[0])+" "+str(q5[1])+" "+str(q5[2])+"\n") 
    f.close()

if __name__ == '__main__': 
    gen_octahedral_random()   
   
