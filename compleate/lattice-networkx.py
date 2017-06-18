# -*- coding: utf-8 -*- 
# creating a new network 
import networkx as nx 
import matplotlib.pyplot as plt
import numpy as np
L = 5
G = nx.Graph()
#N = np.arange(L*L)
#E=[]
for i in range(L):
    for j in range(L):
        i_plus = (i+1)%L 
        j_plus = (j+1)%L
        G.add_node( i*L+j, posxy = (i,j) )
        G.add_edge( i*L+j, i_plus*L+j )
        G.add_edge( i*L+j, i*L+j_plus )

pos = nx.get_node_attributes(G,'posxy') # posxy is a id of key.
nx.draw_networkx(G, pos, node_color='b')
plt.axis('off')
plt.show()
