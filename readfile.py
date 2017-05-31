#!/usr/bin/env python
import numpy as np
i=0
line_max =11 
table = np.zeros((line_max,line_max)) 
for line in open('xy-3.dat'):
    itemList = line.split('  ')
    del itemList[-1]
    if(i==0):
        len_x = len(itemList)
    for j in range(len_x):
        table[i][j] = itemList[j]
    i+=1
print table

