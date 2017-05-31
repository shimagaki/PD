import numpy as np


if __name__ == '__main__':
    fin_name = input("Input file name:")
    fin = open(str(fin_name),"r") 
    data = fin.readlines()
    fin.close()
    foutname = "output.txt" 
    fout = open(foutname,"w")
    for line_data in data:
        line = [x for x in line_data.strip().split(' ')] 
        print(line[2],"  ",line[3],"  ",line[4],"\n")
        fout.write(str(line[2])+"  "+str(line[3])+"  "+str(line[4])+"  "+"\n")
    fout.close()

