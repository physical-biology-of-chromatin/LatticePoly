import os
import vtkReader
import numpy as np
import os
import sys
import matplotlib.pyplot as plt


n = input("Please enter number of simulations:\n")
n=int(n)

def Rg(pos):
    Rg2=[]
    #take one frame
    for i in range(0,len(pos)):
        rg=[]
        #calculate center of mass
        rcm=sum(pos[i])/len(pos[i])
        #calculate squared difference for each monomer in the frame
        for j in range(0,len(pos[i])):
            rg.append((pos[i][j]-rcm).dot((pos[i][j]-rcm).T))
        #sum and and normalize
        Rg2.append(sum(rg)/len(pos[i]))
    return Rg2


rg2_collection=[]
for j in range(0,n):
	pos=[]
	reader=vtkReader.vtkReader('data/output/'+str(j),readLiq=0)
	for i in range(0,reader.N):
		pos.append(reader.polyPos)
		reader=next(reader)
	rg2_collection.append(Rg(pos))

np.save('data/gyrlattice/'+str(int(reader.boxDim[0]))+'_lattice',rg2_collection)

meanRg=sum(np.array(rg2_collection))/len(rg2_collection)
plt.plot(np.array(range(0,len(meanRg))), meanRg, 'r')
plt.show()
