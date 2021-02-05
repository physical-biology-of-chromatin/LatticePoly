import os
import vtkReader
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import math
forkpos=0.0001*500*np.arange(1000)
vmax=0
averagematrices1=[]
for step in range(0,1):
	print(step)
	matrices_step=np.zeros((100,100))
	for i in range(0,50):
		data = np.loadtxt('/Volumes/MyPassport/Simulations_monday23/Set3/en1/'+str(i)+'/Matrixcis'+str(step)+'.txt')
		matrices_step+=data
	vmax=np.max(matrices_step)
	averagematrices1.append(matrices_step)
	averagematrices1.append(matrices_step/50)

os.makedirs('/Volumes/MyPassport/Yeast360/et3/en1/averagecis2/')


for i in range(len(averagematrices1)):
	fig, ax = plt.subplots()
	im = ax.imshow(np.log(averagematrices1[i]+0.1),vmin=math.log(0.1), vmax=math.log(vmax))
	fig.colorbar(im)
	rect = patches.Rectangle((50-int(forkpos[i]),50-int(forkpos[i])),0+2*int(forkpos[i]),0+2*int(forkpos[i]),linewidth=1,edgecolor='r',facecolor='none')
	ax.add_patch(rect)
	plt.suptitle('Simulation step = '+str(i),size=15)
	plt.savefig('/Volumes/MyPassport/Yeast360/et3/en1/averagecis2/%04d.png' %i)


