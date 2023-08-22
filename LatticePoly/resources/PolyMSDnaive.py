##
##  PolyMSD.py
##  LatticePoly
##
##  Created by mtortora on 15/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import os
import sys
import psutil

import numpy as np

from utils import msdFFT
from vtkReader import vtkReader
import PolyMSD
import matplotlib
import matplotlib.pyplot as plt
import time





def msd_straight_forward(poHist):
    comulativemsds=[]
    for j in range(0,len(posHist[0])):
        r=[]
        for k in range(0,len(posHist)):
            r.append(posHist[k][j])
        r=np.array(r)
        shifts = np.arange(len(r))
        msds = np.zeros(shifts.size)    

        for i, shift in enumerate(shifts):
            diffs = r[:-shift if shift else None] - r[shift:]
            sqdist = np.square(diffs).sum(axis=1)
            msds[i] = sqdist.mean()
        comulativemsds.append(msds)
    return sum(comulativemsds)/len(posHist[0])





	
if __name__ == "__main__":
	if len(sys.argv) not in [3, 4]:
		print("\033[1;31mUsage is %s outputDir initFrame [idxTad]\033[0m" % sys.argv[0])
		sys.exit()

	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	
	msd = PolyMSD.PolyMSD(outputDir, initFrame=initFrame)
   
	fig, axs = plt.subplots(1,2, figsize=(15, 4), facecolor='w', edgecolor='k')
	axs = axs.ravel()
	fig.suptitle('MSD',size=20)
	
	posHist=msd.ReadHist()
	start_time = time.time()
	msds=msd_straight_forward(posHist)
	print("--- %s seconds ---" % (time.time() - start_time))
	

	axs[0].set_title('MSD(t) naive timecost='+str((time.time() - start_time)),size=15)
	axs[0].plot(np.arange(len(msds)), msds, 'ro')
	
	msd1 = PolyMSD.PolyMSD(outputDir, initFrame=initFrame)
	start_time = time.time()
	msd1.Compute()
	print("--- %s seconds ---" % (time.time() - start_time))
	msdHom = msd1.cumulDistHom / msd1.reader.nEuc


	axs[1].set_title('MSD(t) timecost='+str((time.time() - start_time)),size=15)
	axs[1].plot(np.arange(len(msdHom)), msdHom, 'ro')
	
	plt.savefig('MSD.png')    
	
	
	
