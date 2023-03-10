##
# radAveraging.py
# LatticePoly
##
# 
# Copyright © 2021 ENS Lyon. All rights reserved.
##

import os
import numpy as np

def fetchFiles(dirPath, fileName):
	fileList = []

	for root, dirs, files in os.walk(dirPath):
		if not dirs:
			if fileName in files:
				filePath = os.path.join(root, fileName)
				fileList.append(filePath)
	
	return fileList

fileList_1 = fetchFiles('/home/adminamith/3D/LatticePoly/LatticePoly/data/output/Living/Nucleation/', 'ContactFracLiq.res')
Nmeas    = 2032
rgyr_avg = np.zeros((Nmeas))
for i in range(0,len(fileList_1)):
    dt = np.loadtxt(fileList_1[i])
    #rgyr_avg[:,i] = np.sum(dt*dt, axis = 1)
    rgyr_avg= rgyr_avg + dt
rgyr_avg = rgyr_avg/len(fileList_1)   

np.savetxt('/home/adminamith/3D/LatticePoly/LatticePoly/data/output/Living/Nucleation/Liqcontact_Nucleation.res', rgyr_avg)
