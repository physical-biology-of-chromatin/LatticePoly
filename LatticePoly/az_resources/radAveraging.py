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

fileList_1 = fetchFiles('/home/aabdulla/3D/LatticePoly/LatticePoly/data/longtime_Jns20MB/', 'endtoend_10MBdom.res')
Nmeas    = 1000
rgyr_avg = np.zeros((Nmeas+1),256)
for i in range(0,len(fileList_1)):
    dt = np.loadtxt(fileList_1[i])
    rgyr_avg[:,i] = np.sum(dt[9000:10001]*dt[9000:10001], axis = 1)

rgyr_avg = rgyr_avg   

np.savetxt('/home/aabdulla/3D/LatticePoly/LatticePoly/data/longtime_Jns20MB/rendtoend_10MB.res', rgyr_avg)
