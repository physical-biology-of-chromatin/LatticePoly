##
# radavg.py
# LatticePoly
##
# 
# Copyright © 2021 ENS Lyon. All rights reserved.
##

import os
import sys
import numpy as np



dpath = '/home/aabdulla/3D/Chap3/LatticePoly/LatticePoly/data/output/1MB_domain/L1MB/'
Nchain = 1078
Nmeas = 1000

def fetchFiles(dirPath, fileName):
	fileList = []

	for root, dirs, files in os.walk(dirPath):
		if not dirs:
			if fileName in files:
				filePath = os.path.join(root, fileName)
				fileList.append(filePath)
	
	return fileList

fileList_1 = fetchFiles(dpath, 'endtoend_L1MB.res')
fileList_2 = fetchFiles(dpath, 'msdEEvect_L1MB.res')
fileList_3 = fetchFiles(dpath, 'monomerdistmat_L1MB.res')
fileList_4 = fetchFiles(dpath, 'contactProb_L1MB.res')
fileList_5 = fetchFiles(dpath, 'rgyration_L1MB.res')

endtoendFile   = os.path.join(dpath, "endtoend_Avg.res")
msdFile        = os.path.join(dpath, "msdEE_Avg.res") 
distmatFile    = os.path.join(dpath, "distmat_Avg.res") 
cProbFile      = os.path.join(dpath, "contactProb_Avg.res")
rgyrFile       = os.path.join(dpath, "rgyr_Avg.res")  

#fileList_5 = fetchFiles('/home/aabdulla/3D/LatticePoly/LatticePoly/data/longtime_Jns20MB/','g2polyMSD_600kb_2.res')

##
rgyr_avg = np.zeros((Nmeas+1))
j = 0
for i in range(0,len(fileList_1)):
   dt = np.loadtxt(fileList_1[i])
   if len(dt[:,0])==Nmeas+1:
      rgyr_avg = rgyr_avg + np.sum(dt*dt, axis = 1)
      j = j+1
   #rgyr_avg = rgyr_avg + np.loadtxt(fileList_1[i])*np.loadtxt(fileList_1[i])
   
rgyr_avg = rgyr_avg/j

np.savetxt(endtoendFile, rgyr_avg)
print("\033[1;32mPrinted end to end squared tracks to '%s'\033[0m" %endtoendFile)
print(j)
##

j = 0
msd_avg3 = np.zeros((101))
for i in range(0,len(fileList_2)):
   dt = np.loadtxt(fileList_2[i])
   if len(dt)==101:
        msd_avg3 = msd_avg3 + dt
        j = j+1
   
msd_avg3 = msd_avg3/j
np.savetxt(msdFile , msd_avg3)
print(j)
##

msd_avg4 = np.zeros((Nchain,Nchain))
for i in range(0,len(fileList_3)):
   msd_avg4 = msd_avg4 + np.loadtxt(fileList_3[i])
   
msd_avg4 = msd_avg4/len(fileList_3)
np.savetxt(distmatFile , msd_avg4)
##

msd_avg5 = np.zeros((Nchain,Nchain))
for i in range(0,len(fileList_4)):
   msd_avg5 = msd_avg5 + np.loadtxt(fileList_4[i])
   
msd_avg5 = msd_avg5/len(fileList_4)
np.savetxt(cProbFile, msd_avg5)
##

j = 0
msd_avg6 = np.zeros((Nmeas+1))
for i in range(0,len(fileList_5)):
   dt = np.loadtxt(fileList_5[i])
   if len(dt)==Nmeas+1:
           msd_avg6 = msd_avg6 + np.loadtxt(fileList_5[i])
           j = j+1
msd_avg6 = msd_avg6/j
np.savetxt(rgyrFile , msd_avg6)
print(j)
print("\033[1;32mDone! \033[0m")



