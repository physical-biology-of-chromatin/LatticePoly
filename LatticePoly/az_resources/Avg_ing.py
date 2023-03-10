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



dpath = '/home/aabdulla/3D/TopoSwap/LatticePoly/LatticePoly/data/output/Twochain_Topo/DenseChain_1016x2/'
Nchain = 2032
Nmeas = 1000

def fetchFiles(dirPath, fileName):
	fileList = []

	for root, dirs, files in os.walk(dirPath):
		if not dirs:
			if fileName in files:
				filePath = os.path.join(root, fileName)
				fileList.append(filePath)
	
	return fileList

fileList_1 = fetchFiles(dpath, 'r_0_endtoend.res')
fileList_11 = fetchFiles(dpath, 'r_1016_endtoend.res')
fileList_2 = fetchFiles(dpath, 'r0msdEEvect.res')
fileList_22 = fetchFiles(dpath, 'r1016msdEEvect.res')
fileList_3 = fetchFiles(dpath, 'r_monomerdistmatTopo.res')
fileList_4 = fetchFiles(dpath, 'r_contactProbTopo.res')
fileList_5 = fetchFiles(dpath, 'r0gyration.res')
fileList_55 = fetchFiles(dpath, 'r1016gyration.res')


endtoendFile1   = os.path.join(dpath, "endtoend_r0_Avg.res")
endtoendFile2   = os.path.join(dpath, "endtoend_r1016_Avg.res")

msdFile1        = os.path.join(dpath, "msdEE_r0_Avg.res")
msdFile2        = os.path.join(dpath, "msdEE_r1016_Avg.res")

distmatFile    = os.path.join(dpath, "distmat_Avg.res") 
cProbFile      = os.path.join(dpath, "contactProb_Avg.res")

rgyrFile1       = os.path.join(dpath, "rgyr_r0_Avg.res")
rgyrFile2       = os.path.join(dpath, "rgyr_r1016_Avg.res")

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

np.savetxt(endtoendFile1, rgyr_avg)
print("\033[1;32mChain1:Printed end to end squared tracks to '%s'\033[0m" %endtoendFile1)
print(j)
##

##                                                                                                                                            
rgyr_avg = np.zeros((Nmeas+1))
j = 0
for i in range(0,len(fileList_11)):
   dt = np.loadtxt(fileList_11[i])
   if len(dt[:,0])==Nmeas+1:
      rgyr_avg = rgyr_avg + np.sum(dt*dt, axis = 1)
      j = j+1
   #rgyr_avg = rgyr_avg + np.loadtxt(fileList_1[i])*np.loadtxt(fileList_1[i])                                                                   

rgyr_avg = rgyr_avg/j

np.savetxt(endtoendFile2, rgyr_avg)
print("\033[1;32mChain2:Printed end to end squared tracks to '%s'\033[0m" %endtoendFile2)
print(j)
## 

j = 0
msd_avg3 = np.zeros((Nmeas+1))
for i in range(0,len(fileList_2)):
   dt = np.loadtxt(fileList_2[i])
   if len(dt)==Nmeas+1:
        msd_avg3 = msd_avg3 + dt
        j = j+1
   
msd_avg3 = msd_avg3/j
np.savetxt(msdFile1 , msd_avg3)
print("\033[1;32mChain1:Printed MSD of end to end vector to '%s'\033[0m" %msdFile1)
print(j)
##

j = 0
msd_avg3 = np.zeros((Nmeas+1))
for i in range(0,len(fileList_22)):
   dt = np.loadtxt(fileList_22[i])
   if len(dt)==Nmeas+1:
        msd_avg3 = msd_avg3 + dt
        j = j+1

msd_avg3 = msd_avg3/j
np.savetxt(msdFile2 , msd_avg3)
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
print("\033[1;32mPrinted distmat and contactProb \033[0m")
##

j = 0
msd_avg6 = np.zeros((Nmeas+1))
for i in range(0,len(fileList_5)):
   dt = np.loadtxt(fileList_5[i])
   if len(dt)==Nmeas+1:
           msd_avg6 = msd_avg6 + dt
           j = j+1
msd_avg6 = msd_avg6/j
np.savetxt(rgyrFile1 , msd_avg6)
print("\033[1;32mChain1:Printed Rg squared tracks to '%s'\033[0m" %rgyrFile1)
print(j)

j = 0
msd_avg6 = np.zeros((Nmeas+1))
for i in range(0,len(fileList_55)):
   dt = np.loadtxt(fileList_55[i])
   if len(dt)==Nmeas+1:
           msd_avg6 = msd_avg6 + dt
           j = j+1
msd_avg6 = msd_avg6/j
np.savetxt(rgyrFile2 , msd_avg6)
print(j)




