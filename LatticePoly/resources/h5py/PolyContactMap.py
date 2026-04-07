##
##  PolyContactMap.py
##  LatticePoly
##
##  Created by ppuel on 13/11/2025.
##  Based on a script by mtortora on 12/12/2019.
##
##  Copyright © 2025 ENS Lyon. All rights reserved.
##

import h5py, os, sys, numba, psutil

import math
import numpy as np

from hdf5Reader import hdf5Reader

# from scipy.spatial.distance import squareform



class ContactMap():

        def __init__(self, inputDir, outputDir, binSize, initFrame, cutoff=1/2**0.5 + 1e-3):
                print(f"ContactMap : Init {inputDir} {outputDir} {binSize}")
                self.reader = hdf5Reader(inputDir, "traj.h5", -1, readLiq=False, readPoly=True, backInBox=False)
                self.filePath = os.path.join(outputDir, "process.h5")

                self.binSize = int(binSize)
                self.cutoff = 5 * cutoff
                self.initFrame = int(initFrame)
                
                self.binNum = self.reader.nTad // self.binSize
                self.pruneNum = self.reader.nTad % self.binSize

                

        def Compute(self):
                print("\n")

                self.contactMap = np.zeros((self.binNum*(self.binNum-1)//2), dtype=np.float32)
                
                vMem = 4e9
                binSize_min = int(np.ceil(self.reader.nTad * 2 / vMem**0.5))
                print(binSize_min)
                if binSize >= binSize_min:
                        for i in range(self.reader.N):
                                self.ProcessFrame(i)
                        
                                if (i+1) % 10 == 0:
                                        print(" : Processed %d out of %d configurations" % (i+1, self.reader.N))

                else:
                        print("Memory overflow likely - increase chosen binSize (minimal value: %d)" % binSize_min)
                        sys.exit()
                        
                        
        def ProcessFrame(self, i):
                data = next(self.reader)
                
                if i>self.initFrame:
                        
                        polyPos = data.polyPos[:self.reader.nTad-self.pruneNum].astype(np.float32)
                        
                        self._sqDistPBC(polyPos, self.cutoff, self.binNum, self.binSize, self.contactMap)
                        

        def Print(self):
                
                with h5py.File(self.filePath, 'a') as hfile:
                        self.PrintDataset(hfile, "contactMap", self.contactMap)

        
        def PrintDataset(self, processFile, dataset_name, data, groupFolder = None, verbosity = True):
                
                tmpFolder = processFile[groupFolder] if groupFolder else processFile
                if dataset_name in list(tmpFolder.keys()):
                        del tmpFolder[dataset_name]
                
                tmpFolder.create_dataset(dataset_name, data = data)
            
                if verbosity:
                        print(f"Dataset {dataset_name} printed")

                           
        @staticmethod
        @numba.njit("void(f4[:,:], f4, i4, i4, f4[:])")
        def _sqDistPBC(pts, cutoff, binNum, binSize, contactMap):
                cnt = 0
                                
                for i in range(binNum-1):
                        for j in range(i+1, binNum):
                                for k in range(binSize):
                                        for l in range(binSize):
                                                pDist = 0.

                                                for m in range(3):
                                                        delta = pts[i*binSize+k, m] - pts[j*binSize+l, m]
                                                         
                                                        pDist += delta**2
                                                
                                                if pDist < cutoff**2:
                                                        contactMap[cnt] += 1
                                
                                cnt += 1

if __name__ == "__main__":
        if len(sys.argv) != 5:
                print("\033[1;31mUsage is %s inputDir outputDir binSize initFrame \033[0m" % sys.argv[0])
                sys.exit()

        inputDir = sys.argv[1]
        outputDir = sys.argv[2]
        binSize = int(sys.argv[3])
        initFrame = int(sys.argv[4])
        
        distMap = ContactMap(inputDir, outputDir, binSize=binSize, initFrame=initFrame)
        
        distMap.Compute()
        distMap.Print()
