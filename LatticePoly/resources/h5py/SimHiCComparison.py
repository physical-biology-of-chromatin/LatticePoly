##
##  NormalizeData.py
##  LatticePoly
##
##  Created by ppuel on 2/12/2025.
##  Copyright © 2025 ENS Lyon. All rights reserved.
##

import os
import sys

import numpy as np

import h5py

from hdf5Reader import hdf5Reader

import numba

class HiCComparison():

        def __init__(self, inputDir, outputDir):
                print(f"HiCComparison : Init {inputDir} {outputDir}")
                
                initPath = os.path.join(inputDir, "N/0")
                aggPath = os.path.join(inputDir, "aggregated_process.h5")

                self.inputPath = os.path.join(inputDir, "normalized_process.h5")
                self.outputPath = os.path.join(outputDir, "normalized_process.h5")
                
                with h5py.File(aggPath, 'r') as aggFile:
                        self.binNum = np.int32(len(aggFile["contactProb"]) + 1)
                
                reader = hdf5Reader(initPath, 'traj.h5', initFrame=-1, readLiq=False)
                
                nTad = reader.nTad
                binSize = nTad // self.binNum
                pruneNum = nTad % binSize
                
                polyType = reader.polyType[:reader.nTad-pruneNum]
                 
                reader.Close() 
                

                averageType = np.mean(np.reshape(polyType, (self.binNum, binSize)), axis = 1)

                self.boolType = averageType > 0
                sumType = np.sum(self.boolType)

                self.pcgContact = np.zeros((sumType)*(sumType-1)//2, dtype=np.float32)
                self.interContact = np.zeros((sumType)*(nTad - sumType), dtype=np.float32)


        def Compute(self):
                print("\n")
                with h5py.File(self.inputPath, 'r') as normFile:
                        
                        contactHiC = np.array(normFile["contactHiC"][:]).astype(np.float64)
                        
                        self._contactPCG(contactHiC, self.boolType, self.binNum, self.pcgContact, self.interContact)
                                        
                           
        @staticmethod
        @numba.njit("void(f8[:], b1[:], i4, f4[:], f4[:])")
        def _contactPCG(contactHiC, boolType, binNum, pcgContact, interContact):
                
                cntPCG = 0
                cntInter = 0
                cnt = 0
                for i in range(binNum-1):
                        for j in range(i+1, binNum):
                                if np.logical_and(boolType[i], boolType[j]):
                                        pcgContact[cntPCG] = contactHiC[cnt]
                                        cntPCG += 1
                                        cnt += 1
                                elif np.logical_or(boolType[i], boolType[j]):
                                        interContact[cntInter] = contactHiC[cnt]
                                        cntInter += 1
                                        cnt += 1
                                else:
                                        cnt += 1



        def Print(self):
                with h5py.File(self.outputPath, 'a') as normFile:
                        self.PrintDataset(normFile, "pcgContact", self.pcgContact)
                        self.PrintDataset(normFile, "interContact", self.interContact)
                
        def PrintDataset(self, processFile, dataset_name, data):
                
                if dataset_name in processFile.keys():
                        tmp = processFile[dataset_name]
                        tmp[:] = data
                else:
                        processFile.create_dataset(dataset_name, data = data)
                        
                print(f"Dataset {dataset_name} printed")
         

if __name__ == "__main__":
        if len(sys.argv) != 3:
                print("\033[1;31mUsage is %s inputDir outputDir\033[0m" % sys.argv[0])
                sys.exit()

        inputDir = sys.argv[1]
        outputDir = sys.argv[2]

        hiC = HiCComparison(inputDir, outputDir)

        hiC.Compute()
        hiC.Print()
        
        print("\n")
        print("NormalizeData : Done\n\n")
