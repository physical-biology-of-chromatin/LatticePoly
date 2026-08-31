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

class NormalizeData():

        def __init__(self, inputDir, outputDir):
                print(f"NormalizeData : Init {inputDir} {outputDir}")
                self.inputPath = os.path.join(inputDir, "aggregated_process.h5")
                self.outputPath = os.path.join(outputDir, "normalized_process.h5")
                

        def Compute(self):
                print("\n")
                with h5py.File(self.inputPath, 'r') as aggfile:
                        baseline = np.array(aggfile["contactProb"][:])
                        self.datasetArray = np.array(aggfile["contactHiC"][:])

                        binNum = len(baseline) + 1
                                
                        if len(self.datasetArray) == (binNum)*(binNum-1)//2:
                                        cnt = 0
                                
                                        for i in range(binNum-1):
                                                for j in range(i+1, binNum):
                                                        self.datasetArray[cnt] /= (baseline[j-i-1]/(binNum-(j-i))) if baseline[j-i-1] > 0 else 1
                                                        cnt += 1
                 
       
        def Print(self):
                with h5py.File(self.outputPath, 'a') as normFile:
                        self.PrintDataset(normFile, "contactHiC", self.datasetArray)
                
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

        norm = NormalizeData(inputDir, outputDir)

        norm.Compute()
        norm.Print()
        
        print("\n")
        print("NormalizeData : Done\n\n")
