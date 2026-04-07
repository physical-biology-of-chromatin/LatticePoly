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
import h5py

import numpy as np

from utils import msdFFT
from hdf5Reader import hdf5Reader

class PolyMSD():
    def __init__(self, inputDir : str):
            
        print(f"PolyMSD : Init {inputDir}")
        self.reader = hdf5Reader(inputDir, 'traj.h5', -1, read_liq=False, read_poly=True)
        self.processPath = os.path.join(inputDir, "process.h5")
        

    def Compute(self):
        vMem = psutil.virtual_memory()
        sizeTot = self.reader.n_frame * self.reader.poly_pos.nbytes
        
        if sizeTot < vMem.available:
            self.cumulDistHet = 0
            self.cumulDistHom = 0
            self.cumulDistPRE = 0
        
            posHist = self.ReadHist()
            
            for idxTad in range(self.reader.n_tad):
                if self.reader.poly_painter[idxTad] == 1:
                    self.cumulDistPRE += msdFFT(posHist[:, idxTad])
                elif self.reader.poly_type[idxTad] == 1:
                    self.cumulDistHet += msdFFT(posHist[:, idxTad])
                else:
                    self.cumulDistHom += msdFFT(posHist[:, idxTad])
                            
                if (idxTad+1) % 1000 == 0:
                    print("Processed %d out of %d TADs" % (idxTad+1, self.reader.nTad))
                    
        else:
            print("Memory overflow likely - reduce chosen number of frames")
            sys.exit()


    def ComputeTad(self, idxTad):
        tadPosHist = np.zeros((self.reader.n_frame, 3), dtype=np.float32)

        for i in range(self.reader.n_frame):
            data = next(self.reader)
            tadPosHist[i] = data.poly_pos[idxTad]
            
        self.distTad = msdFFT(tadPosHist)

    def ReadHist(self):
        posHist = np.zeros((self.reader.n_frame, self.reader.n_tad, 3), dtype=np.float32)
        
        for i in range(self.reader.n_frame):
            data = next(self.reader)
            posHist[i] = data.poly_pos
            
        return posHist
    

    def Print(self):
        print("\n")
        
        self.reader.close()
        
        with h5py.File(self.processPath, 'a') as processFile:
            if np.count_nonzero(self.reader.poly_painter == 1) > 0:
                msdPRE = self.cumulDistPRE / np.count_nonzero(self.reader.poly_painter == 1)
                self.PrintDataset(processFile, "polyPREMSD", data = msdPRE)
            
            if self.reader.n_het > 0:
                msdHet = self.cumulDistHet /  self.reader.n_het
                self.PrintDataset(processFile, "polyHetMSD", data = msdHet)

            if self.reader.n_euc > 0:
                msdHom = self.cumulDistHom / self.reader.n_euc
                self.PrintDataset(processFile, "polyHomMSD", data = msdHom)
    

    def PrintTad(self, idxTad):
        print("\n")
        
        self.reader.close()
        
        with h5py.File(self.processPath, 'a') as processFile:
            self.PrintDataset(processFile, "msdTad%05d" % idxTad, data = self.distTad)
    

    def PrintDataset(self, processFile, dataset_name, data):
    
        if dataset_name in processFile.keys():
            del processFile[dataset_name]
            
        else:
            processFile.create_dataset(dataset_name, data = data)
            
        print(f"Dataset {dataset_name} printed")


if __name__ == "__main__":
    if len(sys.argv) not in [2, 3]:
        print("\033[1;31mUsage is %s inputDir [idxTad]\033[0m" % sys.argv[0])
        sys.exit()

    inputDir = sys.argv[1]

    msd = PolyMSD(inputDir)

    if len(sys.argv) == 2:
        msd.Compute()
        msd.Print()
    
        print("\n")
        print("PolyMSD : Done\n\n")
        
    elif len(sys.argv) == 3:
        idxTad = int(sys.argv[2])
    
        msd.ComputeTad(idxTad)
        msd.PrintTad(idxTad)
    
        print("\n")
        print("PolyMSD[idxTad] : Done\n\n")