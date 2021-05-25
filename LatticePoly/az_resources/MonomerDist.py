##
# MonomerDist.py
# LatticePoly
##
# Based on mtortora's code.
# Copyright © 2021 ENS Lyon. All rights reserved.
##

import os
import sys

import numpy as np

from vtkReader import vtkReader


class MonomerMonomer():

    def __init__(self, outputDir, initFrame):
        self.reader = vtkReader(outputDir, initFrame,
                                readLiq=False, readPoly=True)

        #self.anisoFile = os.path.join(self.reader.outputDir, "polyAniso.res")
        self.monomerFile = os.path.join(self.reader.outputDir, "monomerdist.res")

        if os.path.exists(self.monomerFile):
            print("Files %s' already exist - aborting" % (self.monomerFile))
            sys.exit()

    def Compute(self):
        #self.polyAniso = np.zeros((self.reader.N, self.reader.nDom), dtype=np.float32)
        self.monomer = np.zeros((self.reader.N, self.reader.nTad), dtype=np.float32)

        for i in range(self.reader.N):
            self.ProcessFrame(i)

            if (i+1) % 10 == 0:
                print("Processed %d out of %d configurations" %
                      (i+1, self.reader.N))

    def ProcessFrame(self, i):
        data = next(self.reader)

        for j in range(0,self.reader.nTad):
            for k in range(0,self.reader.nTad):
                dist = data.polyPos[k] - data.polyPos[j]
                self.monomer[i,np.abs(k-j)] = np.sum(dist*dist) +  self.monomer[i,np.abs(k-j)] #Square of the monomer - monomer distance



    def Print(self):
        #np.savetxt(self.anisoFile, self.polyAniso)
        np.savetxt(self.monomerFile, self.monomer)

        print("\033[1;32mPrinted squared monomer-monomer distance with time to '%s'\033[0m" %self.monomerFile)
        #print("\033[1;32mPrinted domain anisotropy factors to '%s'\033[0m" % self.anisoFile)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("\033[1;31mUsage is %s outputDir initFrame\033[0m" % sys.argv[0])
        sys.exit()

    outputDir = sys.argv[1]
    initFrame = int(sys.argv[2])

    monom = MonomerMonomer(outputDir, initFrame=initFrame)

    monom.Compute()
    monom.Print()

