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

from scipy.spatial import cKDTree

from vtkReader import vtkReader


class MonomerDmap():

    def __init__(self, outputDir, initFrame):
        self.reader = vtkReader(outputDir, initFrame,
                                readLiq=False, readPoly=True, backInBox=True)

        #self.anisoFile = os.path.join(self.reader.outputDir, "polyAniso.res")
        self.monomerFile = os.path.join(self.reader.outputDir, "monomerdistmat.res")
        self.contactFile = os.path.join(self.reader.outputDir, "contactProb.res")

        if os.path.exists(self.monomerFile):
            print("Files %s' already exist - aborting" % (self.monomerFile))
            sys.exit()

    def Compute(self):
        #self.polyAniso = np.zeros((self.reader.N, self.reader.nDom), dtype=np.float32)
        self.monomer     = np.zeros((self.reader.nTad, self.reader.nTad), dtype=np.float32)
        self.contactProb = np.zeros((self.reader.nTad, self.reader.nTad), dtype=np.float32)
  

        for i in range(initFrame, self.reader.N):
            self.ProcessFrame(i)

            if (i+1) % 10 == 0:
                print("Processed %d out of %d configurations" %
                      (i+1, self.reader.N))

          

    def ProcessFrame(self, i):
        data = next(self.reader)

        tree1   = cKDTree(data.polyPos, boxsize = data.boxDim)
        sparse1 = tree1.sparse_distance_matrix(tree1,data.boxDim[0]) 
        sparse1 = sparse1.toarray()
        self.monomer = self.monomer + sparse1
        
        sparse2 = tree1.sparse_distance_matrix(tree1,1) #cutoff 1 Lu, for nearest neighbour FCC it should be (1/root(2))*1Lu
        sparse2 = sparse2.toarray()
        sparse2[np.where(sparse2 != 0)] = 1
        self.contactProb = self.contactProb + sparse2
     


    def Print(self):
        #np.savetxt(self.anisoFile, self.polyAniso)
        np.savetxt(self.monomerFile, self.monomer )
        np.savetxt(self.contactFile, self.contactProb )

        print("\033[1;32mPrinted avg.contact probability (normalize with number of frames) to '%s'\033[0m" %self.contactFile)
        print("\033[1;32mPrinted distance map (normalize with number of frames) to '%s'\033[0m" % self.monomerFile)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("\033[1;31mUsage is %s outputDir initFrame\033[0m" % sys.argv[0])
        sys.exit()

    outputDir = sys.argv[1]
    initFrame = int(sys.argv[2])

    monom = MonomerDmap(outputDir, initFrame=initFrame)

    monom.Compute()
    monom.Print()

