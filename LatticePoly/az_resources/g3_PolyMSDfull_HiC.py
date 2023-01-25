##
# g3_cmMSD.py
# LatticePoly
##
# Based on mtortora's code
# Copyright © 2021 ENS Lyon. All rights reserved.
##

import os
import sys
import psutil

import numpy as np

from utils import msdFFT
from vtkReader import vtkReader


class g3_PolyMSD():

    def __init__(self, outputDir, initFrame):
        self.reader = vtkReader(outputDir, initFrame,
                                readLiq=False, readPoly=True)

        self.g3msdFile = os.path.join(
            self.reader.outputDir, "g3polyMSD_full.res")
        self.g2msdFile = os.path.join(
            self.reader.outputDir, "g2polyMSD_full.res")

        if os.path.exists(self.g3msdFile) & os.path.exists(self.g3msdFile):
            print("Files '%s' and '%s' already exist - aborting" %
                  (self.g3msdFile, self.g3msdFile))
            sys.exit()

    def Compute(self):
       
        self.DistHet = 0
        self.DistHom = 0

        g2, posCM = self.ReadHist()

        self.Distg2 = msdFFT(g2)
        self.DistCM = msdFFT(posCM)



    # def ComputeTad(self, idxTad):
    #	tadPosHist = np.zeros((self.reader.N, 3), dtype=np.float32)

    #	for i in range(self.reader.N):
    #		data = next(self.reader)
    #		tadPosHist[i] = data.polyPos[idxTad]

    #	self.distTad = msdFFT(tadPosHist)              

    def ReadHist(self):
        posCM  = np.zeros((self.reader.N, 3), dtype=np.float32)
        g2     = np.zeros((self.reader.N, self.reader.nTad, 3), dtype=np.float32)
        g2_r   = np.zeros((self.reader.N, 3), dtype=np.float32) 
        for i in range(self.reader.N):
            data = next(self.reader)
            g2[i]      = data.polyPos

            for idxTad in range(self.reader.nTad):
                    posCM[i] += data.polyPos[idxTad]
            g2[i]      = np.abs(g2[i] - posCM[i]/self.reader.nTad)         
            g2_r[i]    = np.mean(g2[i], axis = 0)

        posCM  = posCM/(self.reader.nTad)   


        return g2_r, posCM

    def Print(self):
            msdg2 = self.Distg2
            np.savetxt(self.g2msdFile, msdg2)

            print(
                "\033[1;32mPrinted g2 MSD to '%s'\033[0m" % self.g2msdFile)

            msdCM = self.DistCM
            np.savetxt(self.g3msdFile, msdCM)

            print("\033[1;32mPrinted CM MSD to '%s'\033[0m" %
                  self.g3msdFile)

    # def PrintTad(self, idxTad):
    #	msdFile = self.reader.outputDir + "/msdTad%05d.res" % idxTad
    #	np.savetxt(msdFile, self.distTad)

    #	print("\033[1;32mPrinted TAD MSD to '%s'\033[0m" % msdFile)


if __name__ == "__main__":
    if len(sys.argv) not in [3]:
        print("\033[1;31mUsage is %s outputDir initFrame \033[0m" %
              sys.argv[0])
        sys.exit()

    outputDir = sys.argv[1]
    initFrame = int(sys.argv[2])

    msd = g3_PolyMSD(outputDir, initFrame=initFrame)

    if len(sys.argv) == 3:
        msd.Compute()
        msd.Print()

    # elif len(sys.argv) == 4:
    #	idxTad = int(sys.argv[3])

    #	msd.ComputeTad(idxTad)
    #	msd.PrintTad(idxTad)
