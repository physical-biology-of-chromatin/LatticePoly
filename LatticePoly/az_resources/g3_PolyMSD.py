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

        self.g3msdHetFile = os.path.join(
            self.reader.outputDir, "g3polyHetMSD.res")
        self.g3msdHomFile = os.path.join(
            self.reader.outputDir, "g3polyHomMSD.res")

        if os.path.exists(self.g3msdHetFile) & os.path.exists(self.g3msdHomFile):
            print("Files '%s' and '%s' already exist - aborting" %
                  (self.g3msdHetFile, self.g3msdHomFile))
            sys.exit()

    def Compute(self):
        vMem = psutil.virtual_memory()
        sizeTot = self.reader.N * self.reader.polyPos.nbytes

        if sizeTot < vMem.available:
            self.DistHet = 0
            self.DistHom = 0

            posHet, posHomo = self.ReadHist()

            self.DistHet = msdFFT(posHet)
            self.DistHom = msdFFT(posHomo)

        else:
            print("Memory overflow likely - reduce chosen number of frames")
            sys.exit()

    # def ComputeTad(self, idxTad):
    #	tadPosHist = np.zeros((self.reader.N, 3), dtype=np.float32)

    #	for i in range(self.reader.N):
    #		data = next(self.reader)
    #		tadPosHist[i] = data.polyPos[idxTad]

    #	self.distTad = msdFFT(tadPosHist)              

    def ReadHist(self):
        posHet  = np.zeros((self.reader.N, 3), dtype=np.float32)
        posHomo = np.zeros((self.reader.N, 3), dtype=np.float32)

        for i in range(self.reader.N):
            data = next(self.reader)

            for idxTad in range(self.reader.nTad):
                if self.reader.polyType[idxTad] == 1:
                    posHet[i] += data.polyPos[idxTad]
                else:
                    posHomo[i] += data.polyPos[idxTad]

            #if (i+1) % 10 == 0:
            #    print("Read %d out of %d frames" % (i+1, self.reader.N))

        if self.reader.nHet > 0:
            posHet  = posHet/self.reader.nHet

        if self.reader.nEuc > 0:
            posHomo = posHomo/self.reader.nEuc            


        return posHet, posHomo

    def Print(self):
        if self.reader.nHet > 0:
            msdHet = self.DistHet
            np.savetxt(self.g3msdHetFile, msdHet)

            print(
                "\033[1;32mPrinted heterochromatic CM MSD to '%s'\033[0m" % self.g3msdHetFile)

        if self.reader.nEuc > 0:
            msdHom = self.DistHom
            np.savetxt(self.g3msdHomFile, msdHom)

            print("\033[1;32mPrinted euchromatic CM MSD to '%s'\033[0m" %
                  self.g3msdHomFile)

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
