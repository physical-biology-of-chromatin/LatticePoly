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
            self.reader.outputDir, "g3polyMSD_1.res")

        self.g2msdHetFile = os.path.join(
            self.reader.outputDir, "g2polyMSD_1.res")    


        if os.path.exists(self.g3msdHetFile):
            print("Files '%s' and '%s' already exist - aborting" %
                  (self.g3msdHetFile))
            sys.exit()

    def Compute(self, domainstart, domainend):

            self.DistHet = 0
            self.Distg2  = 0
        

            posHet,g2 = self.ReadHist(domainstart, domainend)

            self.DistHet = msdFFT(posHet)
            self.Distg2  = msdFFT(g2)

    # def ComputeTad(self, idxTad):
    #	tadPosHist = np.zeros((self.reader.N, 3), dtype=np.float32)

    #	for i in range(self.reader.N):
    #		data = next(self.reader)
    #		tadPosHist[i] = data.polyPos[idxTad]

    #	self.distTad = msdFFT(tadPosHist)              

    def ReadHist(self, domainstart, domainend):
        posHet  = np.zeros((self.reader.N, 3), dtype=np.float32)
        g2      = np.zeros((self.reader.N, domainend-domainstart, 3), dtype=np.float32)
        g2_r    = np.zeros((self.reader.N, 3), dtype=np.float32)
        for i in range(self.reader.N):
            data = next(self.reader)
            g2[i]      = data.polyPos[domainstart:domainend]

            for idxTad in range(domainstart, domainend):
                    posHet[i] += data.polyPos[idxTad]
            g2[i]      = np.abs(g2[i] - posHet[i]/(domainend-domainstart))    
            g2_r[i]    = np.mean(g2[i], axis = 0)     

            #if (i+1) % 10 == 0:
            #    print("Read %d out of %d frames" % (i+1, self.reader.N))
        posHet  = posHet/(domainend-domainstart)   
        

        return posHet,g2_r

    def Print(self):

            msdHet = self.DistHet
            msdg2  = self.Distg2
            np.savetxt(self.g3msdHetFile, msdHet)
            np.savetxt(self.g2msdHetFile, msdg2)
            print("\033[1;32mPrinted heterochromatic CM MSD to '%s'\033[0m" % self.g3msdHetFile)
            print("\033[1;32mPrinted heterochromatic g2 MSD to '%s'\033[0m" % self.g2msdHetFile)

    # def PrintTad(self, idxTad):
    #	msdFile = self.reader.outputDir + "/msdTad%05d.res" % idxTad
    #	np.savetxt(msdFile, self.distTad)

    #	print("\033[1;32mPrinted TAD MSD to '%s'\033[0m" % msdFile)


if __name__ == "__main__":
    if len(sys.argv) not in [5]:
        print("\033[1;31mUsage is %s outputDir initFrame domainstart domainend \033[0m" %
              sys.argv[0])
        sys.exit()

    outputDir = sys.argv[1]
    initFrame = int(sys.argv[2])
    domainstart = int(sys.argv[3])
    domainend = int(sys.argv[4])

    msd = g3_PolyMSD(outputDir, initFrame=initFrame)

    if len(sys.argv) == 5:
        msd.Compute(domainstart,domainend)
        msd.Print()

 
