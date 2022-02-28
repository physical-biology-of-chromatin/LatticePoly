##
# EndtoEnd.py -> MSD end to end vector
# LatticePoly
##
# Based on mtortora's code.
# Copyright © 2021 ENS Lyon. All rights reserved.
##

from az_resources.g1_PolyMSD_HiC import PolyMSD
import os
import sys

import numpy as np

from utils import msdFFT
from vtkReader import vtkReader


class msdEE():

    def __init__(self, outputDir, initFrame):
        self.reader = vtkReader(outputDir, initFrame, 
                                readLiq=False, readPoly=True)

        #self.anisoFile = os.path.join(self.reader.outputDir, "polyAniso.res")
        self.endtoendFile = os.path.join(self.reader.outputDir, "msdEEvect_160kb.res")

        if os.path.exists(self.endtoendFile):
            print("Files %s' already exist - aborting" % (self.endtoendFile))
            sys.exit()     

    def Compute(self, domainlen):

        if 3*domainlen > self.reader.nTad: 
            print("3 times domain length %d > polymer length %d - aborting" % (3*domainlen,self.reader.nTad))
            sys.exit() 

             
        self.distTad  = np.zeros((self.reader.N), dtype=np.float32)

        for j in range(self.reader.nTad - 3*domainlen):    
            PosHist       = np.zeros((self.reader.N, 3), dtype=np.float32)
            
            for i in range(self.reader.N):
                data = next(self.reader)
               
                posI         = data.polyPos[j + domainlen]
                posN         = data.polyPos[j+2*domainlen-1]
                PosHist[i]   = posN - posI
            
            self.reader.frame = initFrame
            self.distTad += (msdFFT(PosHist)) # / (2*np.mean(np.sum(PosHist*PosHist, axis = 1)))) 
        #self.polyAniso[i, id] = aniso
        self.distTad = self.distTad / (self.reader.nTad - 3*domainlen)

    def ComputeDomainend(self, domainstart, domainend):
        
        #self.polyAniso = np.zeros((self.reader.N, self.reader.nDom), dtype=np.float32)
        self.msdEndtoEnd = np.zeros((self.reader.N), dtype=np.float32)

        PosHist       = np.zeros((self.reader.N, 3), dtype=np.float32)
        for i in range(self.reader.N):
            data = next(self.reader)

            posI         = data.polyPos[domainstart]
            posN         = data.polyPos[domainend-1]
            PosHist[i]   = posN - posI    
        
        self.msdEndtoEnd = msdFFT(PosHist) #/ (2*np.mean(np.sum(PosHist*PosHist, axis = 1)))

    def PrintDomainlen(self):
        #np.savetxt(self.anisoFile, self.polyAniso)
        np.savetxt(self.endtoendFile, self.distTad)

        print("\033[1;32mPrinted msd of end to end vector (domain averaged) to '%s'\033[0m" %self.endtoendFile)
        #print("\033[1;32mPrinted domain anisotropy factors to '%s'\033[0m" % self.anisoFile)

    def PrintDomainend(self):
            #np.savetxt(self.anisoFile, self.polyAniso)
            np.savetxt(self.endtoendFile, self.msdEndtoEnd)

            print("\033[1;32mPrinted msd of end to end vector to '%s'\033[0m" %self.endtoendFile)
            #print("\033[1;32mPrinted domain anisotropy factors to '%s'\033[0m" % self.anisoFile)

if __name__ == "__main__":
    if len(sys.argv) not in [4,5]:
        print("\033[1;31mUsage is %s outputDir initFrame domainlen/start domainend\033[0m" % sys.argv[0])
        sys.exit()

    outputDir = sys.argv[1]
    initFrame = int(sys.argv[2])
    

    msd = msdEE(outputDir, initFrame=initFrame)

    if len(sys.argv) == 4:
        domainlen = int(sys.argv[3])
        msd.Compute(domainlen)
        msd.PrintDomainlen()

    if len(sys.argv) == 5:
        domainstart = int(sys.argv[3])
        domainend  = int(sys.argv[4])
        msd.ComputeDomainend(domainstart, domainend)
        msd.PrintDomainend()    


