##
# EndtoEnd.py
# LatticePoly
##
# Based on mtortora's code.
# Copyright © 2021 ENS Lyon. All rights reserved.
##

import os
import sys

import numpy as np

from vtkReader import vtkReader


class EndtoEnd():

    def __init__(self, outputDir, initFrame):
        self.reader = vtkReader(outputDir, initFrame, 
                                readLiq=False, readPoly=True)

        #self.anisoFile = os.path.join(self.reader.outputDir, "polyAniso.res")
        self.endtoendFile = os.path.join(self.reader.outputDir, "endtoend.res")

        if os.path.exists(self.endtoendFile):
            print("Files %s' already exist - aborting" % (self.endtoendFile))
            sys.exit()     

    def Compute(self, domainlen):

        if 3*domainlen > self.reader.nTad: 
            print("3 times domain length %d > polymer length %d - aborting" % (3*domainlen,self.reader.nTad))
            sys.exit() 

        #self.polyAniso = np.zeros((self.reader.N, self.reader.nDom), dtype=np.float32)
        self.endtoend = np.zeros((self.reader.N), dtype=np.float32)
        
        for i in range(self.reader.N):
            self.ProcessFrame(i, domainlen)

            if (i+1) % 10 == 0:
                print("Processed %d out of %d configurations" %
                      (i+1, self.reader.N))

    def ProcessFrame(self, i, domainlen):
        data = next(self.reader)
        rij  = np.zeros((self.reader.nTad - 3*domainlen,3), dtype=np.float32 )
        for j in range(self.reader.nTad - 3*domainlen):
               
                posI     = data.polyPos[j + domainlen]
                posN     = data.polyPos[j+2*domainlen-1]
                rij[j]   = posN - posI
        
        Rij = np.mean(np.sum(rij*rij, axis = 1))
        #self.polyAniso[i, id] = aniso
        self.endtoend[i] = Rij

    def Print(self):
        #np.savetxt(self.anisoFile, self.polyAniso)
        np.savetxt(self.endtoendFile, self.endtoend)

        print("\033[1;32mPrinted end to end vector squared to '%s'\033[0m" %self.endtoendFile)
        #print("\033[1;32mPrinted domain anisotropy factors to '%s'\033[0m" % self.anisoFile)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("\033[1;31mUsage is %s outputDir initFrame domainlen\033[0m" % sys.argv[0])
        sys.exit()

    outputDir = sys.argv[1]
    initFrame = int(sys.argv[2])
    domainlen = int(sys.argv[3])


    end = EndtoEnd(outputDir, initFrame=initFrame)

    end.Compute(domainlen)
    end.Print()
