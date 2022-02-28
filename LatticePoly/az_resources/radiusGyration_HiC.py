##
# radiusGyration.py
# LatticePoly
##
# Based on mtortora's code.
# Copyright © 2021 ENS Lyon. All rights reserved.
##

import os
import sys

import numpy as np

from vtkReader import vtkReader


class radiusGyration():

    def __init__(self, outputDir, initFrame):
        self.reader = vtkReader(outputDir, initFrame,
                                readLiq=False, readPoly=True)

        # self.anisoFile = os.path.join(self.reader.outputDir, "polyAniso.res")
        self.rgyrationFile = os.path.join(
            self.reader.outputDir, "rgyration.res")

        if os.path.exists(self.rgyrationFile):
            print("Files %s' already exist - aborting" % (self.rgyrationFile))
            sys.exit()

    def Compute(self, domainstart, domainend):
        # self.polyAniso = np.zeros((self.reader.N, self.reader.nDom), dtype=np.float32)
        self.rgyra = np.zeros(self.reader.N, dtype=np.float32)

        for i in range(self.reader.N):
            self.ProcessFrame(i, domainstart, domainend)

            if (i+1) % 10 == 0:
                print("Processed %d out of %d configurations" %
                      (i+1, self.reader.N))

    def ProcessFrame(self, i, domainstart, domainend):

        data   = next(self.reader)
        posd   = data.polyPos[domainstart:domainend]
        pos_cm = np.mean(posd, axis=0)
        posd   = posd - pos_cm 
        posd   = posd*posd
        pos    = np.sum(posd,axis =1)
        pos    = np.mean(pos)

        #for d in range(0, self.reader.nTad):
        #
        #    posd  = data.polyPos[d]
        #    posd -= pos_cm
        #    pos   = np.sum(posd*posd)
        
        self.rgyra[i] = pos

    def Print(self):
        # np.savetxt(self.anisoFile, self.polyAniso)
        np.savetxt(self.rgyrationFile, self.rgyra)

        print("\033[1;32mPrinted Rg square to '%s'\033[0m" %
              self.rgyrationFile)
        # print("\033[1;32mPrinted domain anisotropy factors to '%s'\033[0m" % self.anisoFile)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("\033[1;31mUsage is %s outputDir initFrame domainstart domainend\033[0m" % sys.argv[0])
        sys.exit()

    outputDir = sys.argv[1]
    initFrame = int(sys.argv[2])
    domainstart = int(sys.argv[3])
    domainend = int(sys.argv[4])

    rgyr = radiusGyration(outputDir, initFrame=initFrame)

    rgyr.Compute(domainstart, domainend)
    rgyr.Print()
