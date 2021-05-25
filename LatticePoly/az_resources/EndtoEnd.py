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

    def Compute(self):
        #self.polyAniso = np.zeros((self.reader.N, self.reader.nDom), dtype=np.float32)
        self.endtoend = np.zeros((self.reader.N, 3), dtype=np.float32)

        for i in range(self.reader.N):
            self.ProcessFrame(i)

            if (i+1) % 10 == 0:
                print("Processed %d out of %d configurations" %
                      (i+1, self.reader.N))

    def ProcessFrame(self, i):
        data = next(self.reader)

        posI = data.polyPos[0]
        posN = data.polyPos[np.max(self.reader.nTad)-1]
        rN   = posN - posI

        #self.polyAniso[i, id] = aniso
        self.endtoend[i] = rN

    def Print(self):
        #np.savetxt(self.anisoFile, self.polyAniso)
        np.savetxt(self.endtoendFile, self.endtoend)

        print("\033[1;32mPrinted end to end vector to '%s'\033[0m" %self.endtoendFile)
        #print("\033[1;32mPrinted domain anisotropy factors to '%s'\033[0m" % self.anisoFile)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("\033[1;31mUsage is %s outputDir initFrame\033[0m" % sys.argv[0])
        sys.exit()

    outputDir = sys.argv[1]
    initFrame = int(sys.argv[2])

    end = EndtoEnd(outputDir, initFrame=initFrame)

    end.Compute()
    end.Print()
