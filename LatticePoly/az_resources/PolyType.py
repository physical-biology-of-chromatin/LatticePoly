##
##  PolyType.py
##  LatticePoly
##
##  Created by amithzafal based on mtortora's code on 02/03/2022.
##  Copyright © 2022 ENS Lyon. All rights reserved.
##

import os
import sys

import numpy as np

from vtkReader import vtkReader


class PolyType():
    
    def __init__(self, outputDir, initFrame):
        self.reader = vtkReader(outputDir, initFrame, readLiq=False, readPoly=True)	
        self.typeFile = os.path.join(self.reader.outputDir, "polyType.res")

        if os.path.exists(self.typeFile):
            print("Files '%s' already exist - aborting" % (self.typeFile))
            sys.exit()

    def Compute(self):
        
        self.polyType = np.zeros(self.reader.nTad, dtype=np.float32)
        j = 0
        for i in range(self.reader.N):
            data = next(self.reader)
            self.polyType = self.polyType + data.polyType
            j = j + 1 
        
        self.polyType = self.polyType/j
        
    def Print(self):

        np.savetxt(self.typeFile, self.polyType)
        
        print("\033[1;32mPrinted Polytype to '%s'\033[0m" % self.typeFile)


        
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("\033[1;31mUsage is %s outputDir initFrame\033[0m" % sys.argv[0])
        sys.exit()

    outputDir = sys.argv[1]
    initFrame = int(sys.argv[2])

    gyr = PolyType(outputDir, initFrame=initFrame)

    gyr.Compute()
    gyr.Print()
