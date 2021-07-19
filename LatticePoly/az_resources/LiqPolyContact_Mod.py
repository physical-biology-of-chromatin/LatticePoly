# -*- coding: utf-8 -*-
##
# LiqPolyContact_Mod.py
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


class LiqContact():

    def __init__(self, outputDir, initFrame, cutoff = 1/2**0.5 + 1e-3 ):
        self.reader = vtkReader(outputDir, initFrame,
                                readLiq=True, readPoly=True, backInBox=True)

        self.cutoff = cutoff

        #self.anisoFile = os.path.join(self.reader.outputDir, "polyAniso.res")
        self.contactFile = os.path.join(self.reader.outputDir, "ContactFracLiq.res")

        if os.path.exists(self.contactFile):
            print("Files %s' already exist - aborting" % (self.contactFile))
            sys.exit()

    def Compute(self):
        #self.polyAniso = np.zeros((self.reader.N, self.reader.nDom), dtype=np.float32)
        self.contacts = np.zeros(self.reader.nTad, dtype=np.float32)
        

        for i in range(initFrame, self.reader.N):
            self.ProcessFrame(i)

            if (i+1) % 10 == 0:
                print("Processed %d out of %d configurations" %
                      (i+1, self.reader.N))

        self.contacts = self.contacts/( float(self.reader.nLiq)*(self.reader.N-initFrame) )


    def ProcessFrame(self, i):
        
        data = next(self.reader)

        liqTree = cKDTree(data.liqPos, boxsize=data.boxDim)
        polyTree = cKDTree(data.polyPos, boxsize=data.boxDim)

        polyPolyIds = polyTree.query_ball_tree(liqTree, self.cutoff)
		
        
        for i in range(0,len(polyPolyIds)):
            self.contacts[i] = len(polyPolyIds[i]) + self.contacts[i] 
			
        
    def Print(self):
        #np.savetxt(self.anisoFile, self.polyAniso)
        np.savetxt(self.contactFile, self.contacts )

        print("\033[1;32mPrinted avg.fraction of liqbinders in contact with the polymer to '%s'\033[0m" %self.contactFile)
       


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("\033[1;31mUsage is %s outputDir initFrame\033[0m" % sys.argv[0])
        sys.exit()

    outputDir = sys.argv[1]
    initFrame = int(sys.argv[2])

    monom = LiqContact(outputDir, initFrame=initFrame)

    monom.Compute()
    monom.Print()

