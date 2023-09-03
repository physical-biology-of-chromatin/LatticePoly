##
##  LiqPolyCoincidence.py
##  LatticePoly
##
##  Created by mtortora on 17/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import os
import sys

import numpy as np

from vtkReader import vtkReader
from scipy.spatial import cKDTree


class LiqPolyCoincidence():

    def __init__(self, outputDir, initFrame, cutoff=1e-3):
        self.reader = vtkReader(outputDir, initFrame, readLiq=True, readPoly=True, backInBox=True)

        self.cutoff = cutoff

        self.liq_Het_File = os.path.join(self.reader.outputDir, "liq_Het_Coincidence.res")
        self.poly_Het_File = os.path.join(self.reader.outputDir, "poly_Het_Coincidence.res")
        self.liq_PRE_File = os.path.join(self.reader.outputDir, "liq_PRE_Coincidence.res")
        self.poly_PRE_File = os.path.join(self.reader.outputDir, "poly_PRE_Coincidence.res")


        if os.path.exists(self.liq_Het_File) & os.path.exists(self.poly_Het_File) & os.path.exists(self.liq_PRE_File) & os.path.exists(self.poly_PRE_File):
            print("Files '%s' & '%s' & '%s' and '%s' already exist - aborting" % (self.liq_Het_File, self.poly_Het_File, self.liq_PRE_File, self.poly_PRE_File))
            sys.exit()


    def Compute(self):
        self.liq_Het_Cont = np.zeros(self.reader.N, dtype=np.float32)
        self.poly_Het_Cont = np.zeros(self.reader.N, dtype=np.float32)
        self.liq_PRE_Cont = np.zeros(self.reader.N, dtype=np.float32)
        self.poly_PRE_Cont = np.zeros(self.reader.N, dtype=np.float32)

        for i in range(self.reader.N):
            self.ProcessFrame(i)

        if (i+1) % 10 == 0:
            print("Processed %d out of %d configurations" % (i+1, self.reader.N))


    def ProcessFrame(self, i):
        data = next(self.reader)
        nPRE = np.count_nonzero(self.reader.polyPainter == 1)

        if nPRE > 0:
            liqTree = cKDTree(data.liqPos, boxsize=data.boxDim)
            polyTree = cKDTree(data.polyPos[data.polyPainter == 1], boxsize=data.boxDim)

            liqPolyIds = liqTree.query_ball_tree(polyTree, self.cutoff)
            liqPolyIds = list(filter(None, liqPolyIds))

            polyIds = [id for ids in liqPolyIds for id in ids]
            polyIds = set(polyIds)

            numLiqCont = len(liqPolyIds)
            numPolyCont = len(polyIds)

            self.liq_PRE_Cont[i] = numLiqCont / float(data.nLiq)
            self.poly_PRE_Cont[i] = numPolyCont / float(nPRE)
        else:
            self.liq_PRE_Cont[i] = 0.
            self.poly_PRE_Cont[i] = 0.


        if data.nHet > 0:
            liqTree = cKDTree(data.liqPos, boxsize=data.boxDim)
            polyTree = cKDTree(data.polyPos[data.polyType == 1], boxsize=data.boxDim)

            liqPolyIds = liqTree.query_ball_tree(polyTree, self.cutoff)
            liqPolyIds = list(filter(None, liqPolyIds))

            polyIds = [id for ids in liqPolyIds for id in ids]
            polyIds = set(polyIds)

            numLiqCont = len(liqPolyIds)
            numPolyCont = len(polyIds)

            self.liq_Het_Cont[i] = numLiqCont / float(data.nLiq)
            self.poly_Het_Cont[i] = numPolyCont / float(data.nHet)

        else:
            self.liq_Het_Cont[i] = 0.
            self.poly_Het_Cont[i] = 0.


    def Print(self):
        np.savetxt(self.liq_Het_File, self.liq_Het_Cont)
        np.savetxt(self.poly_Het_File, self.poly_Het_Cont)
        np.savetxt(self.liq_PRE_File, self.liq_PRE_Cont)
        np.savetxt(self.poly_PRE_File, self.poly_PRE_Cont)

        print("\033[1;32mPrinted mean liquid Coincidence fractions to '%s'\033[0m" % self.liq_Het_File)
        print("\033[1;32mPrinted mean polymer Coincidence fractions to '%s'\033[0m" % self.poly_Het_File)
        print("\033[1;32mPrinted mean liquid Coincidence fractions to '%s'\033[0m" % self.liq_PRE_File)
        print("\033[1;32mPrinted mean polymer Coincidence fractions to '%s'\033[0m" % self.poly_PRE_File)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("\033[1;31mUsage is %s outputDir initFrame\033[0m" % sys.argv[0])
        sys.exit()

    outputDir = sys.argv[1]
    initFrame = int(sys.argv[2])

    Coincidence = LiqPolyCoincidence(outputDir, initFrame=initFrame)

    Coincidence.Compute()
    Coincidence.Print()
