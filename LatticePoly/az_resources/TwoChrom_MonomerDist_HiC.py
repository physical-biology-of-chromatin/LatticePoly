# -*- coding: utf-8 -*-
##
# MonomerDist.py
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

from scipy.spatial.distance import pdist, squareform


class MonomerDmap():

    def __init__(self, outputDir, chrom1, chrom2, chrom3 initFrame):
        self.reader1 = vtkReader(outputDir, chrom1, initFrame,
                                readLiq=False, readPoly=True)
        self.reader2 = vtkReader(outputDir, chrom2, initFrame,
                                readLiq=False, readPoly=True)
        self.reader3 = vtkReader(outputDir, chrom3, initFrame,
                                readLiq=False, readPoly=True)                        

        #self.anisoFile = os.path.join(self.reader.outputDir, "polyAniso.res")
        self.monomerFile = os.path.join(self.reader1.outputDir, "r"+str(initFrame)+"_distmatTopo.res")
        self.contactFile = os.path.join(self.reader1.outputDir, "r"+str(initFrame)+"_cProbTopo.res")

        if os.path.exists(self.monomerFile):
            print("Files %s' already exist - aborting" % (self.monomerFile))
            sys.exit()

    def Compute(self, domainstart, domainend):
        #self.polyAniso = np.zeros((self.reader.N, self.reader.nDom), dtype=np.float32)
        self.monomer     = np.zeros((domainend-domainstart, domainend-domainstart), dtype=np.float32)
        self.contactProb = np.zeros((domainend-domainstart, domainend-domainstart), dtype=np.float32)
        

        #for i in range(initFrame, initFrame+11):
        self.ProcessFrame(domainstart, domainend)

            #if (i+1) % 100 == 0:
            #    print("Processed %d out of %d configurations" %
            #          (i+1, self.reader1.N))

        #self.monomer     = self.monomer/11
        #self.contactProb = self.contactProb/11
        np.fill_diagonal(self.contactProb, self.contactProb.diagonal() + 1)  

    def ProcessFrame(self, domainstart, domainend):
        data1 = next(self.reader1)
        data2 = next(self.reader2)
        data3 = next(self.reader3)
        
        d1 = data1.polyPos
        d2 = data2.polyPos
        d3 = data3.polyPos
        pos = np.concatenate((d1, d2, d3), axis=0)

        dists       = pdist(pos[domainstart:domainend])
        distanceMap = squareform(dists)
        self.monomer = self.monomer + distanceMap
        
        tree1   = cKDTree(pos[domainstart:domainend], boxsize = None)
        pairs = tree1.query_pairs(r = 1) # NN distance FCC lattice 1/np.sqrt(2) = 0.71
        for (i,j) in pairs:
            self.contactProb[i,j] = self.contactProb[i,j] + 1



    def Print(self):
        #np.savetxt(self.anisoFile, self.polyAniso)
        np.savetxt(self.monomerFile, self.monomer )
        np.savetxt(self.contactFile, self.contactProb )

        print("\033[1;32mPrinted avg.contact probability to '%s'\033[0m" %self.contactFile)
        print("\033[1;32mPrinted distance map to '%s'\033[0m" % self.monomerFile)


if __name__ == "__main__":
    if len(sys.argv) != 8:
        print("\033[1;31mUsage is %s outputDir chrom1 chrom2 chrom3 initFrame domainstart domainend\033[0m" % sys.argv[0])
        sys.exit()

    outputDir = sys.argv[1]
    chrom1    = sys.argv[2]
    chrom2    = sys.argv[3]
    chrom3    = sys.argv[4]
    initFrame = int(sys.argv[5])
    domainstart = int(sys.argv[6])
    domainend = int(sys.argv[7])
    monom = MonomerDmap(outputDir, chrom1, chrom2, chrom3 initFrame=initFrame)

    monom.Compute(domainstart, domainend)
    monom.Print()

