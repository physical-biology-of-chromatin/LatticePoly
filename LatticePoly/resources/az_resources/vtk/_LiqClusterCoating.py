##
##  LiqCluster_Coating.py
##  LatticePoly
##
##  Created by ppuel on ../../2023.
##  Copyright © 2023 ENS Lyon. All rights reserved.
##

import os
import sys
import numba

import math as m
import numpy as np
import networkx as nx

from vtkReader import vtkReader


class LiqCluster():

    def __init__(self, outputDir, initFrame, threshold=0.24, nMax=10, cutoff=1/2**0.5 + 1e-3):
        self.reader = vtkReader(outputDir, initFrame, readLiq=True, readPoly=True)

        self.nMax = nMax
        self.cutoff = cutoff
        self.threshold = threshold

        self.FirstZone = os.path.join(self.reader.outputDir, "polyFirstZone.res")
        self.SecondZone = os.path.join(self.reader.outputDir, "polySecondZone.res")
        self.ThirdZone = os.path.join(self.reader.outputDir, "polyThirdZone.res")

        if os.path.exists(self.FirstZone) & os.path.exists(self.SecondZone) & os.path.exists(self.ThirdZone):
            print("Files '%s', '%s' and '%s' already exist - aborting" % (self.FirstZone, self.SecondZone, self.ThirdZone))
            sys.exit()


    def Compute(self):
        self.polyFirstZone = np.zeros(self.reader.N,dtype=np.float32)
        self.polySecondZone = np.zeros(self.reader.N,dtype=np.float32)
        self.polyThirdZone = np.zeros(self.reader.N,dtype=np.float32)

        for i in range(self.reader.N):
            self.ProcessFrame(i)

        if (i+1) % 1000 == 0:
            print("Processed %d out of %d configurations" % (i+1, self.reader.N))


    def ProcessFrame(self, i):
        data = next(self.reader)
        nLiq = self.reader.nLiq

        dropPos = data.liqPos[data.liqDens > self.threshold]
        dropDens = data.liqDens[data.liqDens > self.threshold]
        HetPos = data.polyPos[self.reader.polyType == 1]
        connect = self._connectPBC(data.boxDim, dropPos, self.cutoff)

        graph = nx.from_numpy_array(connect)
        clusters = nx.connected_components(graph)
        clusters = sorted(clusters, key=len, reverse=True)

        list_clusters = [np.array(list(c)) for c in sorted(clusters, key=len, reverse=True)]
        for cluster in list_clusters:
            center_of_mass = np.mean(dropPos[cluster, :], 0)
            envelop = dropPos[cluster, :][dropDens[cluster, :] < 0.7]
            envelop_connect = self._connectPBC(data.boxDim, envelop, self.cutoff)
            envelop_graph = nx.from_numpy_array(envelop_connect)
            envelop_clusters = [np.array(list(c)) for c in sorted(nx.connected_components(envelop_graph), key=len, reverse=True)][0]
            for poly in HetPos:
                dist = np.min(np.sum((envelop[envelop_clusters, :] - poly)**2))
                ratio = np.sqrt(dist)/np.sqrt(np.sum((center_of_mass - poly)**2))
                if ratio < 0.6667:
                    self.polyFirstZone[i] += 1
                elif ratio < 1.3333:
                    self.polySecondZone[i] += 1
                elif ratio < 2:
                    self.polyThirdZone[i] += 1


    def Print(self):
        np.savetxt(self.FirstZone, self.polyFirstZone)
        np.savetxt(self.SecondZone, self.polySecondZone)
        np.savetxt(self.ThirdZone, self.polyThirdZone)
        print("\033[1;32mPrinted droplet numbers to '%s'\033[0m" % self.numFile)
        print("\033[1;32mPrinted individual/mean droplet radii to '%s' and '%s'\033[0m" % (self.radFile, self.meanFile))
        print("\033[1;32mPrinted liquid fraction to '%s'\033[0m" % self.fractionFile)
        print("\033[1;32mPrinted poly zone to '%s'\033[0m" % self.FirstZone)


    @staticmethod
    @numba.jit("i4[:,:](f4[:], f4[:,:], f4)", nopython=True)
    def _connectPBC(dims, pts, cutoff):
        nPoints = pts.shape[0]
        connect = np.zeros((nPoints, nPoints), dtype=np.int32)

        for i in range(nPoints):
            for j in range(i+1, nPoints):
                pDist = 0.

                for k in range(3):
                    delta = pts[j, k] - pts[i, k]

                    while abs(delta) > dims[k] / 2.:
                        shift = m.copysign(dims[k], delta)

                        pts[j, k] -= shift
                        delta -= shift

                    pDist += delta**2

                if pDist < cutoff**2:
                    connect[i, j] = 1
                    connect[j, i] = 1

        return connect


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("\033[1;31mUsage is %s outputDir initFrame\033[0m" % sys.argv[0])
        sys.exit()

    outputDir = sys.argv[1]
    initFrame = int(sys.argv[2])

    cluster = LiqCluster(outputDir, initFrame=initFrame)

    cluster.Compute()
    cluster.Print()
