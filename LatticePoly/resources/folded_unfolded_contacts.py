# -*- coding: utf-8 -*-### MonomerDist.py# LatticePoly### Based on mtortora's code.# Copyright © 2021 ENS Lyon. All rights reserved.##import osimport sysimport numpy as npfrom scipy.spatial import cKDTreefrom vtkReader import vtkReaderfrom scipy.spatial.distance import pdist, squareformclass MonomerDmap():    def __init__(self, outputDir, initFrame):        self.reader = vtkReader(outputDir, initFrame,readLiq=False, readPoly=True)        self.reader_folded = vtkReader(outputDir, initFrame,readLiq=False, readPoly=True, backInBox=True)        #self.anisoFile = os.path.join(self.reader.outputDir, "polyAniso.res")        self.ratioFile = os.path.join(self.reader.outputDir, "intra_inter_collisions.res")    def Compute(self):        self.intra_collisions = []        self.inter_collisions = []                for i in range(self.reader.N):            self.ProcessFrame(i)    def ProcessFrame(self, i):        data = next(self.reader)        tree1   = cKDTree(data.polyPos[:], boxsize = None)        pairs = tree1.query_pairs(r = 2*0.71) # NN distance FCC lattice 1/np.sqrt(2) = 0.71                data = next(self.reader_folded)        self.intra_collisions.append(len(pairs))        tree1   = cKDTree(data.polyPos[:], boxsize = None)        pairs = tree1.query_pairs(r = 2*0.71) # NN distance FCC lattice 1/np.sqrt(2) = 0.71        self.inter_collisions.append(len(pairs))    def Print(self):        np.savetxt(self.ratioFile, 1/(np.array(self.inter_collisions)/np.array(self.intra_collisions)) )if __name__ == "__main__":    if len(sys.argv) !=3:        print("\033[1;31mUsage is %s outputDir initFrame\033[0m" % sys.argv[0])        sys.exit()    outputDir = sys.argv[1]    initFrame = int(sys.argv[2])    monom = MonomerDmap(outputDir, initFrame=initFrame)    monom.Compute()    monom.Print()