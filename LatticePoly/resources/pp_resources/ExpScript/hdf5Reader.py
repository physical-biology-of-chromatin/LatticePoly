##
##  hdf5Reader.py
##  LatticePoly
##
##  Created by ppuel on 18/10/2024.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import os
# import numba

import numpy as np

import h5py


class hdf5Reader():

        def __init__(self, outputDir, fileName, initFrame=-1, readLiq=False, readPoly=True, backInBox=False):
                self.outputDir = os.path.realpath(outputDir)
                self.fileName = fileName
                self.filePath = os.path.join(self.outputDir,self.fileName)

                if os.path.exists(self.outputDir):
                        self._h5File  = os.path.join(self.outputDir, "test.h5")
                        # print("\033[1;34mParsing directory '%s'\033[0m" % self.outputDir)
                        
                else:
                        raise IOError("Directory '%s' does not exist" % self.outputDir)
                
                self.liqPos = None
                self.polyPos = None

                self.liqDisp = None
                self.liqDens = None
                
                self.polyType = None
                self.polyPainter = None
                
                self.boxDim = None
                
                self.N = 0
                self.frame = self.initFrame = initFrame
                self.nLiq = self.nTad = self.nEuc = self.nHet = 0
                
                self._readLiq = readLiq
                self._readPoly = readPoly
                self._backInBox = backInBox

                self.file = h5py.File(self.filePath,'r')
                
                self.InitReader(initFrame)

        
        def __iter__(self):
                return self
        
        
        def __next__(self):
                if self.frame < self.N + self.initFrame:
                        if self._readLiq:
                                self._readLiqFrame()
                                
                        if self._readPoly:
                                self._readPolyFrame()
                                
                        self.frame += 1

                        return self
                
                else:
                        raise StopIteration
                        
        def Close(self):
                self.file.close()

        def InitReader(self, initFrame):
                try:                        
                        self._checkRange(initFrame)

                        self.boxDim = [self.file["Pol"].attrs['L'][0] for i in range(3)] #/!\
                        
                        # print("Box linear dimensions: (%.0f,%.0f,%.0f)" % tuple(self.boxDim))

                        if self._readLiq:
                                self._readLiqFrame()
                                
                                self.nLiq = self.liqDens.size
                                
                                # print("Initial liquid state: %d occupied sites" % self.nLiq)
                                
                        if self._readPoly:
                                self._readPolyFrame()
                                
                                self.nTad = self.polyType.size
                                                                
                                # print("Initial chromatin state: %d TADs inc. %d heterochromatic loci" % (self.nTad, self.nHet))
                        
                except IOError:
                        raise


        def _readLiqFrame(self):                
                self.liqPos = np.array(self.file['Liq']['{:05d}_position_dataset'.format(self.frame)])
                
                self.liqDens = np.array(self.file['Liq']['{:05d}_density_dataset'.format(self.frame)])
                self.liqDisp = np.array(self.file['Liq']['{:05d}_displacement_dataset'.format(self.frame)])
                
                
        def _readPolyFrame(self):
               
                self.polyPos = np.array(self.file['Pol']['{:05d}_position_dataset'.format(self.frame)])
                
                self.polyType = np.array(self.file['Pol']['{:05d}_type_dataset'.format(self.frame)])
                self.polyPainter = np.array(self.file['Pol']['{:05d}_painter_dataset'.format(self.frame)])

                self.nEuc = np.count_nonzero(self.polyType == 0)
                self.nHet = np.count_nonzero(self.polyType == 1)

                hetDomains = np.nonzero(self.polyType == 1)[0]
                self.domains = np.split(hetDomains, np.where(np.diff(hetDomains) != 1)[0] + 1)
                
                self.nDom = len(self.domains)
                
                if self._backInBox:
                        self._fixPBCs(self.boxDim, self.polyPos)

                        
                        
        def _checkRange(self, initFrame):
                self._parseFileSeqs()
        
                if (self._readPoly & self._readLiq):
                        minFrame = max(self._minFrameLiq, self._minFramePoly)
                        maxFrame = min(self._maxFrameLiq, self._maxFramePoly)
                                
                        self.frame = self.initFrame = minFrame if initFrame == -1 else initFrame
                                
                        if not minFrame <= self.initFrame <= maxFrame:
                                raise IOError("Frame not in range (%d, %d)" % (minFrame, maxFrame))
                                        
                        self.N = maxFrame - self.initFrame + 1
                        
                elif self._readPoly:
                        self.frame = self.initFrame = self._minFramePoly if initFrame == -1 else initFrame
                
                        if not self._minFramePoly <= self.initFrame <= self._maxFramePoly:
                                raise IOError("Frame not in range (%d, %d)" % (self._minFramePoly, self._maxFramePoly))
                                        
                        self.N = self._maxFramePoly - self.initFrame + 1
                                
                elif self._readLiq:
                        self.frame = self.initFrame = self._minFrameLiq if initFrame == -1 else initFrame

                        if not self._minFrameLiq <= self.initFrame <= self._maxFrameLiq:
                                raise IOError("Frame not in range (%d, %d)" % (self._minFrameLiq, self._maxFrameLiq))
                                
                        self.N = self._maxFrameLiq - self.initFrame + 1
                        
                        
        def _parseFileSeqs(self):
                if self._readPoly:
                        try:
                                polySeq = list(self.file["Pol"].keys())
                                polySeq.sort()

                                self._minFramePoly = int(polySeq[0].split('_')[0])
                                self._maxFramePoly = int(polySeq[-1].split('_')[0])

                        except:
                                raise IOError("Could not locate any polymer configuration dataset in '%s'" % self.filePath)
                if self._readLiq:
                        try:
                                liqSeq = list(self.file["Liq"].keys())
                                liqSeq.sort()

                                self._minFrameLiq = int(liqSeq[0].split('_')[0])
                                self._maxFrameLiq = int(liqSeq[-1].split('_')[0])

                        except:
                                raise IOError("Could not locate any liquid configuration dataset in '%s'" % self.filePath)

      
        # @staticmethod       
        # @numba.jit("void(f4[:], f4[:,:])", nopython=True)
        def _fixPBCs(self, dims, pts):
                nPoints = pts.shape[0]
                
                for i in range(nPoints):
                        for j in range(3):
                                while pts[i, j] < 0:
                                        pts[i, j] += dims[j]
                                
                                while pts[i, j] >= dims[j]:
                                        pts[i, j] -= dims[j]

