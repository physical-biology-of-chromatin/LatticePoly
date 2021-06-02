##
##  vtkReader.py
##  LatticePoly
##
##  Created by mtortora on 12/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import os
import numba

import numpy as np

from vtk import vtkXMLPolyDataReader
from vtk.util import numpy_support as vn

from fileseq import findSequenceOnDisk
from fileseq.exceptions import FileSeqException


class vtkReader():

	def __init__(self, outputDir, initFrame=-1, readLiq=False, readPoly=True, backInBox=False):
		self.outputDir = os.path.realpath(outputDir)

		if os.path.exists(self.outputDir):
			self._liqFile  = os.path.join(self.outputDir, "liq%05d.vtp")
			self._polyFile = os.path.join(self.outputDir, "poly%05d.vtp")
		
			print("\033[1;34mParsing directory '%s'\033[0m" % self.outputDir)
			
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

		self._reader = vtkXMLPolyDataReader()
		
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
			
			
	def InitReader(self, initFrame):
		try:
			boxFile = os.path.join(self.outputDir, "box.vtp")
			boxData = self._read(boxFile)
			
			self._checkRange(initFrame)

			vertPos = vn.vtk_to_numpy(boxData.GetPoints().GetData())
			self.boxDim = vertPos.max(axis=0)
			
			print("Box linear dimensions: (%.0f,%.0f,%.0f)" % tuple(self.boxDim))

			if self._readLiq:
				self._readLiqFrame()
				
				self.nLiq = self.liqDens.size
				
				print("Initial liquid state: %d occupied sites" % self.nLiq)
				
			if self._readPoly:
				self._readPolyFrame()
				
				self.nTad = self.polyType.size
								
				print("Initial chromatin state: %d TADs inc. %d heterochromatic loci" % (self.nTad, self.nHet))
			
		except IOError:
			raise


	def _readLiqFrame(self):
		liqData = self._read(self._liqFile % self.frame)
		
		self.liqPos = vn.vtk_to_numpy(liqData.GetPoints().GetData())
		
		self.liqDens = vn.vtk_to_numpy(liqData.GetPointData().GetArray("Density"))
		self.liqDisp = vn.vtk_to_numpy(liqData.GetPointData().GetArray("Displacement"))
		
		
	def _readPolyFrame(self):
		polyData = self._read(self._polyFile % self.frame)
		
		self.polyPos = vn.vtk_to_numpy(polyData.GetPoints().GetData())
		
		self.polyType = vn.vtk_to_numpy(polyData.GetPointData().GetArray("TAD type"))
		self.polyPainter = vn.vtk_to_numpy(polyData.GetPointData().GetArray("Painter status"))

		self.nEuc = np.count_nonzero(self.polyType == 0)
		self.nHet = np.count_nonzero(self.polyType == 1)
		
		hetDomains = np.nonzero(self.polyType == 1)[0]
		self.domains = np.split(hetDomains, np.where(np.diff(hetDomains) != 1)[0] + 1)
		
		self.nDom = len(self.domains)
		
		if self._backInBox:
			self._fixPBCs(self.boxDim, self.polyPos)


	def _read(self, file):
		if not os.path.exists(file):
			raise IOError("Could not find file '%s'" % file)
			
		self._reader.SetFileName(file)
		self._reader.Update()

		return self._reader.GetOutput()
			
			
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
				polySeq = findSequenceOnDisk(self.outputDir + '/poly@.vtp')
		
				self._minFramePoly = polySeq.start()
				self._maxFramePoly = polySeq.end()
		
			except FileSeqException:
				raise IOError("Could not locate any polymer configuration files in '%s'" % self.outputDir)
		
		if self._readLiq:
			try:
				liqSeq = findSequenceOnDisk(self.outputDir + '/liq@.vtp')
		
				self._minFrameLiq = liqSeq.start()
				self._maxFrameLiq = liqSeq.end()
			
			except FileSeqException:
				raise IOError("Could not locate any liquid configuration files in '%s'" % self.outputDir)


	@staticmethod
	@numba.jit("void(f4[:], f4[:,:])", nopython=True)
	def _fixPBCs(dims, pts):
		nPoints = pts.shape[0]
		
		for i in range(nPoints):
			for j in range(3):
				while pts[i, j] < 0:
					pts[i, j] += dims[j]
				
				while pts[i, j] >= dims[j]:
					pts[i, j] -= dims[j]
