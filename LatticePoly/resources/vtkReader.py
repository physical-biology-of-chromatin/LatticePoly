##
##  vtkReader.py
##  LatticePoly
##
##  Created by mtortora on 12/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import os
import sys
import numba

import numpy as np

from vtk import vtkXMLPolyDataReader
from vtk.util import numpy_support as vn

from fileseq import findSequenceOnDisk
from fileseq.exceptions import FileSeqException


class vtkReader():

	def __init__(self, outputDir, initFrame=0, readLiq=False, readPoly=True):
		self.outputDir = outputDir.strip("/")

		if os.path.exists(self.outputDir):
			print("\033[1;34mParsing directory '%s'...\033[0m" % self.outputDir)
			
		else:
			print("\033[1;31mDirectory '%s' does not exist\033[0m" % self.outputDir)
			sys.exit()
			
		self.liqFile  = os.path.join(self.outputDir, "liq%05d.vtp")
		self.polyFile = os.path.join(self.outputDir, "poly%05d.vtp")
		
		self.InitReader(initFrame, readLiq, readPoly)
			
			
	def InitReader(self, initFrame, readLiq, readPoly):
		self.reader = vtkXMLPolyDataReader()
		self.frame = self.initFrame = initFrame
		
		try:
			boxFile = os.path.join(self.outputDir, "box.vtp")
			boxData = self._read(boxFile)
			
			self._checkRange(readPoly, readLiq)

			vertPos = vn.vtk_to_numpy(boxData.GetPoints().GetData())
			self.boxDim = vertPos.max(axis=0)
			
			print("Box linear dimensions: (%.0f,%.0f,%.0f)" % tuple(self.boxDim))

			if readLiq:
				self.ReadLiqFrame()
				
				self.nLiq = self.liqDens.size
				
				print("Initial liquid state: %d occupied sites" % self.nLiq)
				
			if readPoly:
				self.ReadPolyFrame()
				
				self.nHom = np.count_nonzero(self.polyType == 0)
				self.nHet = np.count_nonzero(self.polyType == 1)
				
				self.nLoc = self.nHom + self.nHet
								
				print("Found %d TADs inc. %d heterochromatic loci" % (self.nLoc, self.nHet))
			
		except IOError as err:
			print("%s - aborting" % err)
			sys.exit()


	def ReadLiqFrame(self, readAttr=True):
		liqData = self._read(self.liqFile % self.frame)
		
		self.liqPos = vn.vtk_to_numpy(liqData.GetPoints().GetData())
		self.liqDisp = vn.vtk_to_numpy(liqData.GetPointData().GetArray("Displacement"))

		if readAttr:
			self.liqDens = vn.vtk_to_numpy(liqData.GetPointData().GetArray("Density"))
		
		
	def ReadPolyFrame(self, readAttr=True, backInBox=False):
		polyData = self._read(self.polyFile % self.frame)
		
		self.polyPos = vn.vtk_to_numpy(polyData.GetPoints().GetData())
		
		if readAttr:
			self.polyType = vn.vtk_to_numpy(polyData.GetPointData().GetArray("TAD type"))
		
		if backInBox:
			self._backInBox(self.boxDim, self.polyPos)


	def _read(self, file):
		if not os.path.exists(file):
			raise IOError("Could not find file '%s'" % file)
			
		self.reader.SetFileName(file)
		self.reader.Update()

		return self.reader.GetOutput()
			
			
	def _checkRange(self, readPoly, readLiq):
		self._parseFileSeqs(readPoly, readLiq)
	
		if (readPoly & readLiq):
			minFrame = max(self._minFrameLiq, self._minFramePoly)
			maxFrame = min(self._maxFrameLiq, self._maxFramePoly)
				
			if not minFrame <= self.initFrame <= maxFrame:
				raise IOError("Frame not in range (%d, %d)" % (minFrame, maxFrame))
					
			self.N = maxFrame - self.initFrame + 1
			
		elif readPoly:
			if not self._minFramePoly <= self.initFrame <= self._maxFramePoly:
				raise IOError("Frame not in range (%d, %d)" % (self._minFramePoly, self._maxFramePoly))
					
			self.N = self._maxFramePoly - self.initFrame + 1
				
		elif readLiq:
			if not self._minFrameLiq <= self.initFrame <= self._maxFrameLiq:
				raise IOError("Frame not in range (%d, %d)" % (self._minFrameLiq, self._maxFrameLiq))
				
			self.N = self._maxFrameLiq - self.initFrame + 1
			
			
	def _parseFileSeqs(self, readPoly, readLiq):
		if readPoly:
			try:
				polySeq = findSequenceOnDisk(self.outputDir + '/poly@.vtp')
		
				self._minFramePoly = polySeq.start()
				self._maxFramePoly = polySeq.end()
		
			except FileSeqException:
				raise IOError("Could not locate any polymer configuration files in '%s'" % self.outputDir)
		
		if readLiq:
			try:
				liqSeq = findSequenceOnDisk(self.outputDir + '/liq@.vtp')
		
				self._minFrameLiq = liqSeq.start()
				self._maxFrameLiq = liqSeq.end()
			
			except FileSeqException:
				raise IOError("Could not locate any liquid configuration files in '%s'" % self.outputDir)


	@staticmethod
	@numba.jit("void(f4[:], f4[:,:])", nopython=True)
	def _backInBox(dims, pts):
		n = pts.shape[0]
		
		for i in range(n):
			for j in range(3):
				while pts[i,j] < 0:
					pts[i,j] += dims[j]
				
				while pts[i,j] >= dims[j]:
					pts[i,j] -= dims[j]
