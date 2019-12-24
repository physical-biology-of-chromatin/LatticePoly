import os
import sys
import numba

import numpy as np

from vtk import vtkXMLPolyDataReader
from vtk.util import numpy_support as vn


class vtkReader():

	def __init__(self, outputDir, frame=0):
		outputDir = outputDir.strip("/")
		
		self.liqFile  = outputDir + "/liq%04d.vtp"
		self.polyFile = outputDir + "/poly%04d.vtp"

		self.reader = vtkXMLPolyDataReader()

		self.frame = frame
		self.outputDir = outputDir
					
		self.boxDims = 0

			
	def ReadBox(self, readPoly=False, readLiq=False):
		try:
			boxFile = "%s/box.vtp" % self.outputDir
			
			self._read(boxFile)

			boxData = self.reader.GetOutput()
			vertPos = vn.vtk_to_numpy(boxData.GetPoints().GetData())

			self.boxDims = vertPos.max(axis=0)
			
			print("Box linear dimensions: (%.0f,%.0f,%.0f)" % tuple(self.boxDims))

			if readPoly:
				self.ReadPolyFrame()
				
				self.nHom = np.count_nonzero(self.polyType == 0)
				self.nHet = np.count_nonzero(self.polyType == 1)
				
				self.nLoc = self.nHom+self.nHet
								
				print("Found %d TADs inc. %d heterochromatic loci" % (self.nLoc, self.nHet))
				
			if readLiq:
				self.ReadLiqFrame()
				
				self.nProt = np.count_nonzero(self.liqType == 0)
				self.nFrap = np.count_nonzero(self.liqType == 1)
				
				self.nLiq = self.nProt+self.nFrap
				
				print("Initial liquid state: %d occupied sites (%d FRAP'ed)" % (self.nLiq, self.nFrap))
			
		except IOError as err:
			print("%s - aborting" % err)
			
			sys.exit()


	def ReadLiqFrame(self):
		self._read(self.liqFile % self.frame)
		
		liqData = self.reader.GetOutput()
		
		self.liqPos = vn.vtk_to_numpy(liqData.GetPoints().GetData())
		self.liqType = vn.vtk_to_numpy(liqData.GetPointData().GetArray("FRAP type"))
		
		self.liqDisp = vn.vtk_to_numpy(liqData.GetPointData().GetArray("Displacement"))
		
		
	def ReadPolyFrame(self, backInBox=False):
		self._read(self.polyFile % self.frame)
		
		polyData = self.reader.GetOutput()
		
		self.polyPos = vn.vtk_to_numpy(polyData.GetPoints().GetData())
		self.polyType = vn.vtk_to_numpy(polyData.GetPointData().GetArray("TAD type"))
		
		if backInBox:
			self._backInBox(self.boxDims, self.polyPos)


	def _read(self, file):
		if not os.path.exists(file):
			raise IOError("Could not find file '%s'" % file)
		
		self.reader.SetFileName(file)
		self.reader.Update()


	@staticmethod
	@numba.jit("void(f4[:], f4[:,:])", nopython=True)
	def _backInBox(dims, pts):
		n = pts.shape[0]
				
		for i in range(n):
			for j in range(3):
				while pts[i,j] < 0: pts[i,j] += dims[j]
				while pts[i,j] >= dims[j]: pts[i,j] -= dims[j]
