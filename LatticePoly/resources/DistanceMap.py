import sys
import numba

import math as m
import numpy as np
import matplotlib.pyplot as plt

from vtkReader import vtkReader

from matplotlib.colors import LogNorm
from scipy.spatial.distance import squareform


class DistanceMap(vtkReader):

	def __init__(self, vtkDir, frameInit=0, printAllFrames=True):
		vtkReader.__init__(self, vtkDir, frameInit)
		
		self.ReadBox(readPoly=True)

		self.mapFile = self.vtkDir + "/dMap.pdf"
		
		self.frameInit = frameInit
		self.printAllFrames = printAllFrames

		self.cumulDist = 0.
		

	def Compute(self):
		while True:
			try:
				self.ProcessFrame()
			
				if (self.frame-self.frameInit) % 10 == 0:
					print("Processed %d configurations" % (self.frame-self.frameInit))
				
			except IOError:
				break
			
			
	def ProcessFrame(self):
		self.ReadPolyFrame()

		sqDist = self._sqDistPBC(self.boxDims, self.polyPos)
		self.cumulDist += np.sqrt(sqDist)
		
		self.frame += 1
			
		if self.printAllFrames:
			self.Print()
	

	def Print(self):
		if self.frame == self.frameInit:
			print("Could not locate any files to process")
						
		else:
			tadDist = self.cumulDist / (self.frame-self.frameInit)

			tadMap = squareform(tadDist)
			
			tadDomains = np.nonzero(self.polyType)[0]
			tadDomains = np.split(tadDomains, np.where(np.diff(tadDomains) != 1)[0]+1)
		
			fig = plt.figure()

			dMap = plt.imshow(tadMap, norm=LogNorm())
			plt.colorbar(dMap)
			
			for domain in tadDomains:
				x = [domain[0], domain[-1], self.nLoc]
				
				y1 = [domain[0], domain[-1], domain[-1]]
				y2 = [domain[0], domain[0], domain[0]]
		
				plt.fill_between(x=x, y1=y1, y2=y2,  color='red', alpha=0.25)
				plt.fill_between(x=x[:2], y1=y1[:2], color='red', alpha=0.25)

			plt.xlim([0, self.nLoc])
			plt.ylim([0, self.nLoc])

			if self.printAllFrames:
				mapFile = self.vtkDir + "/dMap%04d.png" % (self.frame-1)
				
				plt.savefig(mapFile, format="png", dpi=300)
				print("\033[1;32mPrinted figure to '%s'\033[0m" % mapFile)

			else:
				plt.savefig(self.mapFile, format="pdf", transparent=True)
				print("\033[1;32mPrinted figure to '%s'\033[0m" % self.mapFile)

				plt.show()
				
				
	@staticmethod
	@numba.jit("f4[:](f4[:], f4[:,:])", nopython=True)
	def _sqDistPBC(dims, pts):
		cnt = 0
		n = pts.shape[0]
		
		sqDist = np.zeros(n*(n-1)//2, dtype=np.float32)
		
		for i in range(n-1):
			for j in range(i+1, n):
				for k in range(3):
					delta = pts[j,k] - pts[i,k]
					
					while abs(delta) > dims[k]/2.:
						delta -= m.copysign(dims[k], delta)

					sqDist[cnt] += delta**2
					
				cnt += 1
				
		return sqDist


if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("\033[1;31mUsage is %s vtkDir Neq\033[0m" % sys.argv[0])

		sys.exit()

	vtkDir = sys.argv[1]
	Neq    = int(sys.argv[2])
	
	distMap = DistanceMap(vtkDir, frameInit=Neq, printAllFrames=False)
	
	distMap.Compute()
	distMap.Print()
