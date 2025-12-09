
import os
import sys
import pandas as pd
import numpy as np
import cooler
from scipy.spatial import cKDTree
from cooler.create import ArrayLoader
from vtkReader_multi import vtkReader

from scipy.spatial.distance import pdist, squareform

def detect_polymer_series(outputDir):
    """
    Detect all unique polymer prefixes from files like:
    Apoly00000.vtk, Bpoly00100.vtk, Something_poly00300.vtk, etc.
    We split on the substring 'poly' and take everything before it.
    """
    polymers = set()

    for f in os.listdir(outputDir):
        if "poly" in f and f.endswith(".vtp"):
            prefix = f.split("poly")[0]
            polymers.add(prefix)

    polymers = sorted(list(polymers))
    print("Found polymers:", polymers)
    return polymers


class MonomerDmap():
	def __init__(self, outputDir, initFrame):
		self.polymers = detect_polymer_series(outputDir)
		self.readers = []
		for P in self.polymers:
			self.readers.append(vtkReader(outputDir, P, initFrame=initFrame, readLiq=False, readPoly=True ))
	
		self.N_chain=np.sum([(self.readers[0].status==-1)+(self.readers[0].status==0)])
		self.contactFile = os.path.join(outputDir,"r_"+str(r)+"_"+str(initFrame)+ "_"+str(interval)+"_SisterC.cool")


		self.Compute(interval)

		self.Print()
			






		
	#NB here is not the finalFrame but the number of iterations
	def Compute(self,finalFrame):
		#self.polyAniso = np.zeros((self.reader.N, self.reader.nDom), dtype=np.float32)
		n_bins=2*self.N_chain
		self.contactProb = np.zeros((n_bins, n_bins), dtype=np.float32)
		
		

		for i in range(0, finalFrame):
			self.ProcessFrame(i)
			if (i+1) % 10 == 0:
				print("Processed %d out of %d configurations" %
					  (i+1, finalFrame))


#self.contactProb=np.rint(self.contactProb/finalFrame)

	def ProcessFrame(self, i):
		for chrom,reader in enumerate(self.readers):
			data = next(reader)
			positions_sister=np.full((2*self.N_chain,3),np.nan)

			for pos_id,pos in enumerate(data.polyPos):
				if(pos_id<self.N_chain):
					positions_sister[pos_id]=pos
				else:
					positions_sister[self.N_chain+data.SisterID[pos_id]]=pos
					
			full_positions=positions_sister
			# Ensure we have a NumPy array
			full_positions = np.asarray(full_positions, dtype=float)
			# Build mask of valid rows (no NaNs)
			mask = ~np.isnan(full_positions).any(axis=1)
			# If there are no valid positions, skip this frame safely
			if not np.any(mask):
				return
			# Positions without NaN → used to build the tree
			clean_positions = full_positions[mask]
			# Build KD-tree only on valid positions
			tree1 = cKDTree(clean_positions)
			# Map tree indices back to original indices
			orig_index = np.nonzero(mask)[0]
			# Query pairs within distance r*0.71
			pairs_clean = tree1.query_pairs(r=r*0.71)
			# Convert tree indices → original indices in the contact matrix
			pairs = {(orig_index[i], orig_index[j]) for (i, j) in pairs_clean}

			# Accumulate contacts symmetrically
			for (i, j) in pairs:
				self.contactProb[i, j] += 1
				self.contactProb[j, i] += 1
		
		
		
			

					



	def Print(self):
		final_bin_size=int(len(self.contactProb))
		final_map = self.contactProb
		
		chrom_names=["S1","S2"]
		chromsizes=pd.Series([self.N_chain*1000,self.N_chain*1000])
		print(chromsizes)
		chromsizes=chromsizes.rename(lambda x: chrom_names[x])
		chromsizes=chromsizes.astype('int64')	
		bins = cooler.binnify(chromsizes, 1000)
		pixels = ArrayLoader(bins, final_map, chunksize=100000000)
		cooler.create_cooler(self.contactFile,bins,pixels)
		
		
		print("\033[1;32mPrinted avg.contact probability to '%s'\033[0m" %self.contactFile)


if __name__ == "__main__":
	if len(sys.argv) != 5:
		print("\033[1;31mUsage is %s outputDir initFrame r interval \033[0m" % sys.argv[0])
		sys.exit()
	
	outputDir = sys.argv[1]
	initFrame = int(sys.argv[2])
	r=float(sys.argv[3])
	interval=int(sys.argv[4])
	
	#init_time=int(sys.argv[4])


monom = MonomerDmap(outputDir, initFrame=initFrame)
