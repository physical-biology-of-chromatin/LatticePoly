import os
import sys
import psutil
import math
import numpy as np


outputDir = sys.argv[1]



mrt=[]
for folder in os.listdir(outputDir)[:]:
	current_mrt=[]
	print(folder)
	if(folder.endswith('cool')==False and folder.endswith('.gz')==False and folder.endswith('.res')==False):
		for chrom_n in range(1,17):
			if(np.sum(np.array(os.listdir(outputDir+'/'+folder))=="chrom_"+str(chrom_n-1)+"_timing.res")):
				if 0==0:
					file_path = outputDir+'/'+folder+"/chrom_"+str(chrom_n-1)+"_timing.res"
					if  os.path.getsize(file_path) != 0:
						current_mrt.extend(list(np.loadtxt(file_path)[:]))
						break
		mrt.append(current_mrt)
				
print(mrt)                
nsim=len(mrt)
n=6
dp = np.array(np.arange(0, 1+1/(2*n), 1/n)*100, dtype=int)
MRTp = np.zeros_like(mrt[0])

for MRT in mrt:
	for mon in range(len(MRT)):
		if(MRT[mon]==-1.0):
			MRT[mon]=np.nanmax(MRT)
	percentils = np.percentile(MRT, dp)
	for ip, (p1, p2) in enumerate(zip(percentils[:-1], percentils[1:]), 1):
		MRTp[(MRT > p1) & (MRT <= p2)] += ip/n
	MRTp[MRT == 0] += 1/n
MRT_normed = (MRTp/nsim - 1/(2*n))
        
np.savetxt(outputDir+"/mrt.res",1-MRT_normed)
