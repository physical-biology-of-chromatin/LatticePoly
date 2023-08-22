##
##  PolyMSD.py
##  LatticePoly
##
##  Created by mtortora on 15/12/2019.
##  Copyright © 2019 ENS Lyon. All rights reserved.
##

import os
import sys
import psutil
import math
import numpy as np

from utils import msdFFT
from vtkReader import vtkReader
import networkx as nx
import time
import json
from numpy import nanmean



outputDir = sys.argv[1]

cis=[]
trans=[]
time=[]
for folder in os.listdir(outputDir):
	for file_name in os.listdir(outputDir+'/'+folder):
		if file_name.endswith('log.out'):
			file_path = os.path.join(outputDir+'/'+folder, file_name)
			text=open(file_path).readlines()
			i=0
			while(text[-1][15:15+i]!= " "):
				i+=1
			time.append(float(text[-1][15:14+i]))


			


np.savetxt(outputDir+"time.res", time)
