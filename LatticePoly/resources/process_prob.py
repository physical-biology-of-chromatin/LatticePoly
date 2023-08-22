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
allint=[]
for folder in os.listdir(outputDir):
	for file_name in os.listdir(outputDir+'/'+folder):
		if file_name.endswith('01.res'):
			file_path = os.path.join(outputDir+'/'+folder, file_name)
			cis.append(np.loadtxt(file_path))


finalcis=nanmean(cis,axis=0)
np.savetxt(outputDir+"finalprob.res", finalcis)
