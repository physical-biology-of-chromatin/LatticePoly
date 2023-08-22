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
import matplotlib.pyplot as plt
from utils import msdFFT
from vtkReader import vtkReader
import networkx as nx
import time
import json
from numpy import nanmean
import numba
from numba import jit, int32
from cooler.create import ArrayLoader
import cooler
import warnings
warnings.filterwarnings('ignore')
import pandas as pd

import os
import numpy as np
import matplotlib.pyplot as plt
outputDir="./LatticePoly/LatticePoly/data/output/September22/cohesion/keco_1.0/"

#Define all matrices names
for folder in os.listdir(outputDir):
	if folder.endswith('.gz')==False and folder.endswith('.res')==False and folder.endswith('scool')==False :
		for file_name in os.listdir(outputDir+'/'+folder):
			if file_name.endswith('cis_hic.res')==True and file_name.startswith('r_1')==True:
				print(file_name)
				mtx=np.loadtxt(outputDir+'/'+folder+"/"+file_name)
				plt.matshow(np.log2(mtx))
				plt.savefig(outputDir+'/'+folder+"/"folder+file_name+'.png')

				
				
