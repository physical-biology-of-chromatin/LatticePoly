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
from iced import normalization

from utils import msdFFT
from vtkReader import vtkReader
import networkx as nx
import time
import json
from numpy import nanmean
import numba
from numba import jit, int32
import warnings
warnings.filterwarnings('ignore')
import pandas as pd


outputDir = sys.argv[1]

origins=np.full(2490,0)
for folder in os.listdir(outputDir):
                if(folder.endswith('.scool')==False and folder.endswith('.gz')==False and folder.endswith('.res')==False):
                        print(folder)
                        for file_name in os.listdir(outputDir+'/'+folder):
                                if file_name=="loaded_mcm.res":
                                        mcm=np.loadtxt(outputDir+'/'+folder+"/loaded_mcm.res")
                                        for i in mcm[:]:
                                                origins[int(i)]=origins[int(i)]+1
print(origins)
                                                
                                                

np.savetxt(outputDir+"/full_loaded_mcm.res",origins)
