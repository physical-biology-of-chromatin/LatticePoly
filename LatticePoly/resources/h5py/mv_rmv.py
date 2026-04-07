import h5py, os, sys, subprocess
from LiqDroplet import Droplet, Event
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import matplotlib as mpl
import pickle

mpl.rcParams['font.size'] = 16

inputDir = "/scratch/Cascade/ppuel/data/"
outputDir = "/home/ppuel/data/"

exp = 52

for file in os.listdir(inputDir):
        if file.split('_')[0] == f"EXP{exp}":
                exp_file = file

exp_path = os.path.join(inputDir, exp_file)

output_path =  os.path.join(outputDir, exp_file)

os.makedirs(exp_path, exist_ok=True)

jll_list = ["0.5", "1.0", "2.0"]
jlp_list = ["0.2", "0.5", "1.0"]
jll_valency_list = [str(i+1) for i in range(6)]
jpl_valency_list = [str(i+1) for i in range(6)]

c_max = 3*3*6*6
c = 0

for jll in jll_list:
        for jlp in jlp_list:
                for jll_valency in jll_valency_list:
                        for jpl_valency in jpl_valency_list:
                                print(f"{c/c_max*100: 3.1f}%", end='\r')
                                c += 1
                                path = os.path.join(exp_path, f"JLL/{jll}/JLP/{jlp}/JLL_VALENCY/{jll_valency}/JPL_VALENCY/{jpl_valency}/DOMAINPATH/")
                                
                                try:
                                        subprocess.run(f"mv {os.path.join(path, 'data/*')} {path}", shell=True, executable="/bin/bash")
                                        subprocess.run(f"rmdir {os.path.join(path, 'data/')}", shell=True, executable="/bin/bash")
                                except:
                                        print(os.listdir(path))