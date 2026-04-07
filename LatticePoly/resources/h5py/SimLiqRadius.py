##
##  LiqPolyCoincidence.py
##  LatticePoly
##Liq
##  Created by ppuel on ../../2023.
##  Copyright © 2023 ENS Lyon. All rights reserved.
##

import os
import sys

from LiqDroplet import LifeTime, Droplet, Event

import numpy as np

import h5py

import networkx as nx

import pickle

import time

import subprocess


def node_to_tuple(node):
    return(node.replace('(', "").replace(')', "").split(","))

class SimLiqRadius():

    def __init__(self, inputDir, outputDir, N):

        self.inputDir = inputDir
        self.outputDir = outputDir
        self.N = N

 
    def Process(self):

        r_max = -1
      
        for n in range(self.N):
            with open(os.path.join(self.inputDir,f"N/{n}/liq_droplets.pickle"), "rb") as pfile:
                droplet_dict = pickle.load(pfile)
            
            for droplet in droplet_dict.values():
                r_droplet = ((droplet.tau-1)**2+(np.max(droplet.sizes)-2)**2)**(1/2)

                if r_droplet > r_max:
                    r_max = r_droplet
        
        with open(os.path.join(self.outputDir,"r_max.txt"),'w') as tfile:
            tfile.write(f"{r_max:0.4f}")


if __name__ == "__main__":
    if len(sys.argv) != 9:
        print("\033[1;31mUsage is %s inputDir refDir outputDir N Nmeas L is_poly\033[0m" % sys.argv[0])
        sys.exit()

    inputDir = sys.argv[1]

    outputDir = sys.argv[3]
    N = int(sys.argv[4])
    Nmeas = int(sys.argv[5])
    L = int(sys.argv[6])
    lDens = float(sys.argv[7])
    is_poly = sys.argv[8]=="1"

    SimRadius = SimLiqRadius(inputDir, outputDir, N)

    SimRadius.Process()
    
    print("\n")
    print("SimLiqRadius : Done\n\n")