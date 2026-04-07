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

import utils

import h5py

import networkx as nx

import pickle

from hdf5Reader import hdf5Reader

def node_to_tuple(node):
    return(node.replace('(', "").replace(')', "").split(","))

class LiqSimDroplet():

    def __init__(self, experience, reference):

        self.exp = experience
        self.XnfsDir = "/Xnfs/physbiochrom/ppuel/data/"
        
        if int(reference) > 0:
            self.is_ref = True
            self.ref = reference
        else :
            self.is_ref = False

        self.expName = utils.find_exp(self.exp, self.XnfsDir)

        self.expDir = os.path.join(self.XnfsDir, self.expName)

        dict_parameters, _, metaParameterN, metaParameterNmeas = utils.exp_mapping(self.expDir)

        self.metaParameterN = int(metaParameterN)
        self.metaParameterNmeas = int(metaParameterNmeas)

        self.pathList = utils.exp_pathList(dict_parameters, -1, self.expDir)
        
        if self.is_ref:
            self.refName = utils.find_exp(self.ref, self.XnfsDir)

            self.refDir = os.path.join(self.XnfsDir, self.refName)

            for file in os.listdir(f"/home/ppuel/data/{self.ref}"):
                if "r_max_3D" in file:
                    param = file.replace('.npy', '').split("_")[3:]
                    self.list_param_ref = param[::2]
                    break

        else :
            
            self.N_param_ref = len(dict_parameters)



                


    def Compute(self):
        for i in range(len(self.pathList)):
            self.ProcessPath(self.pathList[i])

            if (i + 1) % 1 == 0:
                print(f"Processed {i/len(self.pathList)*100:0.1f}% of configurations",end='\r')


    def ProcessPath(self, path):

        liqFraction = np.zeros(self.metaParameterNmeas)
        liqDropNum = np.zeros(self.metaParameterNmeas)
        liqDropSize = np.zeros(self.metaParameterNmeas)

        r_max = 7

        Nb_traj = self.metaParameterN

        reader = hdf5Reader(os.path.join(path, "N/0"), "traj.h5", readLiq=True, readPoly=False)

        nLiq = reader.nLiq

        reader.Close()
            
        for n in range(self.metaParameterN):
            file = open(os.path.join(path,f"N/{n}/liq_droplets.pickle"), "rb")
            droplet_dict = pickle.load(file)
            file.close()
            for droplet in droplet_dict.values():
                r_droplet = ((droplet.tau-1)**2+(np.max(droplet.sizes)-2)**2)**(1/2)

                if r_droplet > r_max:
                    frame = np.fromiter(map(lambda x: x[0], droplet.nodes), dtype=np.int32)
                    
                    liqFraction[frame] += droplet.sizes
                    liqDropNum[frame] += 1
                    
                    # for t in range(droplet.tau):
                        # frame = droplet.nodes[t][0]
                        # size = droplet.sizes[t]
                        
                        # liqFraction[frame] += droplet.sizes
                        # liqDropNum[frame] += 1
            # print(f"{n/self.metaParameterN*100:0.1f}%",end = '\r')
                
            # except : 
            #     Nb_traj -= 1
            #     print(os.path.join(path,f"N/{n}/liq_droplets.pickle"))

        if Nb_traj != 0:
            liqDropSize = np.where(liqDropNum > 0, liqFraction/liqDropNum, np.nan)/Nb_traj
            liqFraction /= (Nb_traj * nLiq)
            liqDropNum /= Nb_traj
        
        else:
            print(path)

        with h5py.File(os.path.join(path, "aggregated_process.h5"), 'a') as processFile:
                        
            utils.PrintDataset(processFile, "liqFraction", data = liqFraction)
            utils.PrintDataset(processFile, "liqDropNum", data = liqDropNum)
            utils.PrintDataset(processFile, "liqDropSize", data = liqDropSize)
                
        

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("\033[1;31mUsage is %s experience reference\033[0m" % sys.argv[0])
        sys.exit()

    experience = sys.argv[1]
    reference = sys.argv[2]

    SimDroplet = LiqSimDroplet(experience, reference)

    SimDroplet.Compute()
        
    print("\n")
    print("LiqSimDroplet : Done\n\n")
