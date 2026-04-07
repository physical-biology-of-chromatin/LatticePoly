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

class LiqSimDroplet():

    def __init__(self, experience): #, reference):

        self.exp = experience
        # self.ref = reference
        self.outputDir = "/Xnfs/physbiochrom/ppuel/data/"
        self.refDir = "/Xnfs/physbiochrom/ppuel/data/"
        self.scratchDir = "/scratch/Cascade/ppuel/data/"

        self.outputPathList = []

        self.N = 0
        self.Nframe = 0

        # get experience name

        verif = False

        for tmp_experience in os.listdir(self.outputDir):
            if tmp_experience.split("_")[0] == f'EXP{experience}':
                self.exp = tmp_experience
                self.outputDir = os.path.join("/Xnfs/physbiochrom/ppuel/data/", self.exp)
                verif = not(verif)

        if not(verif):
            print(f"The experience {experience} is not found or found multiple time")
            sys.exit()

        print(f"exp = {experience}")
        # verif = False

        # for tmp_reference in os.listdir(self.refDir):
        #     if tmp_reference.split("_")[0] == f'EXP{reference}':
        #         self.ref = tmp_reference
        #         verif = not(verif)

        # if not(verif):
        #     print(f"The reference {reference} is not found or found multiple time")
        #     sys.exit()

        # Construct the arborescence of the experience, and gather metadata

        tmpList = list(map(lambda x: os.path.join(self.outputDir, x), os.listdir(self.outputDir)))
        bool_metadata = 0
        a = 1e5
        c = 0
        while len(tmpList) > 0 and c < a:
            c += 1
            tmpPath = tmpList[0]

            print(f"{c: <30}", end='\r')
            if os.path.isfile(tmpPath) or os.listdir(tmpPath) == []:
                tmpList.remove(tmpPath)
            elif "N" in os.listdir(tmpPath):
                self.outputPathList.append(tmpPath)
                tmpList.remove(tmpPath)
                if not(bool_metadata):
                    print("\nStart gathering meta_data")
                    self.N = len(os.listdir(os.path.join(tmpPath,'N')))-40
                    # tmpH5File = h5py.File(os.path.join(tmpPath,'N/0/process.h5'),'r')
                    self.Nframe = 101
                    bool_metadata = 1
            else:
                for tmpFile in os.listdir(tmpPath):
                    if os.path.isdir(os.path.join(tmpPath,tmpFile)):
                        tmpList.append(os.path.join(tmpPath,tmpFile))
                tmpList.remove(tmpPath)


        self.list_param = self.outputPathList[0].split('/')[6:][::2]
        print(self.list_param)
        
        self.Nparam = len(self.list_param)

        print('sort')

        if 'LDENS' in self.list_param:
            self.outputPathList.sort(key=lambda x: float(x.split('LDENS/')[1].split("/")[0]))
            
        print(f"\ninit_end")



                


    def Compute(self):
        print("compute_start")
        for i in range(len(self.outputPathList)):
            subprocess.run("clear", shell=True, executable="/bin/bash")
            print(f"{i/len(self.outputPathList)*100:0.1f}%")
            
            outputPath = self.outputPathList[i]
            self.ProcessPath(outputPath)



    def ProcessPath(self, outputPath):


        scratchPath = os.path.join(self.scratchDir,'/'.join(outputPath.split("/")[5:]))
        start = time.time()
        
        Nb_traj = self.N

        liqFraction = np.zeros(self.Nframe)
        liqDropNum = np.zeros(self.Nframe)
        liqDropSize = np.zeros(self.Nframe)

        # path_r_max = f"/home/ppuel/data/{self.ref}/r_max_3D"
        # for param in range(self.Nparam):
        #     value = outputPath.split('/')[outputPath.split('/').index(self.list_param[param].upper())+1]
        #     path_r_max += f"_{self.list_param[param]}_{value}"

        # path_r_max += '.npy'

        r_max = 5 # np.load(path_r_max)

        if 'LDENS' in self.list_param:
            #     # tmpH5File = h5py.File(os.path.join(outputPath,'N/0/process.h5'),'r')
            self.nLiq = int(float(outputPath.split('/')[outputPath.split('/').index('LDENS')+1])*24**3*4)
        #     # _, self.nLiq, _ = np.shape(tmpH5File["liqInfo"])
        else:
            self.nLiq = int(0.019*48**3*4)

        print(self.nLiq)

        print(f"Init : {start - time.time():0.2f}s")
        start = time.time()

        for n in range(self.N):
            print(f"{n/self.N*100:0.1f}%", end = "\r")
            try : 
                file = open(os.path.join(outputPath,f"N/{n}/liq_droplets.pickle"), "rb")
                droplet_dict = pickle.load(file)

                file.close()
                for droplet in droplet_dict.values():
                    r_droplet = ((droplet.tau-1)**2+(np.max(droplet.sizes)-2)**2)**(1/2)

                    if r_droplet > r_max:
                        
                        for t in range(droplet.tau):
                            frame = droplet.nodes[t][0]
                            size = droplet.sizes[t]
                            
                            liqFraction[frame] += size
                            liqDropNum[frame] += 1
            
                liqDropSize += np.divide(liqFraction, liqDropNum, out=np.zeros(self.Nframe), where=liqDropNum>.5)
        
            except : 
                Nb_traj -= 1
                print(os.path.join(outputPath,f"N/{n}/liq_droplets.pickle"))

        if Nb_traj != 0:
                
            liqFraction /= (Nb_traj * self.nLiq)
            liqDropNum /= Nb_traj
            liqDropSize /= Nb_traj
        
        else:
            print(outputPath)

        print(f"Compute : {start - time.time():0.2f}s")
        start = time.time()

        with h5py.File(os.path.join(scratchPath, "post_process.h5"), "a") as file:

            if 'liqFraction' in file.keys():
                tmp = file['liqFraction']
                tmp[:] = liqFraction
            else:
                file.create_dataset("liqFraction", data = liqFraction)
            
            if 'liqDropNum' in file.keys():
                tmp = file['liqDropNum']
                tmp[:] = liqDropNum
            else:
                file.create_dataset("liqDropNum", data = liqDropNum)
            
            if 'liqDropSize' in file.keys():
                tmp = file['liqDropSize']
                tmp[:] = liqDropSize
            else:
                file.create_dataset("liqDropSize", data = liqDropSize)
            

        print(f"Print : {start - time.time():0.2f}s")
        start = time.time()
                    


if __name__ == "__main__":
    if len(sys.argv) != 2 and len(sys.argv) != 3:
        print("\033[1;31mUsage is %s experience\033[0m" % sys.argv[0])
        sys.exit()

    experience = sys.argv[1]

    SimDroplet = LiqSimDroplet(experience)

    SimDroplet.Compute()
        
    print("\n")
    print("LiqSimDroplet : Done\n\n")