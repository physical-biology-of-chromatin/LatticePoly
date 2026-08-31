import numpy as np
import os
import sys
import h5py
import scipy.optimize
import scipy.stats
import numba
from LiqCluster_lifeTime import Droplet, Event
import networkx as nx
import pickle
import matplotlib as mpl
import matplotlib.pyplot as plt

plasma = mpl.colormaps['plasma'].resampled(8)

marker = ['*','o', 's', 'd']

exp = "EXP32"

class Post_processing():
        def __init__(self, dataDir, outputDir, initFrame):

                self.dataDir = dataDir

                self.initFrame = int(initFrame)
                self.file = h5py.File(os.path.join(outputDir,'post_process.h5'),'w')

                self.dataset_list = []
                print(list(os.listdir()))
                config_file = open('resources/h5py/input_post_processing.cfg','r')
                for line in config_file.readlines():
                        self.dataset_list.append(line.strip())
                config_file.close()

                os.chdir(dataDir)

                self.parameters_list = []
                self.parameters_range = {}
                self.N = 0
                self.Nmeas = 0

                self.Scan()

                os.chdir(dataDir)

                self.dim = len(self.parameters_list)-1


        def Scan(self):

                i = 1
                while i:
                        list_folder = os.listdir()
                        if "N" in list_folder:
                                i = 0
                                self.parameters_list.append(list_folder[0])
                                list_folder = os.listdir("N/")
                                self.parameters_range[self.parameters_list[-1]] = list(list_folder)
                                self.N = 20 #len(list(list_folder))
                                tmpfile = h5py.File("N/0/process.h5","r")
                                self.Nmeas = len(tmpfile[list(tmpfile.keys())[0]])
                                print(self.Nmeas)
                                tmpfile.close()
                        else:
                                if len(list_folder) == 1:
                                        self.parameters_list.append(list_folder[0])
                                else:
                                        tmp = list(list_folder)
                                        tmp.sort()
                                        self.parameters_range[self.parameters_list[-1]] = tmp

                                os.chdir(list_folder[0])
                                print(list_folder[0])

        def Aggregate(self):
                # data_array = [[[None for i in range(1)] for j in range(len(self.parameters_range[1]))] for k in range(len(self.parameters_range[0]-1))]
                self.parameters_range[self.parameters_list[1]] = [str(i*2) for i in range(1,7)]
                
                for value_0 in self.parameters_range[self.parameters_list[0]]:
                        print(value_0)
                        for value_1 in self.parameters_range[self.parameters_list[1]]:
                                for value_2 in self.parameters_range[self.parameters_list[2]]:
                                        list_droplet = []
                                        fig = plt.figure(figsize=(24,16))
                                        ax = fig.add_subplot()
                        
                                        for i in self.parameters_range['N']:
                                                try : 
                                                        file = open(f"{self.dataDir}/{self.parameters_list[0]}/{value_0}/{self.parameters_list[1]}/{value_1}/{self.parameters_list[2]}/{value_2}/N/{i}/liq_droplets.pickle", "rb")
                                                        droplet_dict = pickle.load(file)
                                                        file.close()
                                                        max_tau = -1
                                                        for droplet in droplet_dict.values():
                                                                max_tau = max(max_tau, droplet.tau)
                                                        rainbow = mpl.colormaps['rainbow'].resampled(max_tau)
                        
                                                        for droplet in droplet_dict.values():
                                                                if droplet.tau > 1 and np.max(droplet.sizes) > 4:
                                                                        list_droplet.append([droplet.tau,np.max(droplet.sizes), np.std(droplet.sizes)/np.mean(droplet.sizes)])
                                                                        ax.scatter([droplet.nodes[0][0]+ i for i in range(len(droplet.sizes))], droplet.sizes, color = rainbow(droplet.tau), s = 100)
                                                                        ax.plot([droplet.nodes[0][0]+ i for i in range(len(droplet.sizes))], droplet.sizes, color = rainbow(droplet.tau), linewidth = 10)
                                                        list_droplet = np.array(list_droplet)
                                                        print(value_0, value_1, len(list_droplet))
                                                except:
                                                        print(value_0, value_1, value_2)

                                        # data_array[self.parameters_range[self.parameters_list[0]].index(value_0)-1,
                                        #                                 self.parameters_range[self.parameters_list[1]].index(value_1),
                                        #                                 self.parameters_range[self.parameters_list[2]].index(value_2)] = list_droplet
                                        ax.set_xlabel("Time")
                                        ax.set_ylabel("Size")
                                        # ax.set_xlim(0,10)
                                        # ax.set_ylim(0,10)
                                        fig.savefig(f"/home/ppuel/data/{exp}/droplet_tracking/fig_droplet_tracking_jll_{value_0}_jll_valency_{value_1}_ldens_{value_2}.png")
                                        plt.close(fig)


if __name__ == "__main__":
        if len(sys.argv) != 4:
                print("\033[1;31mUsage is %s dataDir outputDir initFrame\033[0m" % sys.argv[0])
                sys.exit()

        dataDir = sys.argv[1]
        outputDir = sys.argv[2]
        initFrame = sys.argv[3]

        process = Post_processing(dataDir, outputDir, initFrame)

        process.Aggregate()

