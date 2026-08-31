import numpy as np 
import os
import sys
import matplotlib.pyplot as plt
import h5py
import conversion
import scipy.optimize

class Visualisation():
        def __init__(self, dataDir, dataset, conv = False, dim = 2):
                self.dim = dim
                self.figDir = dataDir
                self.conv = conv
                self.file = h5py.File(os.path.join(dataDir,"post_process.h5"),'r')
                print(list(self.file.keys()))
                self.dataset = dataset
                self.dataset_list = []
                file = open("resources/h5py/input_post_processing.cfg")
                for line in file.readlines():
                        self.dataset_list.append(line.strip())
                file.close()
                print(self.dataset_list)

                if not dataset in self.dataset_list:
                        print(self.dataset_list)
                        sys.exit()

                self.data_array = np.array(self.file[dataset])

                for parameters, range_p in self.file[dataset].attrs.items():
                        if self.conv == True:
                                if parameters.split("_")[0] in ['JLL','JLP','JLPP','EV']:
                                        range_p = [str(np.round(conversion.conversion_kBT_to_kJ_per_mol(float(J)),1)) for J in range_p]
                                elif parameters.split("_")[0] in ['LDENS']:
                                        range_p = [str(np.round(conversion.conversion_density_to_yM(float(D)),0)) for D in range_p]

                        if parameters.split("_")[-1] == "0":
                                self.X_name = parameters.split("_")[:-1]
                                if len(self.X_name) == 1:
                                        self.X_name = self.X_name[0]
                                else:
                                        self.X_name = "_".join(self.X_name)
                                self.X_range = range_p
                        elif parameters.split("_")[-1] == "1":
                                self.Y_name = parameters.split("_")[:-1]
                                if len(self.Y_name) == 1:
                                        self.Y_name = self.Y_name[0]
                                else:
                                        self.Y_name = "_".join(self.Y_name)
                                self.Y_range = range_p
                        elif parameters.split("_")[-1] == "2":
                                self.Z_name = parameters.split("_")[:-1]
                                if len(self.Z_name) == 1:
                                        self.Z_name = self.Z_name[0]
                                else:
                                        self.Z_name = "_".join(self.Z_name)
                                self.Z_range = range_p

                        else:
                                print("error attributs")
                
                if dataset == 'liqFraction':
                        self.data_min = 0
                        self.data_max = 1
                        self.data_array = np.where(self.data_array == -1., np.nan, self.data_array)
                        
                elif 'MSD' in dataset:
                        self.data_min = np.min(self.data_array)
                        self.data_max = np.max(self.data_array)
                elif 'DropNum' in dataset:
                        # self.data_array = np.where(self.data_array > 0 , np.log10(self.data_array), self.data_array)
                        self.data_min = np.min(self.data_array)
                        self.data_max = np.max(self.data_array)
                        # self.data_array = np.where(self.data_array == -1., np.nan, self.data_array)
                        
                elif 'liqrValue' in dataset:
                        self.data_min = np.min(self.data_array)
                        self.data_max = np.max(self.data_array)
                elif 'liqDiff' in dataset:
                        if dim == 2:
                                file = open("data/Diff.txt","r")
                                D0 = float(file.readline().strip())
                                print(D0)
                                file.close()
                                ldens_array = np.zeros((12,12))
                                for i in range(12):
                                        for j in range(12):
                                                ldens_array[i,j] = float(self.file[dataset].attrs['LDENS_1'][j])
                                self.data_array = 1/self.data_array*D0/ldens_array
                                self.data_min = 0 #np.min(self.data_array)
                                self.data_max = np.max(self.data_array)
                        elif dim == 3 and self.Y_name == "JLL_VALENCY":
                                file = open("data/Diff.txt","r")
                                D0 = float(file.readline().strip())
                                file.close()
                                ldens_array = np.zeros((12,12,12))
                                for i in range(12):
                                        for j in range(12):
                                                for k in range(12):
                                                        ldens_array[i,j,k] = float(np.array(self.file[dataset].attrs['LDENS_2'])[k])
                                self.data_array = np.where(self.data_array > 0, np.log10(1/self.data_array*D0/ldens_array), 1)
                                self.data_min = np.min(self.data_array)
                                self.data_max = np.max(self.data_array)
                        else:
                                self.data_min = np.min(self.data_array)
                                self.data_max = np.max(self.data_array)


                else:
                        self.data_array = np.where(self.data_array == -1., np.nan, self.data_array)
                        self.data_min = np.min(self.data_array)
                        self.data_max = np.max(self.data_array)

                if self.dim == 2:
                        self.visualisation_2D()
                elif self.dim == 1:
                        self.visualisation_1D()
                elif self.dim == 3:
                        self.visualisation_3D()
                else:
                        raise("error of dimension")

        def visualisation_2D(self):


                print(self.data_min, self.data_max)

                fig,ax = plt.subplots(1,1)
                im = ax.imshow(np.flip(np.transpose(self.data_array),0), vmin = self.data_min, vmax = self.data_max, cmap = 'plasma' )
                ax.set_xlabel(self.X_name)
                ax.set_ylabel(self.Y_name)
                ax.set_xticks([i for i in range(len(self.X_range))], self.X_range, rotation = 90)
                ax.set_yticks([len(self.Y_range)-i-1 for i in range(len(self.Y_range))], self.Y_range)
                cbar = fig.colorbar(im, ax=ax, shrink=0.7)
                cbar.minorticks_on()
                cbar.set_ticks([np.round(self.data_min + i/4*(self.data_max - self.data_min),1) for i in range(5)],labels = ["{:0.1f}".format(self.data_min + i/4*(self.data_max - self.data_min)) for i in range(5)],fontweight='bold', fontsize=value_fontsize)
                plt.savefig(os.path.join(self.figDir,dataset+".png"))
                plt.show()

        def visualisation_1D(self):    

                def f(x,a):
                        return(-x+a)

                res, pcov = scipy.optimize.curve_fit(f,[np.log(float(i)) for i in self.X_range], np.log(self.data_array))
                print(res)
                
                file = open("data/Diff.txt","w")
                file.write(str(np.exp(res[0])))
                print(str(np.exp(res[0])))
                file.close()
                print(1 - np.sqrt(np.diag(pcov))[0])

                print(self.data_min, self.data_max)
                fig,ax = plt.subplots(1,1)
                im = ax.plot([np.log(float(i)) for i in self.X_range], np.log(self.data_array))
                ax.plot([np.log(float(i)) for i in self.X_range],[-np.log(float(i)) + res[0] for i in self.X_range])
                ax.set_xlabel(self.X_name)
                ax.set_ylabel(self.dataset)
                ax.set_xticks([np.round(np.log(float(i)),2) for i in self.X_range], [np.round(np.log(float(i)),2) for i in self.X_range], rotation = 90)
                plt.savefig(os.path.join(self.figDir,dataset+"_log.png"))
                
                fig,ax = plt.subplots(1,1)
                im = ax.plot([float(i) for i in self.X_range], self.data_array)
                ax.plot([float(i) for i in self.X_range],[np.exp(res[0])/float(i) for i in self.X_range])
                ax.set_xlabel(self.X_name)
                ax.set_ylabel(self.dataset)
                ax.set_xticks([np.round(float(i),3) for i in self.X_range], [np.round(float(i),3) for i in self.X_range], rotation = 90)
                plt.savefig(os.path.join(self.figDir,dataset+".png"))

                
                plt.show()

        def visualisation_3D(self):
                # print(self.X_name, self.Y_name, self.Z_name)
                print(self.X_range, self.Y_range, self.Z_range)
                print(self.data_min,self.data_max)
                # print(self.data_array)
                valency = [0,4,5,6,7,8,9,10,11,1,2,3]
                # print(np.shape(self.data_array))
                fig = plt.figure()
                for i in range(12):
                        ax = fig.add_subplot(3,4,i+1)
                        im = ax.imshow(np.flip(np.transpose(self.data_array[:,valency[i],:]),0), vmin = self.data_min, vmax = self.data_max, cmap = 'plasma' )
                        ax.set_xlabel(self.X_name)
                        ax.set_ylabel(self.Z_name)
                        ax.set_xticks([i for i in range(len(self.X_range))], self.X_range, rotation = 90)
                        ax.set_yticks([len(self.Z_range)-i-1 for i in range(len(self.Z_range))], self.Z_range)
                        cbar = fig.colorbar(im, ax=ax, shrink=0.7)
                        cbar.minorticks_on()
                plt.savefig(os.path.join(self.figDir,dataset+".png"))

                plt.show()



                


if __name__ == "__main__":
        if len(sys.argv) < 6 and len(sys.argv) > 2:
                dataDir = sys.argv[1]
                dataset = sys.argv[2]
                if len(sys.argv) == 4:
                        conv = sys.argv[3]
                else:
                        conv = "False"
                if len(sys.argv) == 5:
                        dim = int(sys.argv[4])
                else:
                        dim = 2
                
                visu = Visualisation(dataDir, dataset, conv, dim)
        else:
                print("\033[1;31mUsage is %s dataDir dataset conv = False dim = 2\033[0m" % sys.argv[0])        
                sys.exit()

        







# def translation(word):
#     if word == "Jlp":
#         return("Liquid-methyl")
#     if word == "Jlpp":
#         return("Liquid-PRE")
#     if word == "Jll":
#         return("Liquid-liquid")
#     if word == "EV":
#         return("Excluded volume")
#     if word in ["Ldens","ldens"]:
#         return("Density liquid")


