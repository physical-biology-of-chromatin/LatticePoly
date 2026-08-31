import numpy as np
import os
import sys
import h5py
import scipy.optimize
import scipy.stats
import numba


class Post_processing():
        def __init__(self, dataDir, outputDir, initFrame):

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

        def Write(self):
                def f(x,a):
                        return(x*a)
                if self.dim == 1:
                        for dataset in self.dataset_list:
                                data_array = np.zeros(tuple([len(self.parameters_range[parameter]) for parameter in self.parameters_list][:-1]))
                                for value_0 in self.parameters_range[self.parameters_list[0]]:
                                        print(value_0)
                                        if dataset == "liqDiff":
                                                liqMSD_CoM = np.zeros(self.Nmeas)
                                                error = 1
                                                for i in self.parameters_range['N']:
                                                        try:
                                                                h5file = h5py.File(f"{self.parameters_list[0]}/{value_0}/N/{i}/process.h5",'r')
                                                                liqMSD_CoM += np.array(h5file["liqMSD_CoM"])/self.N
                                                                h5file.close()
                                                        except:
                                                                data_array[self.parameters_range[self.parameters_list[0]].index(value_0)] = -1
                                                                error = 0
                                                if error:
                                                        res = scipy.optimize.curve_fit(f, xdata = [i for i in range(self.Nmeas)], ydata = liqMSD_CoM)
                                                        data_array[self.parameters_range[self.parameters_list[0]].index(value_0)] = res[0][0]
                                                        

                                        # elif dataset == "liqrValue":
                                        #         liqMSD_CoM = np.zeros(self.Nmeas)
                                        #         error = 1
                                        #         for i in self.parameters_range['N']:
                                        #                 try:
                                        #                         h5file = h5py.File(f"{self.parameters_list[0]}/{value_0}/N/{i}/process.h5",'r')
                                        #                         liqMSD_CoM += np.array(h5file["liqMSD_CoM"])/self.N
                                        #                         h5file.close()
                                        #                 except:
                                        #                         data_array[self.parameters_range[self.parameters_list[0]].index(value_0)] = -1
                                        #                         error = 0
                                        #         if error:
                                        #                 res = scipy.stats.linregress(xdata = [i for i in range(self.Nmeas)], ydata = liqMSD_CoM)
                                        #                 data_array[self.parameters_range[self.parameters_list[0]].index(value_0)] = res.rvalue                                                             

                                        
                                        elif "MSD" in dataset:
                                                for i in self.parameters_range['N']:
                                                        try:
                                                                h5file = h5py.File(f"{self.parameters_list[0]}/{value_0}/{self.parameters_list[1]}/{value_1}/N/{i}/process.h5","r")
                                                                data_array[self.parameters_range[self.parameters_list[0]].index(value_0)] += np.mean(np.array(h5file[dataset][:-self.initFrame]))/self.N
                                                                h5file.close()
                                                        except:
                                                                data_array[self.parameters_range[self.parameters_list[0]].index(value_0)] = -1
                                                                break
                                        else:
                                                for i in self.parameters_range['N']:
                                                        try:
                                                                h5file = h5py.File(f"{self.parameters_list[0]}/{value_0}/{self.parameters_list[1]}/{value_1}/N/{i}/process.h5","r")
                                                                data_array[self.parameters_range[self.parameters_list[0]].index(value_0)] += np.mean(np.array(h5file[dataset][self.initFrame:]))/self.N
                                                                h5file.close()
                                                        except:
                                                                data_array[self.parameters_range[self.parameters_list[0]].index(value_0)] = -1
                                                                break
                                                                
                                self.file.create_dataset(dataset,data=data_array)
                                i = 0
                                for parameter in self.parameters_list[:-1]:
                                        self.file[dataset].attrs[parameter+f"_{i}"] = self.parameters_range[parameter]
                                        i += 1
                                print(f"{dataset} done")

                elif self.dim == 2:

                        for dataset in self.dataset_list:

                                data_array = np.zeros(tuple([len(self.parameters_range[parameter]) for parameter in self.parameters_list][:-1]))

                                for value_0 in self.parameters_range[self.parameters_list[0]]:
                                        print(value_0)
                                        for value_1 in self.parameters_range[self.parameters_list[1]]:
                                                if dataset == "liqDiff":
                                                        liqMSD_CoM = np.zeros(self.Nmeas)
                                                        error = 1
                                                        for i in self.parameters_range['N']:
                                                                try:
                                                                        h5file = h5py.File(f"{self.parameters_list[0]}/{value_0}/{self.parameters_list[1]}/{value_1}/N/{i}/process.h5",'r')
                                                                        liqMSD_CoM += np.array(h5file["liqMSD_CoM"])/self.N
                                                                        h5file.close()
                                                                except:
                                                                       data_array[self.parameters_range[self.parameters_list[0]].index(value_0),
                                                                                  self.parameters_range[self.parameters_list[1]].index(value_1)] = -1
                                                                       error = 0
                                                        if error:
                                                                res = scipy.optimize.curve_fit(xdata = [i for i in range(self.Nmeas)], ydata = liqMSD_CoM)
                                                                data_array[self.parameters_range[self.parameters_list[0]].index(value_0),
                                                                           self.parameters_range[self.parameters_list[1]].index(value_1)] = res[0][0]
                                                                

                                                # elif dataset == "liqrValue":
                                                #         liqMSD_CoM = np.zeros(self.Nmeas)
                                                #         error = 1
                                                #         for i in self.parameters_range['N']:
                                                #                 try:
                                                #                         h5file = h5py.File(f"{self.parameters_list[0]}/{value_0}/{self.parameters_list[1]}/{value_1}/N/{i}/process.h5",'r')
                                                #                         liqMSD_CoM += np.array(h5file["liqMSD_CoM"])/self.N
                                                #                         h5file.close()
                                                #                 except:
                                                #                         data_array[self.parameters_range[self.parameters_list[0]].index(value_0),
                                                #                                    self.parameters_range[self.parameters_list[1]].index(value_1)] = -1
                                                #                         error = 0
                                                #         if error:
                                                #                 res = scipy.stats.linregress(xdata = [i for i in range(self.Nmeas)], ydata = liqMSD_CoM)
                                                #                 data_array[self.parameters_range[self.parameters_list[0]].index(value_0),
                                                #                            self.parameters_range[self.parameters_list[1]].index(value_1)] = res.rvalue                                                             

                                               
                                                elif "MSD" in dataset:
                                                        for i in self.parameters_range['N']:
                                                                try:
                                                                        h5file = h5py.File(f"{self.parameters_list[0]}/{value_0}/{self.parameters_list[1]}/{value_1}/N/{i}/process.h5","r")
                                                                        data_array[self.parameters_range[self.parameters_list[0]].index(value_0),
                                                                                   self.parameters_range[self.parameters_list[1]].index(value_1)] += np.mean(np.array(h5file[dataset][:-self.initFrame]))/self.N
                                                                        h5file.close()
                                                                except:
                                                                        data_array[self.parameters_range[self.parameters_list[0]].index(value_0),
                                                                        self.parameters_range[self.parameters_list[1]].index(value_1)] = -1
                                                                        break
                                                else:
                                                        for i in self.parameters_range['N']:
                                                                try:
                                                                        h5file = h5py.File(f"{self.parameters_list[0]}/{value_0}/{self.parameters_list[1]}/{value_1}/N/{i}/process.h5","r")
                                                                        data_array[self.parameters_range[self.parameters_list[0]].index(value_0),
                                                                                   self.parameters_range[self.parameters_list[1]].index(value_1)] += np.mean(np.array(h5file[dataset][self.initFrame:]))/self.N
                                                                        h5file.close()
                                                                except:
                                                                        data_array[self.parameters_range[self.parameters_list[0]].index(value_0),
                                                                                   self.parameters_range[self.parameters_list[1]].index(value_1)] = -1
                                                                        break
                                                                
                                self.file.create_dataset(dataset,data=data_array)
                                i = 0
                                for parameter in self.parameters_list[:-1]:
                                        self.file[dataset].attrs[parameter+f"_{i}"] = self.parameters_range[parameter]
                                        i += 1
                                print(f"{dataset} done")

                elif self.dim == 3: 

                        for dataset in self.dataset_list:

                                data_array = np.zeros(tuple([len(self.parameters_range[parameter]) for parameter in self.parameters_list][:-1]))

                                for value_0 in self.parameters_range[self.parameters_list[0]]:
                                        print(value_0)
                                        for value_1 in self.parameters_range[self.parameters_list[1]]:
                                                for value_2 in self.parameters_range[self.parameters_list[2]]:
                                                        if dataset == "liqDiff":
                                                                liqMSD_CoM = np.zeros((self.Nmeas))
                                                                error = 1
                                                                for i in self.parameters_range['N']:
                                                                        try:
                                                                                h5file = h5py.File(f"{self.parameters_list[0]}/{value_0}/{self.parameters_list[1]}/{value_1}/{self.parameters_list[2]}/{value_2}/N/{i}/process.h5",'r')
                                                                                liqMSD_CoM += np.array(h5file["liqMSD_CoM"])/self.N
                                                                                h5file.close()
                                                                                
                                                                        except:
                                                                                data_array[self.parameters_range[self.parameters_list[0]].index(value_0),
                                                                                           self.parameters_range[self.parameters_list[1]].index(value_1),
                                                                                           self.parameters_range[self.parameters_list[2]].index(value_2)] = -1
                                                                                error = 0
                                                                if error:
                                                                        res = scipy.optimize.curve_fit(f, xdata = [i for i in range(self.Nmeas)], ydata = liqMSD_CoM)
                                                                        
                                                                        data_array[self.parameters_range[self.parameters_list[0]].index(value_0),
                                                                                   self.parameters_range[self.parameters_list[1]].index(value_1),
                                                                                   self.parameters_range[self.parameters_list[2]].index(value_2)] = res[0][0]
                                                                        
                                 
                                                                                                                                
                                                        # elif dataset == "liqrValue":
                                                        #         liqMSD_CoM = np.zeros((self.Nmeas))
                                                        #         error = 1
                                                        #         for i in self.parameters_range['N']:
                                                        #                 try:
                                                        #                         h5file = h5py.File(f"{self.parameters_list[0]}/{value_0}/{self.parameters_list[1]}/{value_1}/{self.parameters_list[2]}/{value_2}/N/{i}/process.h5",'r')
                                                        #                         liqMSD_CoM += np.array(h5file["liqMSD_CoM"])/self.N
                                                        #                         h5file.close()
                                                                                
                                                        #                 except:
                                                        #                         data_array[self.parameters_range[self.parameters_list[0]].index(value_0),
                                                        #                                    self.parameters_range[self.parameters_list[1]].index(value_1),
                                                        #                                    self.parameters_range[self.parameters_list[2]].index(value_2)] = -1
                                                        #                         error = 0
                                                        #         if error:
                                                        #                 res = scipy.stats.linregress(xdata = [i for i in range(self.Nmeas)], ydata = liqMSD_CoM)
                                                        #                 data_array[self.parameters_range[self.parameters_list[0]].index(value_0),
                                                        #                            self.parameters_range[self.parameters_list[1]].index(value_1),
                                                        #                            self.parameters_range[self.parameters_list[2]].index(value_2)] = res.rvalue                                                             


                                                        elif 'MSD' in dataset:
                                                                for i in self.parameters_range['N']:
                                                                        try:
                                                                                h5file = h5py.File(f"{self.parameters_list[0]}/{value_0}/{self.parameters_list[1]}/{value_1}/{self.parameters_list[2]}/{value_2}/N/{i}/process.h5","r")
                                                                                data_array[self.parameters_range[self.parameters_list[0]].index(value_0),
                                                                                           self.parameters_range[self.parameters_list[1]].index(value_1),
                                                                                           self.parameters_range[self.parameters_list[2]].index(value_2)] += np.mean(np.array(h5file[dataset][:-self.initFrame]))/self.N
                                                                                h5file.close()
                                                                        except:
                                                                                data_array[self.parameters_range[self.parameters_list[0]].index(value_0),
                                                                                           self.parameters_range[self.parameters_list[1]].index(value_1),
                                                                                           self.parameters_range[self.parameters_list[2]].index(value_2)] = -1
                                                                                break
                                                        
                                                        else:
                                                                for i in self.parameters_range['N']:
                                                                        try:
                                                                                h5file = h5py.File(f"{self.parameters_list[0]}/{value_0}/{self.parameters_list[1]}/{value_1}/{self.parameters_list[2]}/{value_2}/N/{i}/process.h5","r")
                                                                                data_array[self.parameters_range[self.parameters_list[0]].index(value_0),
                                                                                           self.parameters_range[self.parameters_list[1]].index(value_1),
                                                                                           self.parameters_range[self.parameters_list[2]].index(value_2)] += np.mean(np.array(h5file[dataset][self.initFrame:]))/self.N
                                                                                h5file.close()
                                                                        except:
                                                                                data_array[self.parameters_range[self.parameters_list[0]].index(value_0),
                                                                                           self.parameters_range[self.parameters_list[1]].index(value_1),
                                                                                           self.parameters_range[self.parameters_list[2]].index(value_2)] = -1
                                                                                break
                                                                
                                self.file.create_dataset(dataset,data=data_array)
                                

                                i = 0
                                for parameter in self.parameters_list[:-1]:
                                        self.file[dataset].attrs[parameter+f"_{i}"] = self.parameters_range[parameter]
                                        i += 1
                                print(f"{dataset} done")



                else:
                        raise("dimension != 2,3")



if __name__ == "__main__":
        if len(sys.argv) != 4:
                print("\033[1;31mUsage is %s dataDir outputDir initFrame\033[0m" % sys.argv[0])
                sys.exit()

        dataDir = sys.argv[1]
        outputDir = sys.argv[2]
        initFrame = sys.argv[3]

        process = Post_processing(dataDir, outputDir, initFrame)

        process.Write()

