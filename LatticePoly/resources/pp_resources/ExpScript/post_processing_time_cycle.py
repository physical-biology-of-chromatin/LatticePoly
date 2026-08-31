import numpy as np
import os
import sys
import h5py
import scipy.stats


class Post_processing():
        def __init__(self, dataDir, outputDir, initFrame):

                self.initFrame = int(initFrame)
                self.file = h5py.File(os.path.join(outputDir,'time_cycle.h5'),'w')

                
                os.chdir(dataDir)

                self.parameters_list = []
                self.parameters_range = {}
                self.N = 0

                self.Scan()

                os.chdir(dataDir)

        def Scan(self):

                i = 1
                while i:
                        list_folder = os.listdir()
                        if "N" in list_folder:
                                i = 0
                                self.parameters_list.append(list_folder[0])
                                list_folder = os.listdir("N/")
                                self.parameters_range[self.parameters_list[-1]] = list(list_folder)
                                self.N = len(list(list_folder))
                        else:
                                if len(list_folder) == 1:
                                        self.parameters_list.append(list_folder[0])
                                else:
                                        tmp = list(list_folder)
                                        tmp.sort()
                                        self.parameters_range[self.parameters_list[-1]] = tmp

                                os.chdir(list_folder[0])
                                print(list_folder[0])


        def Time(self):

                data_array = np.zeros(tuple([len(self.parameters_range[parameter]) for parameter in self.parameters_list][:-1]))

                for value_0 in self.parameters_range[self.parameters_list[0]]:
                        for value_1 in self.parameters_range[self.parameters_list[1]]:
                                for i in self.parameters_range['N']:
                                        total_time = 0
                                        log_file = open(f"{self.parameters_list[0]}/{value_0}/{self.parameters_list[1]}/{value_1}/N/{i}/log.out")
                                        for line in log_file:
                                                if line.split(" ")[0] == 'Total':
                                                        total_time = float(line.split(" ")[2])
                                        log_file.close()
                                        data_array[self.parameters_range[self.parameters_list[0]].index(value_0),
                                                   self.parameters_range[self.parameters_list[1]].index(value_1)] += total_time/self.N
                                        
                
                self.file.create_dataset("time",data=data_array)
                i = 0
                for parameter in self.parameters_list[:-1]:
                        self.file["time"].attrs[parameter+f"_{i}"] = self.parameters_range[parameter]
                        i += 1
                print("time done")


if __name__ == "__main__":
        if len(sys.argv) != 4:
                print("\033[1;31mUsage is %s dataDir outputDir initFrame\033[0m" % sys.argv[0])
                sys.exit()

        dataDir = sys.argv[1]
        outputDir = sys.argv[2]
        initFrame = sys.argv[3]

        process = Post_processing(dataDir, outputDir, initFrame)

        process.Time()

