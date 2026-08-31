import os
import sys
import subprocess
from To_vtk_with_droplet_tracking import Hdf5_to_vtk_with_droplet_tracking
from LiqDropletV2 import Droplet, Event

class Transfer():
        def __init__(self, expNum : str, expData : list[str], isPoly : bool, backInBox : bool):
                print("Transfer init")

                self.expNum = expNum
                self.expData = expData
                self.expName = ""
                self.dict_param = {}
                self.isPoly = isPoly
                self.backInBox = backInBox

                self.input_path = "/Xnfs/physbiochrom/ppuel/data/"

                self.output_path = os.path.join("/home/paulswann/data/")

                self.find_exp_name()

                print(f"Transfer find_exp_name : {self.expName}")
                
                self.input_path = os.path.join(self.input_path,self.expName)
                self.output_path = os.path.join(self.output_path,self.expName)+'/'
                
                self.find_path_name()

                
        def find_exp_name(self):
                
                res = subprocess.run(f'ssh tunnel ls {self.input_path}', shell=True, executable="/bin/bash", capture_output = True)
                
                for exp in res.stdout.decode('utf-8').split('\n'):
                        if exp != "":
                                if 'EXP' in exp:
                                        if self.expNum == exp.split('_')[0].split('EXP')[1]: # and ('liqFraction' in exp.split('_')[1] or 'bithorax' in exp.split('_')[1] or 'Lucy' in exp.split('_')[1] or 'Nazli' in exp.split('_')[1] or 'Lucy' in exp.split('_')[1] or 'Nazly' in exp.split('_')[1]):
                                                self.expName = exp
                if self.expName == "":
                        print(f"Error exp {self.expNum} isn't here")
                        sys.exit()
                

        def find_path_name(self):
                print("Transfer find_path_name")
                
                
                for i in range(len(self.expData)):
                        res = subprocess.run(f'ssh tunnel ls {self.input_path}', shell=True, executable="/bin/bash", capture_output = True).stdout.decode('utf-8').split('\n')
                        
                        print(res)

                        for file in res:
                                if ".h5" not in file and file != "" and '.png' not in file and '.txt' not in file and 'tmp' not in file and 'out' not in file: 
                                        param = file
                        
                        self.input_path = os.path.join(self.input_path,os.path.join(param,self.expData[i]))
                        
                        print(f"self.input_path {self.input_path}")
                        
                        if i == 0:
                                self.output_path = os.path.join(self.output_path,"_".join([param,self.expData[i]]))
                        else:
                                self.output_path = "_".join([self.output_path,param,self.expData[i]])

                        print(f"self.output_path {self.output_path}")
                
                print(f"Transfer {self.output_path}")

        def transfer_file(self):
                os.makedirs(self.output_path,exist_ok=True)
                self.command = f'scp tunnel:{self.input_path}/* {self.output_path}/'
                subprocess.run(self.command, shell=True, executable="/bin/bash")
                Hdf5_to_vtk_with_droplet_tracking(self.output_path, 'traj.h5', initFrame=-1).Print() #, isPoly=self.isPoly, backInBox=self.backInBox, shift=["0","0","0"]).Print()




if __name__ == "__main__":
        if len(sys.argv) < 5:
                print("\033[1;31mUsage is %s expNum 0.0 0.001 2 ... isPoly backInBox\033[0m" % sys.argv[0])
                
        expNum = sys.argv[1]
        expData = sys.argv[2:-2]
        isPoly = sys.argv[-2] == "1"
        backInBox = sys.argv[-1] == "1"

        print(sys.argv)

        move = Transfer(expNum, expData, isPoly, backInBox)

        move.transfer_file()