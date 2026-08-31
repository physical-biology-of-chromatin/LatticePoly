import sys
import os
import numpy

import h5py
import subprocess


class Process():
        def __init__(self, XnfsDir, pyCode, initframe):
                self.XnfsDir = XnfsDir
                self.pyCode = pyCode
                self.initframe = initframe

                self.list_dir = []
                for folder in os.listdir(self.XnfsDir):
                        if not(".h5" in folder):
                                self.list_dir.append(os.path.join(self.XnfsDir, folder))
                print(self.list_dir)

                self.start()

        def start(self):
                i = 0
                while len(self.list_dir) != 0 and i < 1000:
                        list_folder = os.listdir(self.list_dir[0])
                        print(self.list_dir[0])

                        if "traj.h5" in list_folder:
                                command = f".venv/bin/python3 {self.pyCode} {self.list_dir[0]} traj.h5 {self.initframe}"
                                subprocess.run(command, shell=True, executable="/bin/bash") 
                        else:
                                for folder in list_folder:
                                        if not '.h5' in folder:
                                                self.list_dir.append(os.path.join(self.list_dir[0], folder))
                        self.list_dir = self.list_dir[1:]
                        i += 1



if __name__ == "__main__":
        if len(sys.argv) != 4:
                print("\033[1;31mUsage is %s XnfsDir pyCode initframe\033[0m" % sys.argv[0])
                sys.exit()

        XnfsDir = sys.argv[1]
        pyCode = sys.argv[2]
        initframe = sys.argv[3]


        process = Process(XnfsDir, pyCode, initframe)

        process.start()
