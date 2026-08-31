from transfer_from_psmn import Transfer
import subprocess
import sys
import numpy as np
import pickle
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib as mpl
import os
import h5py
import re


from mpl_toolkits import axisartist
from mpl_toolkits.axes_grid1 import host_subplot

import matplotlib as mpl
cmap = mpl.colormaps['viridis']

import asyncio
from desktop_notifier import DesktopNotifier, DEFAULT_SOUND


def find_exp_name(expNum):
        
        input_path = "/Xnfs/physbiochrom/ppuel/data/"
        expName = ""
            
        res = subprocess.run(f'ssh s92node01.psmn.ens-lyon.fr ls {input_path}', shell=True, executable="/bin/bash", capture_output = True)
        
        for exp in res.stdout.decode('utf-8').split('\n'):
                
                if exp != "":
                        if 'EXP' in exp:
                                if expNum == exp.split('_')[0].split('EXP')[1] and ('liqFraction' in exp.split('_')[1] or 'bithorax' in exp.split('_')[1] or 'Lucy' in exp.split('_')[1] or 'Nazli' in exp.split('_')[1] or 'Nazly' in exp.split('_')[1]):
                                        expName = exp
        if expName == "":
                print(f"Error exp {expNum} isn't here")
                sys.exit()
        
        return(expName)

if __name__ == "__main__":
        if len(sys.argv) < 2:
                print("\033[1;31mUsage is %s expNum\033[0m" % sys.argv[0])
                
        expNum = sys.argv[1]

        expName = find_exp_name(expNum=expNum)
        output_path = os.path.join("/home/paulswann/data/", expName)
        os.makedirs(output_path, exist_ok=True)
        
        input_path = os.path.join("/Xnfs/physbiochrom/ppuel/data/", expName)

        # subprocess.run(f'ssh s92node01.psmn.ens-lyon.fr /home/ppuel/Simulation/LatticePoly/LatticePoly/.venv_bis/bin/python3 /home/ppuel/Simulation/LatticePoly/LatticePoly/resources/h5py/Is_finish.py {expNum}', shell=True, executable="/bin/bash")

        res = subprocess.run(
                f'ssh s92node01.psmn.ens-lyon.fr cat {os.path.join(input_path, "process_path.txt")}',
                shell=True, 
                executable="/bin/bash", 
                capture_output = True).stdout.decode('utf-8').split('\n')
        
        
        re_path = re.compile(r"[\D_]+_(\d+\.?\d*)")
        
        list_already_transfer = [re.findall(re_path, path)[:-1] for path in os.listdir(output_path)]
        
        for params in res:
                if params != '' and params.split(",")[:-1] not in list_already_transfer:
                        print(params.split(",")[:-1])
                        move = Transfer(expNum, params.split(","),True, True)
                        move.transfer_file()



notifier = DesktopNotifier()

async def main():
    await notifier.send(title="Transfer fini", message=f"Exp {expNum}", sound=DEFAULT_SOUND)

asyncio.run(main())
