from subprocess import run
import numpy as np

mem_tmp = np.array([line.split('Cascade')[0].replace(" ", "") for line in run("squeue --me", shell=True, executable="/bin/bash", capture_output = True).stdout.decode('utf-8').split("\n")[1:-1]])
c_max = len(mem_tmp)
c = 0
for job in mem_tmp:
        if job != None and "1-48" not in job:
                print(f"{c/c_max*100:0.1f}% {job}", end='\r')
                c+=1
                run(f"scancel {job}", shell=True, executable="/bin/bash", capture_output=True).stdout.decode()

# 1296