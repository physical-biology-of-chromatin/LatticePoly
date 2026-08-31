import time
import subprocess
import numpy as np
import datetime
import matplotlib.pyplot as plt

subprocess.run("clear", shell=True, executable="/bin/bash")

print("\n\n")

while 1:

        squeue = subprocess.run("squeue --me", shell=True, executable="/bin/bash", capture_output = True).stdout.decode('utf-8').split("\n")
        queue = sum(map(lambda x: int("-48]" in x), squeue))
        executing = len(squeue)-2
        if queue != 0:
                print(f"q {queue:^10} e {executing:^10}", end = '\r')
        else :
                if executing != 0:
                        print(f"e {executing:^10}", end = '\r')
                else :
                        print("Finito!")
        

        time.sleep(10)