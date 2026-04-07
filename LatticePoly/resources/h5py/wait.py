from subprocess import run

import time

time.sleep(60*30)

run("python3 resources/submission/submit2d_slurm.py EXP26_liqFraction_pvalue_corrige", shell=True, executable="/bin/bash")