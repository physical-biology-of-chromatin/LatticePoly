import os
import sys

from mpi4py import MPI
from utils import getInputParam


comm = MPI.COMM_WORLD
mpiRank = comm.Get_rank()

if len(sys.argv) != 3:
	if mpiRank == 0:
		print("\033[1;31mUsage is %s executable inputFile\033[0m" % sys.argv[0])
		
	sys.exit()

execFile = sys.argv[1]
inputFile = sys.argv[2]

logFile = "log.out"
errFile = "log.err"

outputDir = getInputParam("outputDir", inputFile)

outputDir = os.path.join(outputDir, str(mpiRank))
outputFile = os.path.join(outputDir, os.path.basename(inputFile))

logFile = os.path.join(outputDir, logFile)
errFile = os.path.join(outputDir, errFile)

if not os.path.isdir(outputDir):
	os.makedirs(outputDir)
	
subst = r"s|\(outputDir[[:space:]]*=[[:space:]]*\)\(.*;\)|\1" + outputDir + " ;|;"

os.system('sed -e "' + subst + '" < ' + inputFile + ' > ' + outputFile)
os.system(execFile + ' ' + outputFile + ' > ' + logFile + ' 2> ' + errFile)
