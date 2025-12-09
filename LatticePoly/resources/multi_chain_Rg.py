import os
import sys
import re
import numpy as np

from vtkReader import vtkReader


def detect_polymer_series(outputDir):
    """
    Detect all unique polymer prefixes from files like:
    Apoly00000.vtk, Bpoly00100.vtk, Something_poly00300.vtk, etc.
    We split on the substring 'poly' and take everything before it.
    """
    polymers = set()

    for f in os.listdir(outputDir):
        if "poly" in f and f.endswith(".vtp"):
            prefix = f.split("poly")[0]
            polymers.add(prefix)

    polymers = sorted(list(polymers))
    print("Found polymers:", polymers)
    return polymers


class RadiusGyration():

    def __init__(self, outputDir):
        self.outputDir = outputDir
        self.polymers = detect_polymer_series(outputDir)
        print(f"Processing polymers: {self.polymers}")

        # Create one reader per polymer
        self.readers = {
            P: vtkReader(outputDir, P, initFrame=initFrame, readLiq=False, readPoly=True, )
            for P in self.polymers
        }

        # Determine number of frames by inspecting one reader
        first_reader = next(iter(self.readers.values()))
        self.Nframes = first_reader.N

        # Allocate Rg table: Npolymers × Nframes
        self.Rg = np.zeros((len(self.polymers), self.Nframes), dtype=np.float32)

        # Output file
        self.outputFile = os.path.join(self.outputDir, "RgAllPolymers.res")
        if os.path.exists(self.outputFile):
            print(f"File {self.outputFile} already exists – aborting")
            sys.exit()

    def compute(self):
        for p_index, (P, reader) in enumerate(self.readers.items()):
            print(f"\nProcessing polymer {P}")
            for i in range(self.Nframes):
                data = next(reader)
                pos = data.polyPos[:]               # shape = n_monomers × 3
                pos_cm = np.mean(pos, axis=0)
                dr = pos - pos_cm
                Rg2 = np.mean(np.sum(dr * dr, axis=1))
                self.Rg[p_index, i] = Rg2

                if (i+1) % 1000 == 0:
                    print(f" Polymer {P}: processed {i+1}/{self.Nframes}")

    def save(self):
        np.savetxt(self.outputFile, self.Rg)
        print(f"\033[1;32mSaved Rg(t) for all polymers to '{self.outputFile}'\033[0m")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} outputDir")
        sys.exit()

    outputDir = sys.argv[1]
    initFrame = int(sys.argv[2])

    rg = RadiusGyration(outputDir)
    rg.compute()
    rg.save()