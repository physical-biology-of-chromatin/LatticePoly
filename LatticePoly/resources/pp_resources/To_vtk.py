import math as m
import os

from typing import cast
import pickle
import sys

import numpy as np

import h5py
from hdf5Reader import hdf5Reader
from Liq_Droplet import Droplet, Event
from vtk import (
    vtkCellArray,
    vtkFloatArray,
    vtkIntArray,
    vtkLine,
    vtkPoints,
    vtkPolyData,
    vtkXMLPolyDataWriter,
    vtkCubeSource
)


class Hdf5_to_vtk_with_droplet_tracking:
    def __init__(self, inputDir: str, fileName: str, initFrame: int = -1):

        self.reader = hdf5Reader(inputDir, 
                     fileName, 
                     initFrame, 
                     readLiq=True, 
                     readPoly=False)

        self.boxDims = np.array(self.reader.boxDim, dtype=np.float32)
        
        self.TimeId = 0
        
        with h5py.File(os.path.join(inputDir, "process.h5"), 'r') as processFile:
            print(list(processFile.keys()))
            self.liqInfo = np.array(processFile["liq_info"], dtype = np.int32)
        
        
        
        with open(os.path.join(inputDir, "liq_droplets.pickle"), "rb") as dropletFile:
            self.droplet_dict = pickle.load(dropletFile)
    
        self.Ndroplet = len(self.droplet_dict)
        self.NstrFormat = len(str(self.Ndroplet)) + 1
        
        # self.dropletNum, self.liqToDroplet, self.dropletTraj = self.InitDroplet() 
        self.liqToDroplet, self.liqToDroplet_tau = self.InitDroplet() 

        os.chdir(inputDir)

        self.outputPath = fileName.split(".")[0] + "_vtk"
        os.makedirs(self.outputPath, exist_ok=True)
        
    def InitDroplet(self):
        
        # dropletNum = [[] for _ in range(self.reader.N)]
        liqToDroplet = np.zeros((self.reader.N, self.reader.nLiq), dtype = np.int32) - 1
        liqToDroplet_tau = np.zeros((self.reader.N, self.reader.nLiq), dtype = np.int32) - 1
        
        # dropletTraj = {} 
        compteur = 0
        for dropletId, droplet in self.droplet_dict.items():
            for node in droplet.nodes:
                # dropletNum[node[0]].append(dropletId)
                if droplet.tau > 10:
                    liqToDroplet[node[0]][self.liqInfo[node[0]] == node[1]] = dropletId
                    liqToDroplet_tau[node[0]][self.liqInfo[node[0]] == node[1]] = droplet.tau
                        
                else:
                    liqToDroplet[node[0]][self.liqInfo[node[0]] == node[1]] = -1
            # if droplet.tau > 10:
            #     print(dropletId, droplet.droplet_id, droplet.tau, droplet.sizes, droplet.nodes
            #     )
            #     sys.exit()
            #     compteur += 1
            
            # dropletTraj[dropletId] = self.ComputeDropletTraj(droplet)
            
            
        return(liqToDroplet, liqToDroplet_tau) #dropletNum, liqToDroplet, dropletTraj)
            
        
    def Print(self):
        self.PrintBox()

        for i in range(self.reader.N):
            
            data = next(self.reader)
            self.PrintLiqFrame(i, data)
    
            if (i + 1) % 10 == 0:
                print("Printed %d out of %d configurations" % (i + 1, self.reader.N))
                
    # pour frame i -> k VTK with all points as particle and color as droplet k in / out. 
    # two file per time, one with particle inside droplet k, one with the other              
    
    
        # fileDropName = "Drop{}.vtp".format(str(i + 1).zfill(self.NstrFormat))
        # fileDropPath = os.path.join(self.outputPath, fileDropName)
    
    
    def PrintLiqFrame(self, i : int, data : hdf5Reader):

        fileLiqName = 'liq{:05d}.vtp'.format(i+self.reader.initFrame)
        fileLiqPath = os.path.join(self.outputPath,fileLiqName)
        
        points = vtkPoints()
        liqDensity = vtkFloatArray()
        liqDisplacement = vtkFloatArray()
        liqDropletId = vtkFloatArray()
        liqDropletTau = vtkFloatArray()
        
        liqDensity.SetName("Density")
        liqDensity.SetNumberOfComponents(1)
        
        liqDisplacement.SetName("Displacement")
        liqDisplacement.SetNumberOfComponents(3)
        
        liqDropletId.SetName("DropletId")
        liqDropletId.SetNumberOfComponents(1)
        
        liqDropletTau.SetName("DropletTau")
        liqDropletTau.SetNumberOfComponents(1)
            
        for j in range(self.reader.nLiq):
        
            aveDensity = data.liqDens[j]

            # if self.shift != [0, 0, 0]:
            #     x = data.liqPos[j][0]+self.shift[0]
            #     if x < 0:
            #         x += self.reader.boxDim[0]
            #     elif x > self.reader.boxDim[0]:
            #         x -= self.reader.boxDim[0]
            #     y = data.liqPos[j][1]+self.shift[1]
            #     if y < 0:
            #         y += self.reader.boxDim[1]
            #     elif y > self.reader.boxDim[1]:
            #         y -= self.reader.boxDim[1]
            #     z = data.liqPos[j][2]+self.shift[2]
            #     if z < 0:
            #         z += self.reader.boxDim[2]
            #     elif z > self.reader.boxDim[2]:
            #         z -= self.reader.boxDim[2]
                

            # else:
            x = data.liqPos[j][0]
            y = data.liqPos[j][1]
            z = data.liqPos[j][2]

            
            dx = data.liqDisp[j][0]
            dy = data.liqDisp[j][1]
            dz = data.liqDisp[j][2]
                    
            points.InsertNextPoint(x, y, z)
        
            liqDensity.InsertNextValue(aveDensity)
            liqDisplacement.InsertNextTuple3(dx, dy, dz)
            liqDropletId.InsertNextValue(self.liqToDroplet[i,j])
            
            liqDropletTau.InsertNextValue(self.liqToDroplet_tau[i,j])
        
        
        polyData = vtkPolyData()
        writer = vtkXMLPolyDataWriter()

        polyData.SetPoints(points)
        
        polyData.GetPointData().AddArray(liqDensity)
        polyData.GetPointData().AddArray(liqDisplacement)
        polyData.GetPointData().AddArray(liqDropletId)
        polyData.GetPointData().AddArray(liqDropletTau)

        writer.SetFileName(fileLiqPath)
        writer.SetInputData(polyData)
        
        writer.Write()


        
        # for k in self.dropletNum[i]:
            
        #     droplet = self.droplet_dict[k]
    
        #     fileDropName = 'drop{:05d}.vtp'.format(self.TimeId)
        #     fileDropPath = os.path.join(self.outputPath,fileDropName)
        
        #     dropPoints = vtkPoints()
            
        #     dropIO = vtkIntArray()
            
        #     dropIO.SetName("IO")
        #     dropIO.SetNumberOfComponents(1)
            
        
        #     for j in range(self.reader.nLiq):
            
        #         x = data.liqPos[j][0]
        #         y = data.liqPos[j][1]
        #         z = data.liqPos[j][2]
                
        #         dropPoints.InsertNextPoint(x, y, z)
                
        #         sep = 2
                
        #         # if droplet.tau + droplet.frame_start > i:
        #         #     if self.distance(droplet.center_of_mass[i - droplet.frame_start - 1], droplet.center_of_mass[i - droplet.frame_start], self.boxDims) > 4:
        #         #         sep = 4
                
        #         if self.liqToDroplet[i][j] == k:    
        #             dropIO.InsertNextValue(int(i<self.reader.N-1 and self.liqToDroplet[i+1][j] == k)+sep)
        #         else:
        #             dropIO.InsertNextValue(int(i<self.reader.N-1 and self.liqToDroplet[i+1][j] == k))
                
            
        #     polyData = vtkPolyData()
        #     writer = vtkXMLPolyDataWriter()
        
        #     polyData.SetPoints(dropPoints)
        #     polyData.GetPointData().AddArray(dropIO)
            
        #     writer.SetFileName(fileDropPath)
        #     writer.SetInputData(polyData)
            
        #     writer.Write()
            
            
        #     fileTrajName = 'CoM{:05d}.vtp'.format(self.TimeId)
        #     fileTrajPath = os.path.join(self.outputPath,fileTrajName)
            
        #     polyDataDropletTraj = self.dropletTraj[k]
        #     writerDropletTraj = vtkXMLPolyDataWriter()
            
        #     writerDropletTraj.SetFileName(fileTrajPath)
        #     writerDropletTraj.SetInputData(polyDataDropletTraj)
        
        #     writerDropletTraj.Write()
            
            
            
        #     self.TimeId += 1
    

    
    
    def ComputeDropletTraj(self, droplet: Droplet) -> vtkPolyData:

        points = vtkPoints()
        lines = vtkCellArray()
    
        displacements = vtkFloatArray()
    
        displacements.SetName("Distance")
        displacements.SetNumberOfComponents(1)
    
        for t in range(droplet.tau):
            
            x = droplet.center_of_mass[t][0]
            y = droplet.center_of_mass[t][1]
            z = droplet.center_of_mass[t][2]
        
            points.InsertNextPoint(x, y, z)
        
        for t in range(droplet.tau - 1):
            displacement, _pbc = self.distance(
                droplet.center_of_mass[t], droplet.center_of_mass[t + 1], self.boxDims
            )
            
            if _pbc:
                line = vtkLine()
            
                line.GetPointIds().SetId(0, t)
                line.GetPointIds().SetId(1, t + 1)
                
                lines.InsertNextCell(line)
            
        
            displacements.InsertNextValue(displacement)
        
            
        # displacements.InsertNextValue(0)

        polyData = vtkPolyData()
    
        polyData.SetPoints(points)
        polyData.SetLines(lines)
    
        polyData.GetCellData().AddArray(displacements)
    
        return(polyData)

    @staticmethod
    def distance(
        CoM_1: tuple[float, float, float],
        CoM_2: tuple[float, float, float],
        BoxDims: np.ndarray[tuple[int, ...], np.dtype[np.float32]],
    ) -> tuple[float, bool]:
        
        _pbc : bool = True
        
        pDist: float = 0.0
    
        for j in range(3):
            delta = CoM_1[j] - CoM_2[j]
        
            while abs(delta) > BoxDims[j] / 2.0:
                shift = m.copysign(BoxDims[j], delta)
        
                delta -= shift
                
                _pbc = False
    
            pDist += delta**2
    
        return pDist ** (1 / 2), _pbc

    def PrintBox(self):
        fileBoxName = 'box.vtp'
        fileBoxPath = os.path.join(self.outputPath,fileBoxName)
    
        cubeSource = vtkCubeSource()
        
        L = self.reader.boxDim[0]
    
        cubeSource.SetCenter((L-0.5)/2., (L-0.5)/2., (L-0.5)/2.)
        
        cubeSource.SetXLength(L+0.5)
        cubeSource.SetYLength(L+0.5)
        cubeSource.SetZLength(L+0.5)
        
        cubeSource.Update()
    
        writer = vtkXMLPolyDataWriter()
        
        writer.SetFileName(fileBoxPath)
        writer.SetInputConnection(cubeSource.GetOutputPort())
        
        writer.Write()



if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("\033[1;31mUsage is %s inputDir fileName initFrame\033[0m" % sys.argv[0])
        sys.exit()

    inputDir = sys.argv[1]
    fileName = sys.argv[2]
    initFrame = int(sys.argv[3])

    print(sys.argv, inputDir, fileName, initFrame)

    vtk_writer = Hdf5_to_vtk_with_droplet_tracking(
    inputDir, fileName, initFrame=initFrame
    )

    vtk_writer.Print()
