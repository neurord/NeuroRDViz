import numpy as np
import h5py as h5
import time
import math
from Tkinter import *
from tkFileDialog   import askopenfilename
from tkMessageBox import *
import vtk


reader = vtk.vtkXMLImageDataReader()
actor = vtk.vtkActor()
renderer = vtk.vtkRenderer()
renderWindow = vtk.vtkRenderWindow()
imageDataGeometryFilter = vtk.vtkImageDataGeometryFilter()
renderer = vtk.vtkRenderer()
renderWindowInteractor = vtk.vtkRenderWindowInteractor()

simData = accessFile()

#might need to be replaced if we don't have regular grid
def Movie(ConcData,molNum, (rows, cols)):
    
    imageData = vtk.vtkImageData()
    imageData.SetDimensions(rows, cols, 1) #This is the x,y,z dimensions of image
    if vtk.VTK_MAJOR_VERSION <= 5:
        imageData.SetNumberOfScalarComponents(1)
        imageData.SetScalarTypeToDouble()
    else:
        imageData.AllocateScalars(vtk.VTK_DOUBLE, 1)
    
    dims = imageData.GetDimensions()

    # Fill every entry of the image data with "2.0"
        #^(len(simData['trial0']['simulation']['dend']['concentrations']))
    for z in range(dims[2]):
       for y in range(dims[1]):
          for x in range(dims[0]):
	     print x+(y*dims[0]), " simData:", ConcData[0,(x+(y*dims[0])),molNum]
             imageData.SetScalarComponentFromDouble(x, y, z, 0, ConcData[0,(x+(y*dims[0])),molNum])
	mapper.SetScalarRange(np.min(ConcData[:,(x+(y*dims[0])),molNum]),np.max(ConcData[:,(x+(y*dims[0])),molNum]))
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(imageDataGeometryFilter.GetOutputPort())
     
    global actor
    actor.SetMapper(mapper)
    actor.GetProperty().SetPointSize(3)
    render()     
    renderOfWindow()

# Setup rendering
def render():
    global renderer
    global actor
    renderer.AddActor(actor)
    renderer.SetBackground(1,1,1)
    renderer.ResetCamera()

def renderOfWindow(): 
    global renderWindow
    renderWindow.AddRenderer(renderer)
 
    global renderWindowInteractor
    renderWindowInteractor.SetRenderWindow(renderWindow)
    renderWindowInteractor.Initialize()
    renderWindowInteractor.Start()

#Replace the temporary getmolfunction with menu/gui (in separate file)
def TempGetMolfunction(Species)
	molChoice=rawinput('Enter Molecule Name')
	molNum = list(Species).index(molChoice)
	return

molNum=TempGetMolfunction(simData['trial0']['output']['dend']['output_species'])

#Replace arraysize spec with a function that deals with the morphology
arraysize=(3,4)

Movie(simData['trial0']['simulation']['dend']['concentrations'],molNum,arraysize)
