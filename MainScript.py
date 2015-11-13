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

def Movie(molNum, (rows, cols)):
    
    simData = accessFile()

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
	     print x+(y*dims[0]), " simData:", simData['trial0']['simulation']['dend']['concentrations'][0,(x+(y*dims[0])),molNum]
             imageData.SetScalarComponentFromDouble(x, y, z, 0, simData['trial0']['simulation']['dend']['concentrations'][0,(x+(y*dims[0])),molNum])
	
    mapper.SetScalarRange(np.min(simData['trial0']['simulation']['dend']
    ['concentrations'][:,(x+(y*dims[0])),molNum]),np.max(simData['trial0']['simulation']['dend']['concentrations'][:,(x+(y*dims[0])),molNum]))
	#time.sleep(0.0001)
	#^Create "function('concarray') 3-concarray = simData[......]

    writeFile(filename, imageData) 
    readFile(filename)
    convertImageToPoly()
    render()     
    renderOfWindow()

#Interprets morphology from HDF5 format to one vtk can use
#def readMorphology(simDataGrid):
    

# Write the file
def writeFile(filename, imageData):
    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(filename)
    if vtk.VTK_MAJOR_VERSION <= 5:
        writer.SetInputConnection(imageData.GetProducerPort())
    else:
        writer.SetInputData(imageData)
    writer.Write()


# Read the file (to test that it was written correctly)
def readFile(filename):
    global reader 
    reader.SetFileName(filename)
    reader.Update()


# Convert the image to a polydata
def convertImageToPoly():
    global imageDataGeometryFilter
    imageDataGeometryFilter.SetInputConnection(reader.GetOutputPort())
    imageDataGeometryFilter.Update()
     
    mapper = vtk.vtkPolyDataMapper()
    #mapper.SetScalarRange(2311,2421) #IMPORTANT!!! Set this (dynamically?) to min/max of simulation. 
    mapper.SetInputConnection(imageDataGeometryFilter.GetOutputPort())
     
    global actor
    actor.SetMapper(mapper)
    actor.GetProperty().SetPointSize(3)


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


def select_molecule(simData):
   master = Tk()

   var = StringVar(master)
   var.set("Select Molecule...") # default value
   ###Selection between dend/soma/"output_species???" needs to be clarified/fixed
   molOptions = list(simData['trial0']['output']['dend']['output_species'])

   w = apply(OptionMenu, (master, var) + tuple(molOptions))
   w.pack()

   def Run():
    #might need preceding 'print' on line below
    locateMolIndex(var.get())
    master.quit()

   button = Button(master, text="Display", command=Run)
   button.pack()


def search_molecule():
   
   #Needs to be local 
   def Run():
      locateMolIndex(e.get())
      master.quit()

   groot = Tk()
   groot.title('Enter Molecule')


   e = Entry(groot)
   e.pack()
   e.focus_set()

   b = Button(groot,text='Display',command=Run)
   b.pack(side='top')

def accessFile():
    simData = h5.File("PurkdifmodelCopy.h5","r")
    return simData

##Takes user input and returns respective index of concentration data
def locateMolIndex(molChoice):

    simData = accessFile()
    molNum = list(simData['trial0']['output']['dend']['output_species']).index(molChoice)
    print "Displaying data for " + simData['trial0']['output']['output_species'][molNum]
    Movie(molNum, (3,4))
    #return molNum
      
    
root = Tk()
menu = Menu(root)
root.config(menu=menu)
moleculeMenu = Menu(menu)
menu.add_cascade(label="Molecule", menu=moleculeMenu)
moleculeMenu.add_command(label="Select", command=select_molecule)
moleculeMenu.add_command(label="Search...", command=search_molecule)
moleculeMenu.add_separator()
#moleculeMenu.add_command(label="Exit", command=root.quit)

'''helpmenu = Menu(menu)
menu.add_cascade(label="Help", menu=helpmenu)'''
