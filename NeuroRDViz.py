from __future__ import print_function, division
import time
import os
os.environ['ETS_TOOLKIT'] = 'qt4'
os.environ['QT_API'] = 'pyqt'
import sip
sip.setapi('Qstring',2)
from pyface.qt import QtGui, QtCore

from traits.api import HasTraits, Instance, on_trait_change
from traitsui.api import View, Item
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, \
        SceneEditor

import numpy as np
from tvtk.api import tvtk
from mayavi.scripts import mayavi2
from mayavi import mlab
from mayavi.sources.vtk_data_source import VTKDataSource

import h5py as h5
import sys

################################################################################
#The actual visualization
class Visualization(HasTraits):
    scene = Instance(MlabSceneModel, ())

    @on_trait_change('scene.activated')
    def update_plot(self):
        # This function is called when the view is opened. We don't
        # populate the scene when the view is not yet open, as some
        # VTK features require a GLContext.
        print("In update_plot")
        ug=create_morphology(simData)
        surf = mlab.pipeline.surface(ug, opacity=0.1)

        self.scene.mlab.pipeline.surface(mlab.pipeline.extract_edges(surf), color=(0, 0, 0), )
        #anim(create_morphology(), 0)

    # the layout of the dialog created
    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=250, width=300),
                resizable=True # We need this to resize with the parent widget
                )


################################################################################

def create_morphology(simData):

    grid = np.array(getMorphologyGrid()).view(np.recarray)                               #change
    xmin = np.min((grid.x0, grid.x1, grid.x2, grid.x3), axis=0)
    ymin = np.min((grid.y0, grid.y1, grid.y2, grid.y3), axis=0)
    zmin = np.min((grid.z0, grid.z1, grid.z2, grid.z3), axis=0)
    xmax = np.max((grid.x0, grid.x1, grid.x2, grid.x3), axis=0)
    ymax = np.max((grid.y0, grid.y1, grid.y2, grid.y3), axis=0)
    zmax = np.max((grid.z0, grid.z1, grid.z2, grid.z3), axis=0)

    points = np.array((
         (xmin, ymin, zmin), (xmax, ymin, zmin), (xmax, ymax, zmin), (xmin, ymax, zmin),
         (xmin, ymin, zmax), (xmax, ymin, zmax), (xmax, ymax, zmax), (xmin, ymax, zmax)))
    points = points.swapaxes(0, 2).swapaxes(1, 2)
    points = points.reshape(-1, 3)

    voxels = np.arange(points.shape[0]).reshape(-1, 8)

    voxel_type = tvtk.Hexahedron().cell_type # VTK_HEXAHEDRON == 12
    ug = tvtk.UnstructuredGrid(points=points) 
    ug.set_cells(voxel_type, voxels)

    #This sets them all to zeros 
    #numVoxels = len(grid.x0)                                                                           #can be redone in setDef
    #concentrations = np.zeros(numVoxels)   # Have this pass less through it ********************

    #ug.point_data.scalars = np.repeat(concentrations, 8)
    #ug.point_data.scalars.name = 'concentrations'

    return ug

def get_voxel_molecule_concs(ms, simData, molecule):
    sets = simData['model']['output'].keys() #just once
    gridpoints=len(getMorphologyGrid())       #justonce
    outputSetSnapshot=np.zeros(0)
    for each in sets[1:]:
        molnum = moleculeToNumber(molecule, each, simData)        #once per molecule (not once per timestep)
        if molnum > -1:
            tempSnapshot = simData['trial0']['output'][each]['population'][ms,:,molnum]
            if len(tempSnapshot) == gridpoints:
                outputSetSnapshot=tempSnapshot
                return outputSetSnapshot
            else:
                outputSetSnapshot=np.concatenate((outputSetSnapshot,tempSnapshot),axis=0)
    if len(outputSetSnapshot)==gridpoints:
        tempSnapshot=np.zeros(gridpoints-len(outputSetSnapshot))
        outputSetSnapshot=np.concatenate((outputSetSnapshot,tempSnapshot),axis=0)
    #if outputSetSnapshot == []:                                                            #change! - account for molecule non-existent in both sets
    #    return __main__set
    return outputSetSnapshot

class MayaviQWidget(QtGui.QWidget):
    animator = None
    currentFrame = 0
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QVBoxLayout(self)
        layout.setContentsMargins(0,0,0,0)
        layout.setSpacing(0)
        self.visualization = Visualization()
        print("In MayaviQWidget init")

        # If you want to debug, beware that you need to remove the Qt
        # input hook.
        #QtCore.pyqtRemoveInputHook()
        #import pdb ; pdb.set_trace()
        #QtCore.pyqtRestoreInputHook()

        # The edit_traits call will generate the widget to embed.
        self.ui = self.visualization.edit_traits(parent=self, kind='subpanel').control
        layout.addWidget(self.ui)
        self.ui.setParent(self)

    def molecule_selected(self, text):
        if text != None: #This is necessary otherwise ~anim-loop-obect will be instantiated initially with whatever conditions passed
            if self.animator != None:
                self.animator.close()
            self.animator = anim(create_morphology(simData), simData, text, self)

    def setCurrentFrame(self, frame):
        self.currentFrame = frame

    def getCurrentFrame(self):
        return self.currentFrame

@mlab.animate(delay=10) #This is a "decorator" - dynamically alters method function w/o need for subclasses
def anim(ug, simData, molecule, frameTracker):
    f = mlab.gcf()
    iterations = len(simData['trial0']['output']['all']['times'])                                             #change! all to chosen outputSets from getSomething
    currentFrame = frameTracker.getCurrentFrame()
    while currentFrame < iterations:
        concentrations = get_voxel_molecule_concs(currentFrame, simData, molecule)
        ug.point_data.scalars = np.repeat(concentrations, 8) #make 8 to static variable of "voxelpts"
        ug.point_data.scalars.name = 'concentrations'
        surf = mlab.pipeline.surface(ug, opacity=0.1)
        mlab.pipeline.surface(mlab.pipeline.extract_edges(surf), color=(0, 0, 0))

        f.scene.render()
        currentFrame += 1
        frameTracker.setCurrentFrame(currentFrame)
        yield


def getMoleculeList(simData):
    return simData['model']['species']
                                            
def get_h5simData(fileName):
    simData = h5.File(fileName,"r")
    return simData
    
def getMorphologyGrid():
    return simData['model']['grid']

#def getSets():
    #sets = 
    #return sets

def moleculeToNumber(molecule, outputSet, simData):
    indices=np.where(simData['model']['output'][outputSet]['species'][:]==molecule)[0]
    if len(indices) == 1:
        return indices[0]
    else:
        return -1              
 
if __name__ == "__main__":
    # Don't create a new QApplication, it would unhook the Events
    # set by Traits on the existing QApplication. Simply use the
    # '.instance()' method to retrieve the existing one.
    app = QtGui.QApplication.instance()
    '''GUI = Window()'''
    container = QtGui.QWidget()
    #container.setWindowTitle("NeuoRD Visualization")
    # define a "complex" layout to test the behaviour
    layout = QtGui.QGridLayout(container)
    comboBox = QtGui.QComboBox(container)
    try:
        fileName=fname
    except NameError:
        fileName = sys.argv[1]
    simData = get_h5simData(fileName)
    moleculeList = getMoleculeList(simData) #simdata.... then past just the list below
    for moleculeType in range(len(moleculeList)):
        comboBox.addItem(moleculeList[moleculeType])
    layout.addWidget(comboBox, 0, 0)  # 0,0 = top left widget location, 0,1 = one to the right of it, etc.
    mayavi_widget = MayaviQWidget(container)
    comboBox.activated[str].connect(mayavi_widget.molecule_selected)
    layout.addWidget(mayavi_widget, 1, 1) # This is visualization of morphology
    container.show()
    window = QtGui.QMainWindow()
    window.setCentralWidget(container)
    window.show()
    print("after window shown")
    app.exec_() # Start the main event loop.
