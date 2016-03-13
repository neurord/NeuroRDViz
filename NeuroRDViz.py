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
        self.ug=create_morphology(simData)  
        surf = mlab.pipeline.surface(self.ug, opacity=0.1)

        self.scene.mlab.pipeline.surface(mlab.pipeline.extract_edges(surf), color=(0, 0, 0), )
        #anim(create_morphology(), 0)
        
    # the layout of the dialog created
    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=250, width=300),
                resizable=True # We need this to resize with the parent widget
                )


################################################################################

def create_morphology(simData):
 
    grid = np.array(simData['trial0']['model']['grid']).view(np.recarray)
    xmin = np.min((grid.x0, grid.x1, grid.x2, grid.x3), axis=0)
    ymin = np.min((grid.y0, grid.y1, grid.y2, grid.y3), axis=0)
    zmin = np.min((grid.z0, grid.z1, grid.z2, grid.z3), axis=0)
    xmax = np.max((grid.x0, grid.x1, grid.x2, grid.x3), axis=0)
    ymax = np.max((grid.y0, grid.y1, grid.y2, grid.y3), axis=0)
    zmax = np.max((grid.z0, grid.z1, grid.z2, grid.z3), axis=0)
    zmax[:]=0.5
    zmin[:]=0.0
    points = np.array((
         (xmin, ymin, zmin), (xmax, ymin, zmin), (xmax, ymax, zmin), (xmin, ymax, zmin),
         (xmin, ymin, zmax), (xmax, ymin, zmax), (xmax, ymax, zmax), (xmin, ymax, zmax)))
    points = points.swapaxes(0, 2).swapaxes(1, 2)
    points = points.reshape(-1, 3)

    voxels = np.arange(points.shape[0]).reshape(-1, 8)

    voxel_type = tvtk.Hexahedron().cell_type # VTK_HEXAHEDRON == 12
    ug = tvtk.UnstructuredGrid(points=points)
    ug.set_cells(voxel_type, voxels)
    
    concentrations = get_voxel_molecule_concs(0, simData, 0)   # Have this pass less through it ********************
    #temperature = np.repeat(colors, 8)

    ug.point_data.scalars = np.repeat(concentrations, 8)
    ug.point_data.scalars.name = 'concentrations'
    # Some vectors.
    #ug1.point_data.vectors = velocity
    #ug1.point_data.vectors.name = 'velocity'
    
    return ug

def get_voxel_molecule_concs(ms, simData, molnum):
    '''Must make the following generic for instances with A. Varying sizes. B.No Soma  etc..''' #################################<=======
    dendSnapshot = simData['trial0']['simulation']['dend']['concentrations'][ms,:,molnum] #takes [milisecond, all voxel's of it's data, molecule of interest)
    somaSnapshot = simData['trial0']['simulation']['soma']['concentrations'][ms,:,molnum] 
    wholeCellSnapshot = np.concatenate((dendSnapshot,somaSnapshot),axis=0)
    return wholeCellSnapshot

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
        #simData = get_simData(get_fileName())       
        molnum = molecule_to_number(text, simData)
        print("in molecule_selected", text, molnum)
        if molnum != None: #This is necessary otherwise ~anim-loop-obect will be instantiated initially with whatever conditions passed
            if self.animator != None:
                self.animator.close()
            self.animator = anim(create_morphology(simData), simData, molnum, self)
            
    def setCurrentFrame(self, frame):
        self.currentFrame = frame

    def getCurrentFrame(self):
        return self.currentFrame
        
@mlab.animate(delay=10) #This is a "decorator" - dynamically alters method function w/o need for subclasses
def anim(ug, simData, molnum, frameTracker): 
    print("in anim: test")
    print("in anim, before loop, molnum = ", molnum)
    f = mlab.gcf()
    #molnum = mayavi_widget.molecule_selected(mayavi_widget.molecule_selected)
    iterations = len(simData['trial0']['simulation']['dend']['times'])
    currentFrame = frameTracker.getCurrentFrame()
    while currentFrame < iterations:
        print("in anim molnum = ",molnum)
        concentrations = get_voxel_molecule_concs(currentFrame, simData, molnum)
        print("concentrations",concentrations)
        ug.point_data.scalars = np.repeat(concentrations, 8) #make 8 to static variable of "voxelpts"
        print("before surf----------------------____")                        
        surf = mlab.pipeline.surface(ug, opacity=0.1)            
        print("before pipeline----------------------____")
        mlab.pipeline.surface(mlab.pipeline.extract_edges(surf), color=(0, 0, 0))
        
        f.scene.render()
        print("currentFrame=",currentFrame)
        currentFrame += 1
        frameTracker.setCurrentFrame(currentFrame)
        yield 

      
def sendOutofWidget(molnum):
    print(molnum)
    #a = anim(create_morphology(simData), simData)
    
def get_molecule_list(simData):
    return simData['trial0']['output']['output_species']

def get_simData(fileName):
    simData = h5.File(fileName,"r")
    return simData

def molecule_to_number(molecule, simData):
    return np.where(simData['trial0']['output']['output_species'][:]==molecule)[0][0]    
    

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
    simData = get_simData(fileName)   
    moleculeList = get_molecule_list(simData) #simdata.... then past just the list below
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
    