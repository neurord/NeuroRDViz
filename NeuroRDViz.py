from __future__ import print_function, division
import os
os.environ['ETS_TOOLKIT'] = 'qt4'
os.environ['QT_API'] = 'pyqt'
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

globalnum = 0
################################################################################
#The actual visualization
class Visualization(HasTraits):
    scene = Instance(MlabSceneModel, ())
    
    @on_trait_change('scene.activated')
    def update_plot(self):
        # This function is called when the view is opened. We don't
        # populate the scene when the view is not yet open, as some
        # VTK features require a GLContext.     
                      
        surf = mlab.pipeline.surface(create_morphology(simData), opacity=0.1)
    
        self.scene.mlab.pipeline.surface(mlab.pipeline.extract_edges(surf), color=(0, 0, 0), )
        #anim(create_morphology(), 0)
        
    # the layout of the dialog created
    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=250, width=300, show_label=False),
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
    dendSnapshot = get_simData(get_fileName())['trial0']['simulation']['dend']['concentrations'][ms,:,molnum] #takes [milisecond, all voxel's of it's data, molecule of interest)
    somaSnapshot = get_simData(get_fileName())['trial0']['simulation']['soma']['concentrations'][ms,:,molnum] 
    wholeCellSnapshot = np.concatenate((dendSnapshot,somaSnapshot),axis=0)
    return wholeCellSnapshot

class MayaviQWidget(QtGui.QWidget):
    
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QVBoxLayout(self)
        layout.setContentsMargins(0,0,0,0)
        layout.setSpacing(0)
        self.visualization = Visualization()

        # If you want to debug, beware that you need to remove the Qt
        # input hook.
        #QtCore.pyqtRemoveInputHook()
        #import pdb ; pdb.set_trace()
        #QtCore.pyqtRestoreInputHook()

        # The edit_traits call will generate the widget to embed.
        self.ui = self.visualization.edit_traits(parent=self,
                                                 kind='subpanel').control
        layout.addWidget(self.ui)
        self.ui.setParent(self)
    
    def molecule_selected(self, text):
        simData = get_simData(get_fileName())
        globalnum = molecule_to_number(text, simData)
        #a = anim(create_morphology(simData), simData, molnum)
        print("in molecule_selected", text,globalnum)
        a = anim(create_morphology(simData), simData, 0)
        
@mlab.animate(delay=1000) #This is a "decorator" - dynamically alters method function w/o need for subclasses
def anim(ug, simData, molnum): 
    print("in anim: test")
    f = mlab.gcf()
    #molnum = mayavi_widget.molecule_selected(mayavi_widget.molecule_selected, )
    for i in range(len(simData['trial0']['simulation']['dend']['times'])):
        try:
            molnum = globalnum
            print("in anim molnum = ",molnum)
            print("in anim global = " ,globalnum)
            concentrations = get_voxel_molecule_concs(i, simData, 1)
            ug.point_data.scalars = np.repeat(concentrations, 8) #make 8 to static variable of "voxelpts"
            print("before surf----------------------____")                        
            surf = mlab.pipeline.surface(ug, opacity=0.1)            
            print("before pipeline----------------------____")
            mlab.pipeline.surface(mlab.pipeline.extract_edges(surf),
                            color=(0, 0, 0), )
            #print(concentrations)
            #f.scene.camera.azimuth(10)  #Rotates camera by 10
            f.scene.render()
        except AttributeError:
            pass
        if True != True: #set to something you want to make cause anim to break - assuming you find out how to have it read changes outside of the loop
                break;
        yield
    print("in anim, outside of loop, global = ", globalnum)
    
        #Increase Start/stop/delay window size
        #Adjust default - x10
        
    '''tip for accelerating animation:
        obj.scene.disable_render = True
        # Do all your scripting that takes ages.
        # ...
        # Once done, do the following:
        obj.scene.disable_render = False
        '''
        

      
def sendOutofWidget(morphology, molname, simData):
    print(molname, molecule_to_number(molname,simData))
    #a = anim(create_morphology(simData), simData)
    
def get_molecule_list(simData):
    return simData['trial0']['output']['output_species']

def get_simData(fileName):
    simData = h5.File(fileName,"r")
    return simData

def get_fileName():
    return sys.argv[1]

def molecule_to_number(molecule, simData):
    return np.where(simData['trial0']['output']['output_species'][:]==molecule)[0][0]    
    

if __name__ == "__main__":
    # Don't create a new QApplication, it would unhook the Events
    # set by Traits on the existing QApplication. Simply use the
    # '.instance()' method to retrieve the existing one.
    app = QtGui.QApplication.instance()
    '''GUI = Window()'''
    container = QtGui.QWidget()
    switch = False
    #container.setWindowTitle("NeuoRD Visualization")
    # define a "complex" layout to test the behaviour
    layout = QtGui.QGridLayout(container)
    comboBox = QtGui.QComboBox(container)
    fileName = sys.argv[1]
    simData = get_simData(fileName)   
    moleculeList = get_molecule_list(simData) #simdata.... then past just the list below
    for moleculeType in range(len(moleculeList)):
        comboBox.addItem(moleculeList[moleculeType])
    layout.addWidget(comboBox, 0, 0)  # 0,0 = top left widget location, 0,1 = one to the right of it, etc.
    label_list = [] 
    label_list.append(comboBox)
    # put some stuff around mayavi
    mayavi_widget = MayaviQWidget(container)
    comboBox.activated[str].connect(mayavi_widget.molecule_selected)

    

    layout.addWidget(mayavi_widget, 1, 1) # This is visualization of morphology
    container.show()
    window = QtGui.QMainWindow()
    window.setCentralWidget(container)
    window.show()
    print("after window shown")
    
    app.exec_() # Start the main event loop.