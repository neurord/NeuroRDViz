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
Avogadro=6.023e14
mol_per_nM_u3=Avogadro*1e-15

################################################################################
#The actual visualization
class Visualization(HasTraits):
    scene = Instance(MlabSceneModel, ())

    @on_trait_change('scene.activated')
    def update_plot(self):
        ug=create_morphology(simData)
        surf = mlab.pipeline.surface(ug, opacity=1)

        self.scene.mlab.pipeline.surface(mlab.pipeline.extract_edges(surf), color=(0, 0, 0))


    # the layout of the dialog created
    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=250, width=300, show_label=False),
                resizable=True  
                )



def create_morphology(simData):


    grid = np.array(getMorphologyGrid()).view(np.recarray)                          
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


    return ug


def get_voxel_molecule_conc(simData, moleculeType, out_location):
    grid_points=len(getMorphologyGrid())
    samples = out_location[moleculeType]['samples']
    #volume = model.grid[:][/vol
    #"snapshot" may now be the wrong name. Consider change:
    outputSet=np.zeros((samples,grid_points))
    for currentSet in out_location[moleculeType]['location'].keys():
        molnum = out_location[moleculeType]['location'][currentSet]['mol_index'] 
        voxels = out_location[moleculeType]['location'][currentSet]['elements'] 
        
        #NEEDS CHANGE: Must divide by avagadro's number to get concentration. and then divide by voxel size to get concentration
        tempSnapshot = simData['trial0']['output'][currentSet]['population'][:,:,molnum] #Check 1st one.
        print(np.shape(tempSnapshot), np.shape(outputSet[:,voxels]))
        outputSet[:,voxels]=tempSnapshot
    outputSetConcs = population_list_to_concentration_list(outputSet, simData['model']['grid']['volume'])
    return outputSetConcs


#Conert molecular population to molecular concentration
def population_list_to_concentration_list(pop_list, voxel_volumes):
    conc_list = np.zeros((len(pop_list),len(pop_list[0])))
    
    #Iterate through one timeframe of pop_list to divide each voxel's population by the ~[grid][voxel volume]
    for z, pop_list_snapshot in enumerate(pop_list):
        for i, (a,b) in enumerate(zip(pop_list_snapshot, voxel_volumes)):
            conc_list[z][i] = (a/b) * mol_per_nM_u3
    return conc_list


class MayaviQWidget(QtGui.QWidget):
    animator = None
    currentFrame = 0
    f = mlab.gcf()
    #unable to call functions within this space, sending methods to anim.
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QVBoxLayout(self)
        layout.setContentsMargins(0,0,0,0)
        layout.setSpacing(0)
        self.visualization = Visualization()

        # The edit_traits call will generate the widget to embed.
        self.ui = self.visualization.edit_traits(parent=self, kind='subpanel').control
        layout.addWidget(self.ui)
        self.ui.setParent(self)
        
    
    def molecule_selected(self, text):
        #This is necessary otherwise ~anim-loop-obect will be instantiated initially 
        #with whatever conditions passed
        if text != None: 
            if self.animator != None:
                self.animator.close()
            self.animator = anim(create_morphology(simData), simData, text, self, self.f)

    def setCurrentFrame(self, frame):
        self.currentFrame = frame

    def getCurrentFrame(self):
        return self.currentFrame
        
    

@mlab.animate(delay=10) 
def anim(ug, simData, moleculeType, frameTracker, f):
    
    out_location,dt,samples = get_mol_info(simData,simData['model']['output']['__main__']['species'][:],getMorphologyGrid())
    molnum = get_mol_index(simData, "all", moleculeType)
    iterations = len(simData['trial0']['output']['all']['times'])  #times[1] - times[0] = dt                                         #change! all to chosen outputSets from getSomething
    currentFrame = frameTracker.getCurrentFrame()
    population = get_voxel_molecule_conc(simData, moleculeType, out_location)
    print(np.shape(population))
    surf = mlab.pipeline.surface(ug, opacity =1, colormap='hot')  # Decide how max/min color values are assigned.
    mlab.pipeline.surface(mlab.pipeline.extract_edges(surf), color=(0, 0, 0)) 
    #mlab.colorbar(title='Concentration', orientation='vertical', nb_labels=7)
   

    while currentFrame < iterations:
        mlab.colorbar(title='Concentration', orientation='vertical', nb_labels=7)
        concentrations = population[currentFrame,:]
        ug.point_data.scalars = np.repeat(concentrations, 8) 
        ug.point_data.scalars.name = 'concentrations' 
        ug.modified()

        currentFrame += 1
        #getQtWindow().label.setText(currentFrame+ "ms")
        print("Progress:",currentFrame)
        frameTracker.setCurrentFrame(currentFrame)
        yield


def getMoleculeList(simData):
    return simData['model']['species']
                                            
def get_h5simData(fileName):
    simData = h5.File(fileName,"r")
    return simData
    
def getMorphologyGrid():
    return simData['model']['grid']
    
def getQtWindow():
    return container 
    
#Searches the list of molecules and returns thes corresponding index. Returns -1 if not found.
def get_mol_index(simData, outputSet, molecule):
    indices=np.where(simData['model']['output'][outputSet]['species'][:]==molecule)[0]
    if len(indices) == 1:
        return indices[0]
    else:
        return -1

#Returns the largest dataset - should contain all species. (breaking on "TypeError:
#string indices must be integers. "Picking battles" and going with always picking 
#__main__ for now.
'''def get_Master_Set(simData):
    biggestOutputSet = []
    for thisOutputSet in simData['model']['output']:
        print(len(thisOutputSet['species'][:]))
        print(len(biggestOutputSet))
        if len(thisOutputSet['species']) > len(biggestOutputSet):
            biggestOutputSet = thisOutputSet['species'][:]
    return biggestOutputSet'''
 
def get_mol_info(simData,plot_molecules,grid_points):
    outputsets=simData['model']['output'].keys()
    dt=np.zeros((len(plot_molecules)))
    samples=np.zeros((len(plot_molecules)),dtype=int)
    out_location={}
    for imol,molecule in enumerate(plot_molecules):
        temp_dict={}
        tot_voxels=0
        for outset in outputsets[1:]:           #better to go backward from last set, and then go to 0 set if mol not found
        #for each in outputsets[-1::-1]:
            mol_index=get_mol_index(simData,outset, molecule)
            if mol_index>-1:
                samples[imol]=len(simData['trial0']['output'][outset]['times'])
                dt[imol]=simData['trial0']['output'][outset]['times'][1]/1000. #convert msec to sec
                tot_voxels=tot_voxels+len(simData['model']['output'][outset]['elements'])
                temp_dict[outset]={'mol_index':mol_index,'elements':simData['model']['output'][outset]['elements'][:]}
        if len(temp_dict)>0:
            out_location[molecule]={'samples':samples[imol],'dt':dt[imol],'voxels': tot_voxels,'location': temp_dict}
        else:
            outset=outputsets[0]
            print("************* MOLECULE",molecule, " NOT IN REGULAR OUTPUT SETS !")
            mol_index=get_mol_index(simData,outset,molecule)
            samples[imol]=len(simData['trial0']['output'][outset]['times'])
            dt[imol]=simData['trial0']['output'][outset]['times'][1]/1000. #convert msec to sec
            temp_dict[outset]={'mol_index':mol_index,'elements':simData['model']['output'][outset]['elements'][:]}
            out_location[molecule]={'samples':samples[imol],'dt':dt[imol],'voxels': grid_points,'location': temp_dict}
    return out_location,dt,samples

def make_temp_dict(simData):
    samples[moleculeType]=len(simData['trial0']['output'][outset]['times'])
    dt[moleculeType]=simData['trial0']['output'][outset]['times'][1]/1000. #convert msec to sec
    temp_dict[outset]={'mol_index':mol_index,'elements':simData['model']['output'][outset]['elements'][:]}
    return temp_dict

 
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
    label = QtGui.QLabel(container)
    label.setText("Progress: ms")
    label.setGeometry(100,100, 100, 100)
    label.setAlignment(QtCore.Qt.AlignHCenter|QtCore.Qt.AlignVCenter)
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
    app.exec_() # Start the main event loop.
