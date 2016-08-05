#Bradley William English - brad.w.english@gmail.com

from __future__ import print_function, division
import os
os.environ['ETS_TOOLKIT'] = 'qt4'
os.environ['QT_API'] = 'pyqt'
import sip
sip.setapi('Qstring',2)
from pyface.qt import QtGui, QtCore
from PyQt4.QtCore import  *
from PyQt4.QtGui import *

from traits.api import HasTraits, Instance, on_trait_change
from traitsui.api import View, Item
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, \
        SceneEditor

import numpy as np
from tvtk.api import tvtk
from mayavi import mlab     ######## This line of code takes MUCH longer than anything else (be more selective in imports to expedite?)

import h5py as h5
import sys
Avogadro=6.023e14
mol_per_nM_u3=Avogadro*1e-15

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
    outputSet=np.zeros((samples,grid_points))
    for currentSet in out_location[moleculeType]['location'].keys():
        molnum = out_location[moleculeType]['location'][currentSet]['mol_index'] 
        voxels = out_location[moleculeType]['location'][currentSet]['elements'] 
        tempSnapshot = simData['trial0']['output'][currentSet]['population'][:,:,molnum] 
        outputSet[:,voxels]=tempSnapshot
    outputSetConcs = population_to_concentration(outputSet, simData['model']['grid']['volume'])
    return outputSetConcs


#Conert molecular population to molecular concentration
def population_to_concentration(pop_list, voxel_volumes):
    conc_list = np.zeros(np.shape(pop_list))
    
    #Iterate through one timeframe of pop_list to divide each voxel's population by the ~[grid][voxel volume]
    for z, pop_snapshot in enumerate(pop_list): #
        for i, (a,b) in enumerate(zip(pop_snapshot, voxel_volumes)):
            conc_list[z][i] = (a/b) / mol_per_nM_u3
            
    #print out a volumes[3] volumes[6] & conc_list[3] [6] with two types 
    return conc_list


class MayaviQWidget(QtGui.QWidget):

    #unable to call functions within this space, sending methods to anim.
    def __init__(self, parent, progress_bar, progress_slider, progress_label, progress_slider_label):
        QtGui.QWidget.__init__(self, parent)
        layout = QtGui.QVBoxLayout(self)
        layout.setContentsMargins(0,0,0,0)
        layout.setSpacing(0)
        self.visualization = Visualization()
        # The edit_traits call will generate the widget to embed.
        self.ui = self.visualization.edit_traits(parent=self, kind='subpanel').control
        layout.addWidget(self.ui)
        self.ui.setParent(self)
        self.home()
    
    def home(self):
        self.animator = None
        self.currentFrame = 0
        self.colorbar_ug = tvtk.UnstructuredGrid()
        self.colorbar_min, self.colorbar_max = 0,0 
        self.colorbar_ug.point_data.scalars = np.linspace(self.colorbar_min, self.colorbar_max,7)  
        self.colorbar_ug.point_data.scalars.name = 'concentrations'
        self.colorBarDummySurf = mlab.pipeline.surface(self.colorbar_ug, opacity =1, colormap='hot')
        self.colorBar = mlab.colorbar(object=self.colorBarDummySurf, title='Concentration', orientation='vertical', nb_labels=7)
        self.colorBar.visible = False
        self.iterations = 0
        self.progress_bar = progress_bar
        self.progress_slider = progress_slider
        self.progress_slider_label = progress_slider_label
        self.progress_label = progress_label
    
    def molecule_selected(self, text):
        #This is necessary otherwise ~anim-loop-obect will be instantiated initially 
        #with whatever conditions passed
        if text != None:
            if self.animator != None:
                self.animator.close()
            self.animator = anim(create_morphology(simData), simData, text, self)
    
    def slider_movement(self):
        position = self.progress_slider.value()
        x = int((position/100)*self.iterations)
        self.setCurrentFrame(x)
        self.progress_slider_label.setText(str(x)+ "ms")

    def setCurrentFrame(self, frame):
        self.currentFrame = frame

    def getCurrentFrame(self):
        return self.currentFrame
    
    def setcolorbar_min(self, min):
        self.colorbar_min = min

    def getcolorbar_min(self):
        return self.colorbar_min
    
    def setcolorbar_max(self, max):
        self.colorbar_max = max

    def getcolorbar_max(self):
        return self.colorbar_max

    
@mlab.animate(delay=10) 
def anim(ug, simData, moleculeType, widgetObject):

    #Simulation Data Gathering
    out_location,dt,samples = get_mol_info(simData,simData['model']['output']['__main__']['species'][:],getMorphologyGrid())
    molnum = get_mol_index(simData, "all", moleculeType)
    population = get_voxel_molecule_conc(simData, moleculeType, out_location)
    widgetObject.iterations = out_location[moleculeType]['samples']
    dt=out_location[moleculeType]['dt']
    
    #Creates mayavi surface to be shown corresponding with the unstructured grid(ug)
    surf = mlab.pipeline.surface(ug, opacity =1, colormap='hot') 
    mlab.pipeline.surface(mlab.pipeline.extract_edges(surf), color=(0, 0, 0)) 
    surf.module_manager.scalar_lut_manager.data_range = [0, np.max(population)]
    
    #Set Colorbar range for this Molecule Type
    widgetObject.colorBarDummySurf.module_manager.scalar_lut_manager.data_range = [0, np.max(population)]
    widgetObject.colorBar.visible = True

    #Actual Animation Loop
    while widgetObject.getCurrentFrame() < widgetObject.iterations:
        concentrations = population[widgetObject.getCurrentFrame(),:]
        ug.point_data.scalars = np.repeat(concentrations, 8)  # Decide how max/min color values are assigned.
        ug.point_data.scalars.name = 'concentrations' 
        ug.modified()
       
        widgetObject.setCurrentFrame(widgetObject.getCurrentFrame()+1)   
        widgetObject.progress_label.setText(str(widgetObject.getCurrentFrame()) + "ms")
        widgetObject.progress_bar.setValue((widgetObject.getCurrentFrame()/widgetObject.iterations)*100)
        #widgetObject.progress_slider.setValue((widgetObject.getCurrentFrame()/iterations)*100)
        yield
        
    #If completed, reset to beginning.
    if widgetObject.getCurrentFrame() >= (iterations-1):
        widgetObject.setCurrentFrame(0)


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

if __name__ == "__main__":
    
    try:
        fileName=fname
    except NameError:
        fileName = sys.argv[1] 
    
    simData = get_h5simData(fileName)
    
    progress_bar = QtGui.QProgressBar()

    progress_slider = QSlider(Qt.Horizontal)
    
    
    app = QtGui.QApplication.instance()
    container = QtGui.QWidget()
    layout = QtGui.QGridLayout(container)
    comboBox = QtGui.QComboBox(container)
    progress_label = QtGui.QLabel(container)
    progress_slider_label = QtGui.QLabel(container)
                                
    mayavi_widget = MayaviQWidget(container, progress_bar, progress_slider, progress_label, progress_slider_label)
    
    
    moleculeList = getMoleculeList(simData) #simdata.... then past just the list below
    for moleculeType in range(len(moleculeList)):
        comboBox.addItem(moleculeList[moleculeType])
    
    comboBox.activated[str].connect(mayavi_widget.molecule_selected)
    
    progress_slider.valueChanged.connect(mayavi_widget.slider_movement)
    
    fileNameLabel = QtGui.QLabel(mayavi_widget)
    fileNameLabel.setText(fileName)
    
    
    layout.addWidget(comboBox, 0, 0)  # 0,0 = top left widget location, 0,1 = one to the right of it, etc.
    layout.addWidget(mayavi_widget, 1, 1) # This is visualization of morphology
    layout.addWidget(fileNameLabel, 0,1)
    layout.addWidget(progress_slider, 3, 1)
    layout.addWidget(progress_bar, 2, 1)
    layout.addWidget(progress_label, 2,0)
    layout.addWidget(progress_slider_label,3, 0)  
    container.show()
   
    app.exec_() # Start the main event loop.
