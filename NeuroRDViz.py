#Developed by Bradley William English - brad.w.english@gmail.com

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

#Combobox that allows for search/autocompletion
class ExtendedCombo( QComboBox ):
    def __init__( self,  parent = None):
        super( ExtendedCombo, self ).__init__( parent )

        self.setFocusPolicy( Qt.StrongFocus )
        self.setEditable( True )
        self.completer = QCompleter(self)

        # always show all completions
        self.completer.setCompletionMode( QCompleter.UnfilteredPopupCompletion )
        self.pFilterModel = QSortFilterProxyModel( self )
        self.pFilterModel.setFilterCaseSensitivity( Qt.CaseInsensitive )
        self.completer.setPopup( self.view() )

        self.setCompleter( self.completer )

        self.lineEdit().textEdited[unicode].connect( self.pFilterModel.setFilterFixedString )
        self.completer.activated.connect(self.setTextIfCompleterIsClicked)

    def setModel( self, model ):
        super(ExtendedCombo, self).setModel( model )
        self.pFilterModel.setSourceModel( model )
        self.completer.setModel(self.pFilterModel)

    def setModelColumn( self, column ):
        self.completer.setCompletionColumn( column )
        self.pFilterModel.setFilterKeyColumn( column )
        super(ExtendedCombo, self).setModelColumn( column )


    def view( self ):
        return self.completer.popup()

    def index( self ):
        return self.currentIndex()

    def setTextIfCompleterIsClicked(self, text):
        if text:
            index = self.findText(text)
            self.setCurrentIndex(index)
        
class Visualization(HasTraits):
    scene = Instance(MlabSceneModel, ())

    @on_trait_change('scene.activated')
    def update_plot(self):
        ug=create_morphology(simData)
        surf = mlab.pipeline.surface(ug, opacity=1)
        self.scene.mlab.pipeline.surface(mlab.pipeline.extract_edges(surf), color=(0, 0, 0)) # @UndefinedVariable - this comments tells Eclipse IDE to ignore "error"
        mlab.axes(surf, nb_labels=7)

    # the layout of the dialog created
    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=250, width=300, show_label=False),
                resizable=True  
                )

def create_morphology(simData):
           
    grid = np.array(getMorphologyGrid()).view(np.recarray)                          
    points = np.array((
         (grid.x0, grid.y0, grid.z0), (grid.x1, grid.y1, grid.z1), (grid.x2, grid.y2, grid.z2), (grid.x3, grid.y3, grid.z3-grid.deltaZ),
         (grid.x0, grid.y0, grid.z0+grid.deltaZ), (grid.x1, grid.y1, grid.z1+grid.deltaZ), (grid.x2, grid.y2, grid.z2+grid.deltaZ), (grid.x3, grid.y3, grid.z3),))
    points = points.swapaxes(0, 2).swapaxes(1, 2)
    points = points.reshape(-1, 3)
    voxels = np.arange(points.shape[0]).reshape(-1, 8)

    voxel_type = tvtk.Hexahedron().cell_type # @UndefinedVariable - this comments tells Eclipse IDE to ignore "error"
    ug = tvtk.UnstructuredGrid(points=points) # @UndefinedVariable - this comments tells Eclipse IDE to ignore "error"
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

class colorBarInputDialog(QWidget):

    def __init__(self):
        super(colorBarInputDialog, self).__init__()
        self.setGeometry(0,50,200,50)
        layout = QFormLayout()
        
        self.minLabel = QLabel(str(mayavi_widget_list[0].colorbar_min))
        self.btnMin = QPushButton("Min:")
        self.btnMin.clicked.connect(self.getMin)
        layout.addRow(self.btnMin,self.minLabel)
        
        self.maxLabel = QLabel(str(mayavi_widget_list[0].colorbar_max))
        self.btnMax = QPushButton("Max:")
        self.btnMax.clicked.connect(self.getMax)
        layout.addRow(self.btnMax,self.maxLabel)
        
        self.btnDefeault = QPushButton("Restore Defaults")
        self.btnDefeault.clicked.connect(self.restoreDefaults)
        layout.addRow(self.btnDefeault)
        
        self.btnScale = QPushButton("Scale:")
        self.btnScale.clicked.connect(self.getItem)
        self.leScale = QLabel()
        self.leScale.setText(str("Linear"))
        layout.addRow(self.btnScale,self.leScale)
        
        self.btnApply = QPushButton("Apply")
        self.btnClose = QPushButton("Close")
        self.btnApply.clicked.connect(self.applyChanges)
        self.btnClose.clicked.connect(self.closePopup)
        layout.addRow(self.btnClose, self.btnApply)
        
        self.setLayout(layout)
        self.setWindowTitle("Colorbar Options")
        
    def getItem(self): 
        items = ("Linear", "Logarithmic (Note: Min cannot be 0!)")
        item, ok = QInputDialog.getItem(self, "select input dialog", "Select Scale:", items, 0, False)
        if ok and item:
            self.leScale.setText(item)
            
    def gettext(self): #Not Used at the moment
        text, ok = QInputDialog.getText(self, 'Text Input Dialog', 'Enter your name:')
        if ok:
            self.le1.setText(str(text))
            
    def getMin(self):
        num,ok = QInputDialog.getText(self,"Double Input Dialog","Enter Min")
        if ok:
            self.minLabel.setText(str(num))
            
    def getMax(self):
        num,ok = QInputDialog.getText(self,"Double Input Dialog","Enter Max")
        if ok:
            self.maxLabel.setText(str(num))
    def restoreDefaults(self):
        try:
            mayavi_widget_list[0].colorBarDummySurf.module_manager.scalar_lut_manager.data_range = [0, np.max(mayavi_widget_list[0].population)]
            mayavi_widget_list[0].surf.module_manager.scalar_lut_manager.data_range = [0, np.max(mayavi_widget_list[0].population)]
            self.minLabel.setText(str(0))
            self.maxLabel.setText(str(np.max(mayavi_widget_list[0].population)))
        except:
            self.msg = QMessageBox()
            self.msg.setIcon(QMessageBox.Information)
            self.msg.setText("Select a Molecule First")
            self.msg.setWindowTitle("No Existing Defaults")
    
    def applyChanges(self):
        try:    
            newMin, newMax = float(self.minLabel.text()), float(self.maxLabel.text()) 
            mayavi_widget_list[0].colorBarDummySurf.module_manager.scalar_lut_manager.data_range = [newMin, newMax]            
            mayavi_widget_list[0].surf.module_manager.scalar_lut_manager.data_range = [newMin, newMax]
            lut = mayavi_widget_list[0].surf.module_manager.scalar_lut_manager.lut
            cbarlut = mayavi_widget_list[0].colorBarDummySurf.module_manager.scalar_lut_manager.lut
            if self.leScale.text() == "Logarithmic (Note: Min cannot be 0!)":
                lut.scale = 'log10'
                cbarlut.scale = 'log10'
            elif self.leScale.text() == "Linear":
                lut.scale = 'linear'
                cbarlut.scale = 'linear'
        except:
            self.msg = QMessageBox()
            self.msg.setIcon(QMessageBox.Information)
            self.msg.setText("Select a Molecule First")
            self.msg.setWindowTitle("No Existing Colorbar")   
    def closePopup(self):
        self.close()

class helpWindow(QWidget):

    def __init__(self):
        super(helpWindow, self).__init__()
        self.setGeometry(0,50,100,50)
        layout = QFormLayout()
        
        self.minLabel = QLabel("Coming Soon!")
        layout.addRow(self.minLabel)
        
        self.setLayout(layout)
        self.setWindowTitle("Help Menu")
        
    def closePopup(self):
        self.close()
        
class Window(QtGui.QMainWindow):
    
    #Main Menus should go here... things which need to appear regardless
    def __init__(self):
        super(Window, self).__init__()
        self.setGeometry(50, 50, 1100, 800)
        self.setWindowTitle("NeuoRD Visualizer")
        
        exitAction = QtGui.QAction("&Exit -", self)
        exitAction.setShortcut("Ctrl+Q")
        exitAction.setStatusTip('Close Application')
        exitAction.triggered.connect(self.close_application) #.triggered = .clicked
        
        minMaxColorBarAction = QtGui.QAction("&Min/Max Range for Colorbar -", self)
        minMaxColorBarAction.setShortcut("Ctrl+M")
        minMaxColorBarAction.setStatusTip('Change the displayed minimum & maximum ranges on the color bar.')
        minMaxColorBarAction.triggered.connect(self.changeMinMaxColorBar)
        
        addAction = QtGui.QAction("&Add a Model -", self)
        addAction.setShortcut("Ctrl+A")
        addAction.setStatusTip('Add Items to Visualizer')
        addAction.triggered.connect(self.add_model) #.triggered = .clicked
        
        helpAction = QtGui.QAction("&How Do I... -", self)
        helpAction.setShortcut("Ctrl+H")
        helpAction.setStatusTip('Learn More About How to Use the Visualizer')
        helpAction.triggered.connect(self.help_action) #.triggered = .clicked
        
        mainMenu = self.menuBar()
        fileMenu = mainMenu.addMenu('&File')
        fileMenu.addAction(exitAction)
        editMenu = mainMenu.addMenu('&Edit')
        editMenu.addAction(minMaxColorBarAction)
        addMenu = mainMenu.addMenu('&Add')
        addMenu.addAction(addAction)
        helpMenu = mainMenu.addMenu('&Help')
        helpMenu.addAction(helpAction)
        
        self.modelCount = 1
        self.currentRow = 4
        self.currentColumn = 1
        
        self.statusBar()
        self.home()
        
    def home(self):
       
        toolBarColorBarMinMax = QtGui.QAction(QtGui.QIcon('colorBarIcon.png'), "Set Min/Max of ColorBar", self)
        toolBarColorBarMinMax.setStatusTip('Change the default min/max range on the color bar.')
        toolBarColorBarMinMax.triggered.connect(self.changeMinMaxColorBar)
        toolBarAddModel = QtGui.QAction(QtGui.QIcon('addModelIcon.png'), "Add a Model", self)
        toolBarAddModel.setStatusTip('Add another model to the visualizer.')
        toolBarAddModel.triggered.connect(self.add_model)
        toolBarHelp = QtGui.QAction(QtGui.QIcon('helpIcon.png'), "Help", self)
        toolBarHelp.setStatusTip("Learn More About How to Use the Visualizer")
        toolBarHelp.triggered.connect(self.help_action)
        
        self.toolBar = self.addToolBar("ToolBar")
        self.toolBar.addAction(toolBarColorBarMinMax)
        self.toolBar.addAction(toolBarAddModel)
        self.toolBar.addAction(toolBarHelp)

        self.show()

    def close_application(self):
        #print("Brad: Put are you sure? back in")
        #Remove following line:
        #sys.exit()
        #Add back in lines below:
        choice = QtGui.QMessageBox.question(self, 'Exit', "Are you sure?", "Yes", "No")
        if choice == 0:
            sys.exit()
    
    def changeMinMaxColorBar(self):
        self.newEditWindow = colorBarInputDialog()
        self.newEditWindow.show()
    
    def help_action(self):
        self.newHelpWindow = helpWindow()
        self.newHelpWindow.show()
    
    def add_model(self): 
        mayavi_widget_list.append(MayaviQWidget(container, progress_bar, progress_slider, progress_label, progress_slider_label))
        if self.modelCount % 2 != 0:
            layout.addWidget(mayavi_widget_list[self.modelCount], self.currentRow, self.currentColumn+1)
            self.currentRow += 1
        else:
            layout.addWidget(mayavi_widget_list[self.modelCount], self.currentRow, self.currentColumn)
        self.modelCount += 1
    def resetAnimation(self):
        mayavi_widget_list[0].setCurrentFrame(0)
        mayavi_widget_list[0].progress_slider.setValue(0)
        mayavi_widget_list[0].progress_bar.setValue(0)
        
        
class MayaviQWidget(QtGui.QWidget):
    #unable to call functions within this space, sending methods to anim
    def __init__(self, parent, progress_bar, progress_slider, 
                 progress_label, progress_slider_label):
        
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
        self.colorbar_ug = tvtk.UnstructuredGrid() # @UndefinedVariable
        self.colorbar_min, self.colorbar_max = 0,0 
        self.colorbar_ug.point_data.scalars = np.linspace(self.colorbar_min, self.colorbar_max,7)  
        self.colorbar_ug.point_data.scalars.name = 'concentrations'
        self.colorBarDummySurf = mlab.pipeline.surface(self.colorbar_ug, opacity =0.8, colormap='hot')
        self.colorBar = mlab.colorbar(object=self.colorBarDummySurf, title='Concentration', orientation='vertical')
        self.iterations = 0
        self.surf = None
        self.progress_bar = progress_bar
        self.progress_slider = progress_slider
        self.progress_slider_label = progress_slider_label
        self.progress_label = progress_label
        self.population = 0
        
    
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
        self.progress_slider_label.setText(str(x/1000)+ "s")

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
    #Note to Self: Can probably change every instance of "mayavi_widget_list[0]" to "mayavi_widget_list[0]" from __main__
        
    #The decorated function will return the Animator instance used and a user may call its Stop method to stop the animation.
    #Simulation Data Gathering
    out_location,dt,samples = get_mol_info(simData,simData
                                           ['model']['output']['__main__']['species'][:],getMorphologyGrid())
    molnum = get_mol_index(simData, "all", moleculeType)
    mayavi_widget_list[0].population = get_voxel_molecule_conc(simData, moleculeType, out_location)
    mayavi_widget_list[0].iterations = out_location[moleculeType]['samples']
    dt=out_location[moleculeType]['dt']
    #Creates mayavi surface to be shown corresponding with the unstructured grid(ug)
    mayavi_widget_list[0].surf = mlab.pipeline.surface(ug, opacity =1, colormap='hot') 
    mlab.pipeline.surface(mlab.pipeline.extract_edges(mayavi_widget_list[0].surf), color=(0, 0, 0))# @UndefinedVariable
    mayavi_widget_list[0].surf.module_manager.scalar_lut_manager.data_range = [0, np.max(mayavi_widget_list[0].population)]
    
    #Set Colorbar range for this Molecule Type
    mayavi_widget_list[0].colorbar_min, mayavi_widget_list[0].colorbar_max = 0, np.max(mayavi_widget_list[0].population)
    mayavi_widget_list[0].colorBarDummySurf.module_manager.scalar_lut_manager.data_range = [mayavi_widget_list[0].colorbar_min, mayavi_widget_list[0].colorbar_max]
    

    #Actual Animation Loop
    while mayavi_widget_list[0].getCurrentFrame() < mayavi_widget_list[0].iterations:       
        
        concentrations = mayavi_widget_list[0].population[mayavi_widget_list[0].getCurrentFrame(),:]
        ug.point_data.scalars = np.repeat(concentrations, 8)  # Decide how max/min color values are assigned.
        ug.point_data.scalars.name = 'concentrations' 
        ug.modified()
       
        mayavi_widget_list[0].setCurrentFrame(mayavi_widget_list[0].getCurrentFrame()+1)   
        mayavi_widget_list[0].progress_label.setText(str(mayavi_widget_list[0].getCurrentFrame()/1000) + "s")
        mayavi_widget_list[0].progress_bar.setValue((mayavi_widget_list[0].getCurrentFrame()/mayavi_widget_list[0].iterations)*100)
        yield
        
    #If completed, reset to beginning.
    if mayavi_widget_list[0].getCurrentFrame() >= (mayavi_widget_list[0].iterations-1):
        mayavi_widget_list[0].setCurrentFrame(0)


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
    outputsets=simData['model']['output'].keys() #Gathers list of outputsets
    dt=np.zeros((len(plot_molecules)))
    samples=np.zeros((len(plot_molecules)),dtype=int)
    out_location={}
    for imol,molecule in enumerate(plot_molecules):
        temp_dict={}
        tot_voxels=0
        for outset in outputsets[1:]: #better to go backward from last set, and then go to 0 set if mol not found
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
        fileName=fname # @UndefinedVariable
    except NameError:
        fileName = sys.argv[1] 
    
    simData = get_h5simData(fileName)
    
    progress_bar = QtGui.QProgressBar()

    progress_slider = QSlider(Qt.Horizontal)
    
    
    app = QtGui.QApplication.instance()
    container = QtGui.QWidget()
    layout = QtGui.QGridLayout(container)
    model = QStandardItemModel() #Required for searchable comboBox
    progress_label = QtGui.QLabel(container)
    progress_slider_label = QtGui.QLabel(container)
    window = Window()
    comboBox = ExtendedCombo()
    mayavi_widget_list = []
    reset_button_list =[]

    mayavi_widget_list.append(MayaviQWidget(container, progress_bar, progress_slider, progress_label, progress_slider_label)) #Initiates first morphology
    reset_button_list.append(QtGui.QPushButton('Reset', window))
    reset_button_list[0].clicked.connect(window.resetAnimation)
    
    animControlWindow_list = []
    
    #mayavi_widget_list[0]3 = MayaviQWidget(container, progress_bar, progress_slider, progress_label, progress_slider_label)
    #mayavi_widget_list[0]4 = MayaviQWidget(container, progress_bar, progress_slider, progress_label, progress_slider_label)
    
    #These go somewhere after layout is created:
    #layout.addWidget(mayavi_widget_list[0]3, 5, 1)
    #layout.addWidget(mayavi_widget_list[0]4, 5, 2)

    
    moleculeList = sorted(getMoleculeList(simData)) #simdata.... then past just the list below
    for i,moleculeType in enumerate(moleculeList):
        item = QStandardItem(moleculeType)
        model.setItem(i, 0, item)
    
    
    comboBox.setModel(model)
    comboBox.setModelColumn(0)
    comboBox.activated[str].connect(mayavi_widget_list[0].molecule_selected)
    
    progress_slider.valueChanged.connect(mayavi_widget_list[0].slider_movement)
    
    fileNameLabel = QtGui.QLabel(mayavi_widget_list[0])
    fileNameLabel.setText(fileName)
    
    layout.addWidget(comboBox, 0, 0)  # 0,0 = top left widget location, 0,1 = one to the right of it, etc.
    layout.addWidget(fileNameLabel, 0,1)
    layout.addWidget(mayavi_widget_list[0], 4, 1) # Visualization of morphology
    layout.addWidget(reset_button_list[0], 5, 1)
    layout.addWidget(progress_label, 2,0)
    layout.addWidget(progress_bar, 2, 1)
    layout.addWidget(progress_slider_label,3, 0)  
    layout.addWidget(progress_slider, 3, 1)
    

    container.show()
   
    
    window.setCentralWidget(container)
    window.show()
    app.exec_() # Start the main event loop.
