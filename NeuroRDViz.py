'''
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version and with attribution to the author.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    Developed by Bradley William English - brad.w.english@gmail.com
'''
from __future__ import annotations

import os
import sys

os.environ.setdefault("ETS_TOOLKIT", "qt")
os.environ.setdefault("QT_API", "pyside6")

from PySide6 import QtCore, QtGui, QtWidgets
from PySide6.QtCore import Qt

# ----------------------------
# Qt4-style compatibility aliases
# ----------------------------

# Widgets
QWidget = QtWidgets.QWidget
QMainWindow = QtWidgets.QMainWindow
QApplication = QtWidgets.QApplication
QComboBox = QtWidgets.QComboBox
QCompleter = QtWidgets.QCompleter
QFormLayout = QtWidgets.QFormLayout
QLabel = QtWidgets.QLabel
QPushButton = QtWidgets.QPushButton
QInputDialog = QtWidgets.QInputDialog
QMessageBox = QtWidgets.QMessageBox
QProgressBar = QtWidgets.QProgressBar
QSlider = QtWidgets.QSlider
QVBoxLayout = QtWidgets.QVBoxLayout
QGridLayout = QtWidgets.QGridLayout

# Models / data
QSortFilterProxyModel = QtCore.QSortFilterProxyModel
QStandardItemModel = QtGui.QStandardItemModel
QStandardItem = QtGui.QStandardItem

# Actions / icons
QAction = QtGui.QAction
QIcon = QtGui.QIcon

# ---- Map QtGui.* widget references used in legacy code ----
QtGui.QWidget = QtWidgets.QWidget
QtGui.QMainWindow = QtWidgets.QMainWindow
QtGui.QApplication = QtWidgets.QApplication
QtGui.QVBoxLayout = QtWidgets.QVBoxLayout
QtGui.QGridLayout = QtWidgets.QGridLayout
QtGui.QFormLayout = QtWidgets.QFormLayout
QtGui.QLabel = QtWidgets.QLabel
QtGui.QPushButton = QtWidgets.QPushButton
QtGui.QInputDialog = QtWidgets.QInputDialog
QtGui.QMessageBox = QtWidgets.QMessageBox
QtGui.QProgressBar = QtWidgets.QProgressBar
QtGui.QSlider = QtWidgets.QSlider


os.environ["QT_API"] = "pyside6"

from pyface.qt import QtCore as _QtCore  # noqa
from pyface.qt import QtGui as _QtGui    # noqa
from pyface.qt import QtWidgets as _QtWidgets  # noqa

from traits.api import HasTraits, Instance, on_trait_change
from traitsui.api import View, Item
from mayavi.core.ui.api import MayaviScene, MlabSceneModel, SceneEditor
from tvtk.api import tvtk
from mayavi import mlab

import numpy as np
import h5py as h5

Avogadro=6.023e14
mol_per_nM_u3=Avogadro*1e-15

'''
This class allows for a Combobox that auto-completes as you enter molecule types.
Likely never needs to be changed.
'''
class ExtendedCombo(QtWidgets.QComboBox):

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

        self.lineEdit().textEdited.connect(self.pFilterModel.setFilterFixedString)

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

'''
This class is responsible for creating the embedded window of Mayavi that you display a model in.
See this link for more: http://docs.enthought.com/mayavi/mayavi/building_applications.html
'''
class Visualization(HasTraits):
    scene = Instance(MlabSceneModel, ())

    @on_trait_change('scene.activated')
    def update_plot(self):
        print('simData:', simData)
        ug = create_morphology(simData)
        surf = mlab.pipeline.surface(ug, opacity=1)
        self.scene.mlab.pipeline.surface(mlab.pipeline.extract_edges(surf), color=(0, 0, 0)) # @UndefinedVariable - this comment tells Eclipse IDE to ignore "error"
        mlab.axes(surf, nb_labels=7)

    # the layout of the dialog created
    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=250, width=300, show_label=False),
                resizable=True
                )
'''
This function maps the dendritic voxels to
8 respective points on a mayavi will group into Hexahedrons
'''
def create_morphology(simData):

    grid = np.array(getMorphologyGrid()).view(np.recarray)
    points = np.array((
         (grid.x0, grid.y0, grid.z0), (grid.x1, grid.y1, grid.z1), (grid.x2, grid.y2, grid.z2), (grid.x3, grid.y3, grid.z3-grid.deltaZ),
         (grid.x0, grid.y0, grid.z0+grid.deltaZ), (grid.x1, grid.y1, grid.z1+grid.deltaZ), (grid.x2, grid.y2, grid.z2+grid.deltaZ),
         (grid.x3, grid.y3, grid.z3),))
    points = points.swapaxes(0, 2).swapaxes(1, 2)
    points = points.reshape(-1, 3)
    voxels = np.arange(points.shape[0]).reshape(-1, 8)

    voxel_type = tvtk.Hexahedron().cell_type # @UndefinedVariable - this comment tells Eclipse IDE to ignore "error"
    ug = tvtk.UnstructuredGrid(points=points) # @UndefinedVariable - this comment tells Eclipse IDE to ignore "error"
    ug.set_cells(voxel_type, voxels)


    return ug

'''
1) Accepts a molecule type
2) Sifts through all output sets to find where it exists
3) Returns the concentrations of the molecule type in all sets it was found
'''
def get_voxel_molecule_conc(simData, moleculeType, out_location):
    # normalize moleculeType to str for consistent lookup
    if isinstance(moleculeType, (bytes, bytearray)):
        moleculeType = moleculeType.decode("utf-8")
    else:
        moleculeType = str(moleculeType)

    grid_points = len(getMorphologyGrid())
    samples = int(out_location[moleculeType]['samples'])
    outputSet = np.zeros((samples, grid_points))

    for currentSet, meta in out_location[moleculeType]['location'].items():
        molnum = meta['mol_index']
        voxels = list(meta['elements'])

        # population array shape is (time, voxels_in_outset, mol_index) or similar
        tempSnapshot = np.array(simData['trial0']['output'][currentSet]['population'][:, :, molnum])

        # If time-length mismatches expected samples, handle defensively:
        t_len = tempSnapshot.shape[0]
        if t_len == samples:
            outputSet[:, voxels] = tempSnapshot
        elif t_len > samples:
            # temp has MORE time points than our target: truncate (keep earliest frames)
            print(f"Warning: truncating {currentSet} population from {t_len} -> {samples} time-steps")
            outputSet[:, voxels] = tempSnapshot[:samples, :]
        else:
            # temp has FEWER time points: pad final frame forward to match samples
            print(f"Warning: padding {currentSet} population from {t_len} -> {samples} time-steps")
            padded = np.zeros((samples, tempSnapshot.shape[1]))
            padded[:t_len, :] = tempSnapshot
            # extend last frame forward
            if t_len > 0:
                padded[t_len:, :] = tempSnapshot[-1, :][None, :].repeat(samples - t_len, axis=0)
            outputSet[:, voxels] = padded

    # convert population -> concentration using robust voxel volumes
    outputSetConcs = population_to_concentration(outputSet, get_voxel_volumes())
    return outputSetConcs


'''
helper function to get the voxel volumes; looks for common places in the
HDF5 structure where this data might be stored, and falls back to uniform 1.0 if not found
'''
def get_voxel_volumes():
    # try common places to locate voxel volumes in the HDF5
    grid = simData['model']['grid']
    # If grid is structured and has a 'volume' field:
    try:
        arr = np.array(grid)
        if arr.dtype.names and 'volume' in arr.dtype.names:
            return arr['volume'].astype(float)
    except Exception:
        pass
    # fallback: look for model/grid/volume dataset
    try:
        return np.array(simData['model']['grid']['volume'])
    except Exception:
        # last-resort: use uniform 1.0 for all voxels (not ideal)
        n = len(getMorphologyGrid())
        return np.ones(n)


'''
Converts raw number of molecules(population) to an actual molecular concentration density
'''
def population_to_concentration(pop_list, voxel_volumes):
    conc_list = np.zeros(np.shape(pop_list))

    #Iterate through one timeframe of pop_list to divide each voxel's population by the ~[grid][voxel volume]
    for z, pop_snapshot in enumerate(pop_list): #
        for i, (a,b) in enumerate(zip(pop_snapshot, voxel_volumes)):
            conc_list[z][i] = (a/b) / mol_per_nM_u3

    #print out a volumes[3] volumes[6] & conc_list[3] [6] with two types
    return conc_list

'''
Class for setting options on the scale bar to the left of the model view
'''
class colorBarInputDialog(QWidget):

    def __init__(self):
        super(colorBarInputDialog, self).__init__()
        self.setGeometry(0,50,200,50)
        layout = QFormLayout()

        self.minLabel = QLabel(str(mayavi_widget_list[window.viewIndex-1].colorbar_min))
        self.btnMin = QPushButton("Min:")
        self.btnMin.clicked.connect(self.getMin)
        layout.addRow(self.btnMin,self.minLabel)

        self.maxLabel = QLabel(str(mayavi_widget_list[window.viewIndex-1].colorbar_max))
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
            viewer = mayavi_widget_list[window.viewIndex-1]
            viewer.surf.module_manager.scalar_lut_manager.data_range = [0, np.max(viewer.population)]
            if viewer.colorBar is not None:
                viewer.colorBar.visible = True
            self.minLabel.setText("0")
            self.maxLabel.setText(str(np.max(viewer.population)))
        except Exception:
            self.msg = QMessageBox()
            self.msg.setIcon(QMessageBox.Information)
            self.msg.setText("Select a Molecule First")
            self.msg.setWindowTitle("No Existing Defaults")

    def applyChanges(self):
        try:
            newMin, newMax = float(self.minLabel.text()), float(self.maxLabel.text())
            viewer = mayavi_widget_list[window.viewIndex - 1]
            viewer.surf.module_manager.scalar_lut_manager.data_range = [newMin, newMax]
            lut = viewer.surf.module_manager.scalar_lut_manager.lut
            # If a colorbar object exists, try to keep its LUT in sync (mayavi colorbar usually uses same LUT)
            if getattr(viewer, "colorBar", None) is not None:
                try:
                    viewer.colorBar.module_manager.scalar_lut_manager.lut = lut
                except Exception:
                    pass
            if self.leScale.text().startswith("Logarithmic"):
                lut.scale = 'log10'
            else:
                lut.scale = 'linear'
        except Exception:
            self.msg = QMessageBox()
            self.msg.setIcon(QMessageBox.Information)
            self.msg.setText("Select a Molecule First")
            self.msg.setWindowTitle("No Existing Colorbar")
    def closePopup(self):
        self.close()
'''
Help window to explain how to optimally use visualizer; needs work.
'''
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

'''
This function runs the animation portion of the visualizer

The "@mlab.animate" code above it indicates that anim
is a decorator function of the original mayavi function named animate
Decorators essentially work as wrappers, modifying the behavior of the code
before and after the target function, augmenting the original functionality.
In short, "anim" does what "animate" does, but with its own specifications.

(delay=x) sets the speed where x is # of miliseconds between each frame.
'''
@mlab.animate(delay=10)
def anim(simData, moleculeType, viewer_index):
    # --- Normalize molecule name so GUI strings match HDF5 keys ---
    if isinstance(moleculeType, (bytes, bytearray)):
        moleculeType = moleculeType.decode("utf-8")
    else:
        moleculeType = str(moleculeType)

    if viewer_index < 0 or viewer_index >= len(mayavi_widget_list):
        return
    viewer = mayavi_widget_list[viewer_index]
    # --- Gather molecule metadata and output locations for the entire model ---

    molecule_list = getMoleculeList(simData)
    grid_pts_len = len(getMorphologyGrid())
    out_location, dt_arr, samples_arr = get_mol_info(simData, molecule_list, grid_pts_len)
    # --- Ensure selected molecule exists (fallback to case-insensitive match) ---

    if moleculeType not in out_location:
        for key in out_location:
            if key.lower() == moleculeType.lower():
                moleculeType = key
                break
    if moleculeType not in out_location:
        return
    # --- Load voxel-wise concentration data and animation length ---

    viewer.population = get_voxel_molecule_conc(simData, moleculeType, out_location)
        # sanity check: population columns vs voxel volumes
    if viewer.population.ndim == 2:
        vv = get_voxel_volumes()
        if viewer.population.shape[1] != len(vv):
            print("Warning: population voxel count != voxel_volumes count", viewer.population.shape, len(vv))

    viewer.iterations = int(out_location[moleculeType]['samples'])
    dt = float(out_location[moleculeType]['dt'])
    # --- Build the 3D morphology mesh (unstructured grid) ---

    viewer.ug = create_morphology(simData)
    # --- Extract voxel volumes for population->concentration conversion (robust fallback) ---

    try:
        grid_obj = simData['model']['grid']
        try:
            arr = np.array(grid_obj)
            if getattr(arr.dtype, "names", None) and 'volume' in arr.dtype.names:
                voxel_volumes = arr['volume'].astype(float)
            else:
                voxel_volumes = np.array(simData['model']['grid']['volume'])
        except Exception:
            voxel_volumes = np.array(simData['model']['grid']['volume'])
    except Exception:
        voxel_volumes = np.ones(grid_pts_len)

    if viewer.population is None or viewer.population.size == 0:
        return
    # --- Compute global min/max across the ENTIRE animation for the colorbar ---

    global_min = float(np.min(viewer.population))
    global_max = float(np.max(viewer.population))
    viewer.colorbar_min, viewer.colorbar_max = global_min, global_max

    # --- Initialize visualization with first animation frame ---

    first_frame = viewer.population[0, :]
    viewer.ug.point_data.scalars = np.repeat(first_frame, 8)
    viewer.ug.point_data.scalars.name = 'concentrations'
    viewer.ug.modified()
    # --- Create surface and bind color mapping to global data range ---

    viewer.surf = mlab.pipeline.surface(viewer.ug, opacity=1, colormap='hot')

    # Lock the LUT to the global min/max so the colorbar is stable across all frames
    viewer.surf.module_manager.scalar_lut_manager.use_default_range = False
    viewer.surf.module_manager.scalar_lut_manager.data_range = [global_min, global_max]

    # Attach the initial scalars explicitly to the surface's mlab_source so Mayavi updates the colors
    try:
        init_scalars = np.repeat(first_frame, 8)
        viewer.surf.mlab_source.set(scalars=init_scalars)
        viewer.surf.mlab_source.dataset.point_data.scalars.name = 'concentrations'
        viewer.surf.mlab_source.update()
    except Exception:
        viewer.ug.point_data.scalars = np.repeat(first_frame, 8)
        viewer.ug.point_data.scalars.name = 'concentrations'
        viewer.ug.modified()

    # Re-lock the range after setting scalars (mlab_source.set can auto-adjust it)
    viewer.surf.module_manager.scalar_lut_manager.data_range = [global_min, global_max]

    # --- Attach a colorbar directly to the surface ---

    try:
        if getattr(viewer, "colorBar", None) is not None:
            try:
                viewer.colorBar.visible = False
            except Exception:
                pass
        viewer.colorBar = mlab.colorbar(object=viewer.surf, title='Concentration', orientation='vertical')
        viewer.colorBar.visible = True
    except Exception:
        viewer.colorBar = None

    if viewer.getCurrentFrame() is None:
        viewer.setCurrentFrame(0)
    # --- Main animation loop: update scalars, advance frame, update UI ---

    while viewer.getCurrentFrame() < viewer.iterations:
        if viewer_index >= len(mayavi_widget_list):
            break
        v = mayavi_widget_list[viewer_index]

        frame_idx = v.getCurrentFrame()
        frame_idx = max(0, min(frame_idx, v.iterations - 1))

        concentrations = v.population[frame_idx, :]
        scalars_pts = np.repeat(concentrations, 8)

        # Preferred: update the surface's mlab_source so the display updates
        try:
            v.surf.mlab_source.set(scalars=scalars_pts)
            try:
                v.surf.mlab_source.dataset.point_data.scalars.name = 'concentrations'
                v.surf.mlab_source.update()
            except Exception:
                pass
        except Exception:
            v.ug.point_data.scalars = scalars_pts
            v.ug.point_data.scalars.name = 'concentrations'
            v.ug.modified()

        # Re-lock the colorbar to global range after every frame update
        try:
            v.surf.module_manager.scalar_lut_manager.data_range = [v.colorbar_min, v.colorbar_max]
        except Exception:
            pass

        v.setCurrentFrame(frame_idx + 1)

        # Update the time label for whichever viewer is currently selected
        if viewer_index == window.viewIndex - 1:
            window.progress_label.setText(f"{v.getCurrentFrame()/1000:.3f}s")

            try:
                pct = int((v.getCurrentFrame() / max(1, v.iterations)) * 100)
                progress_bar.setValue(pct)

                try:
                    window.progress_slider.blockSignals(True)
                    window.progress_slider.setValue(pct)
                    window.progress_slider.blockSignals(False)
                except Exception:
                    pass
            except Exception:
                pass

        yield
    # --- Reset animation to beginning after completion ---

    if viewer.getCurrentFrame() >= (viewer.iterations - 1):
        viewer.setCurrentFrame(0)


'''A view embedded in the window to contain an instance of the model'''
class MayaviQWidget(QtGui.QWidget):
    def __init__(self, parent):

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
        self.iterations = 0

    def home(self):
        self.animator = None
        self.currentFrame = 0
        self.ug = None
        self.colorbar_min, self.colorbar_max = 0, 0
        self.colorBar = None
        self.surf = None
        self.population = None
        self.iterations = 0


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
'''
Main window of QtGui;
Overall layout of widgets are organized here and
most operations will address objects created here.
'''
class Window(QtGui.QMainWindow):

    #Main Menus should go here - things which appear at startup.
    def __init__(self):
        super(Window, self).__init__()
        self.setGeometry(50, 50, 1100, 800)
        self.setWindowTitle("NeuoRD Visualizer" + " - " + fileName)

        #Main Menu details for "Add a Viewer" button
        addAction = QtGui.QAction("&Add a Viewer -", self)
        addAction.setShortcut("Ctrl+A")
        addAction.setStatusTip('Add Items to Visualizer')
        addAction.triggered.connect(self.add_view) #.triggered = .clicked

        #Main Menu details for "Exit" button
        exitAction = QtGui.QAction("&Exit -", self)
        exitAction.setShortcut("Ctrl+Q")
        exitAction.setStatusTip('Close Application')
        exitAction.triggered.connect(self.close_application) #.triggered = .clicked

        #Main Menu details for "Min/Max Range for Colorbar" button
        minMaxColorBarAction = QtGui.QAction("&Min/Max Range for Colorbar -", self)
        minMaxColorBarAction.setShortcut("Ctrl+M")
        minMaxColorBarAction.setStatusTip('Change the displayed minimum & maximum ranges on the color bar.')
        minMaxColorBarAction.triggered.connect(self.changeMinMaxColorBar)

        #Main Menu details for "Select a Model" button
        ## note: this should become obsolete eventually once viewer specific molecule section & start/stop buttons are added
        selectModelAction = QtGui.QAction("&Select a Model -", self)
        selectModelAction.setShortcut("Ctrl+S")
        selectModelAction.setStatusTip('Select another view to simulate.')
        selectModelAction.triggered.connect(self.select_view)

        #Main Menu details for "Help" button
        helpAction = QtGui.QAction("&Help -", self)
        helpAction.setShortcut("Ctrl+H")
        helpAction.setStatusTip('Learn More About How to Use the Visualizer')
        helpAction.triggered.connect(self.help_action) #.triggered = .clicked

        mainMenu = self.menuBar()

        #Adds Main Menu Toolbar "File" & assigns items, created above, to its dropdown.
        fileMenu = mainMenu.addMenu('&File')
        fileMenu.addAction(addAction)
        fileMenu.addAction(exitAction)

        #Adds Main Menu Toolbar "Edit" & assigns items, created above, to its dropdown.
        editMenu = mainMenu.addMenu('&Edit')
        editMenu.addAction(minMaxColorBarAction)
        editMenu.addAction(selectModelAction)

        #Adds Main Menu Toolbar "Help" & assigns items, created above, to its dropdown.
        helpMenu = mainMenu.addMenu('&Help')
        helpMenu.addAction(helpAction)

        #Create Progress bar
        self.progress_slider_label = progress_slider_label
        self.progress_label = progress_label

        #Index of current Viewer
        self.viewIndex = 0
        #Indexes for current row and column positions.
        self.rowIndex = 0
        self.columnIndex = 0
        #Tally of total Viewers added
        self.viewTally = 1

        self.statusBar()
        self.home()

    #Similar to init; home loads objects at startup.
    def home(self):

        toolBarColorBarMinMax = QtGui.QAction(QtGui.QIcon('colorBarIcon.png'), "Set Min/Max of ColorBar", self)
        toolBarColorBarMinMax.setStatusTip('Change the default min/max range on the color bar.')
        toolBarColorBarMinMax.triggered.connect(self.changeMinMaxColorBar)
        toolBarAddView = QtGui.QAction(QtGui.QIcon('addModelIcon.png'), "Add a Viewer", self)
        toolBarAddView.setStatusTip('Add another window to the visualizer.')
        toolBarAddView.triggered.connect(self.add_view)
        toolBarSelectView = QtGui.QAction(QtGui.QIcon('selectViewIcon.png'), "Select a View", self)
        toolBarSelectView.setStatusTip('Select an Existing Model to Visualize.')
        toolBarSelectView.triggered.connect(self.select_view)
        toolBarHelp = QtGui.QAction(QtGui.QIcon('helpIcon.png'), "Help", self)
        toolBarHelp.setStatusTip("Learn More About How to Use the Visualizer")
        toolBarHelp.triggered.connect(self.help_action)

        self.toolBar = self.addToolBar("ToolBar")
        self.toolBar.addAction(toolBarColorBarMinMax)
        self.toolBar.addAction(toolBarAddView)
        self.toolBar.addAction(toolBarSelectView)
        self.toolBar.addAction(toolBarHelp)

        self.show()

    def close_application(self):
        choice = QMessageBox.question(
            self,
            "Exit",
            "Are you sure?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
        )
        if choice == QMessageBox.StandardButton.Yes:
            sys.exit()


    def changeMinMaxColorBar(self):
        if not mayavi_widget_list or self.viewIndex < 1 or self.viewIndex > len(mayavi_widget_list):
            QMessageBox.information(self, "No viewer", "Add/select a viewer first.")
            return
        self.newEditWindow = colorBarInputDialog()
        self.newEditWindow.show()

    def help_action(self):
        self.newHelpWindow = helpWindow()
        self.newHelpWindow.show()

    #Adds new molecule visualization view
    def add_view(self):
        viewer_index = self.viewTally - 1  # 0-based index for the new viewer
        mayavi_widget_list.append(MayaviQWidget(container))
        if self.viewTally % 2 != 0:
            self.columnIndex=0
            layout.addWidget(populate_comboBox(viewer_index), self.rowIndex, self.columnIndex)
            layout.addWidget(mayavi_widget_list[viewer_index], self.rowIndex+1, self.columnIndex)
        else:
            self.columnIndex=1
            layout.addWidget(populate_comboBox(viewer_index), self.rowIndex, self.columnIndex)
            layout.addWidget(mayavi_widget_list[viewer_index], self.rowIndex+1, self.columnIndex)
            self.rowIndex += 2


        self.viewTally += 1
        self.viewIndex += 1

    def select_view(self):
        text, ok = QInputDialog.getText(self, 'Molecule Selection', 'Enter Window # to Simulate in  (1-' +str(window.viewTally) +")" )
        if ok:
            self.viewIndex = int(text)
    #This is where the animation portion of the program is called.

    def molecule_selected_for_viewer(self, text, viewer_index):
        """Start animation for a specific viewer. Each viewer runs independently."""
        if not text:
            return
        if viewer_index < 0 or viewer_index >= len(mayavi_widget_list):
            return

        viewer = mayavi_widget_list[viewer_index]

        # close previous animator for THIS viewer (if any)
        if getattr(viewer, "animator", None) is not None:
            try:
                viewer.animator.close()
            except Exception:
                pass

        # start a fresh animator bound to this viewer
        viewer.animator = anim(simData, text, viewer_index)

        # also set this viewer as the selected one for slider/progress bar
        self.viewIndex = viewer_index + 1

        print(f"Started animation for '{text}' in viewer {viewer_index + 1}")

    def slider_movement(self):
        if not mayavi_widget_list:
            return

        # ensure viewIndex maps to a valid list index
        idx = max(0, min(self.viewIndex - 1, len(mayavi_widget_list) - 1))
        viewer = mayavi_widget_list[idx]

        try:
            iterations = int(getattr(viewer, "iterations", 0))
        except Exception:
            return
        if iterations <= 0:
            return

        position = self.progress_slider.value()
        x = int((position / 100.0) * iterations)
        x = max(0, min(x, iterations - 1))

        viewer.setCurrentFrame(x)

        try:
            self.progress_slider_label.setText(f"{x/1000:.3f}s")
        except Exception:
            self.progress_slider_label.setText(str(x/1000) + "s")



    def resetAnimation(self, resetButtonNumber):
        # clamp index for safety
        idx = max(0, min(self.viewIndex - 1, len(mayavi_widget_list) - 1))
        viewer = mayavi_widget_list[idx]

        # reset viewer internal frame counter
        try:
            viewer.setCurrentFrame(0)
        except Exception:
            viewer.currentFrame = 0

        # update the global UI controls (single slider + progressbar)
        try:
            window.progress_slider.blockSignals(True)
            window.progress_slider.setValue(0)
            window.progress_slider.blockSignals(False)
        except Exception:
            pass

        try:
            progress_bar.setValue(0)
        except Exception:
            pass







'''
Returns the list of molecule types available in the h5 simulation file
'''
def getMoleculeList(simData):
    raw = simData['model']['species'][:]   # <-- NOTE the [:] to get the array
    return [m.decode("utf-8") if isinstance(m, (bytes, bytearray)) else str(m) for m in raw]

'''
Returns the h5 simulation file
'''
def get_h5simData(fileName):
    simData = h5.File(fileName,"r")
    return simData
'''
Returns the grid of the h5 simulation file itself
'''
def getMorphologyGrid():
    grid = simData['model']['grid']
    # if grid is a dataset containing a structured array, convert to numpy recarray
    try:
        arr = np.array(grid)
        # If arr is structured array, return it as recarray
        if arr.dtype.names:
            return arr.view(np.recarray)
        return arr
    except Exception:
        # fallback: try to read elements list
        try:
            return grid[:]
        except Exception:
            raise

'''
Returns the container portion of the QtWindow so it can be accessed during runtime.
'''
def getQtWindow():
    return container

'''
Searches the list of molecules and returns the corresponding index. Returns -1 if not found.
'''
def get_mol_index(simData, outputSet, molecule):
    species = simData['model']['output'][outputSet]['species'][:]
    # If species are bytes, compare bytes
    if len(species) > 0 and isinstance(species[0], (bytes, bytearray)):
        if isinstance(molecule, str):
            molecule = molecule.encode("utf-8")
    else:
        if isinstance(molecule, (bytes, bytearray)):
            molecule = molecule.decode("utf-8")

    indices = np.where(species == molecule)[0]
    return int(indices[0]) if len(indices) == 1 else -1

'''
This function returns various information for a molecule type: its concentrations(samples), time intervals(dt), locations(out_location)
'''
def get_mol_info(simData, plot_molecules, grid_points):
    # Normalize plot_molecules to python str
    plot_molecules = [
        m.decode("utf-8") if isinstance(m, (bytes, bytearray)) else str(m)
        for m in plot_molecules
    ]

    # If grid_points is a dataset/array, convert to length (int)
    try:
        if not isinstance(grid_points, int):
            grid_points = len(grid_points)
    except Exception:
        grid_points = int(grid_points)

    outputsets = list(simData['model']['output'].keys())
    dt = np.zeros((len(plot_molecules)))
    samples = np.zeros((len(plot_molecules)), dtype=int)
    out_location = {}

    # iterate through the provided plot_molecules (now str)
    for imol, molecule in enumerate(plot_molecules):
        temp_dict = {}
        tot_voxels = 0
        samples_seen = []
        dt_seen = []

        # iterate outputsets (search last->first to mirror older behavior)
        for outset in outputsets[::-1]:
            mol_index = get_mol_index(simData, outset, molecule)
            if mol_index > -1:
                # record samples and dt for this outset
                try:
                    n_samples = len(simData['trial0']['output'][outset]['times'])
                    samples_seen.append(int(n_samples))
                except Exception:
                    samples_seen.append(0)
                try:
                    # convert msec->sec like original code
                    dt_seen.append(float(simData['trial0']['output'][outset]['times'][1]) / 1000.0)
                except Exception:
                    dt_seen.append(0.0)

                tot_voxels += len(simData['model']['output'][outset]['elements'])
                temp_dict[outset] = {
                    'mol_index': mol_index,
                    'elements': simData['model']['output'][outset]['elements'][:].tolist()
                }

        if len(temp_dict) > 0:
            # choose the maximum samples seen (robust when outputs have different lengths)
            samples[imol] = int(max(samples_seen) if samples_seen else 0)
            # choose a dt: prefer the dt corresponding to the largest sample count, else first nonzero
            dt_val = 0.0
            if samples_seen:
                # pick the index of max samples if possible
                try:
                    idx_max = int(np.argmax(samples_seen))
                    if idx_max < len(dt_seen):
                        dt_val = float(dt_seen[idx_max])
                except Exception:
                    pass
            if dt_val == 0.0:
                # fallback to first non-zero dt
                for d in dt_seen:
                    if d:
                        dt_val = float(d)
                        break

            out_location[molecule] = {
                'samples': int(samples[imol]),
                'dt': float(dt_val),
                'voxels': int(tot_voxels),
                'location': temp_dict
            }
        else:
            # fallback: search the first outset (old behavior)
            outset = outputsets[0]
            print("************* MOLECULE", molecule, " NOT IN REGULAR OUTPUT SETS !")
            mol_index = get_mol_index(simData, outset, molecule)
            samples[imol] = len(simData['trial0']['output'][outset]['times'])
            dt[imol] = simData['trial0']['output'][outset]['times'][1] / 1000.0
            temp_dict[outset] = {
                'mol_index': mol_index,
                'elements': simData['model']['output'][outset]['elements'][:] .tolist()
            }
            out_location[molecule] = {
                'samples': int(samples[imol]),
                'dt': float(dt[imol]),
                'voxels': int(grid_points),
                'location': temp_dict
            }

    return out_location, dt, samples



'''
Fills the dropdown window with choices for molecules
'''
def populate_comboBox(viewer_index):
    comboBoxItemModel = QStandardItemModel() #Required for searchable comboBox
    comboBox = ExtendedCombo()
    moleculeList = getMoleculeList(simData)  # HDF5 bytes
    moleculeList = [
        (m.decode("utf-8") if isinstance(m, (bytes, bytearray)) else str(m))
        for m in moleculeList
    ]
    moleculeList = sorted(moleculeList)

    for i, moleculeType in enumerate(moleculeList):
        item = QStandardItem(moleculeType)
        comboBoxItemModel.setItem(i, 0, item)

    comboBox.setModel(comboBoxItemModel)
    comboBox.setModelColumn(0)

    # Bind this combobox to its specific viewer using a lambda with captured index
    idx = viewer_index  # capture for closure
    comboBox.textActivated.connect(lambda text, vi=idx: window.molecule_selected_for_viewer(text, vi))

    return comboBox


if __name__ == "__main__":

    #Passes filename argument to NeuroRDViz.py to find the desired model to visualize. Must be in the same folder.
    #For example, the following could be entered into the command prompt:
    # "python NeuroRDViz.py Model_CamKIInew_pDglUchi5s-dhpg5.h5
    try:
        fileName=fname # @UndefinedVariable
    except NameError:
        fileName = sys.argv[1]

    simData = get_h5simData(fileName)

    #Creating instances of mayavi UI

    app = QApplication.instance() or QApplication(sys.argv)
    container = QtGui.QWidget()
    layout = QtGui.QGridLayout(container)

    progress_label = QtGui.QLabel(container)
    progress_slider_label = QtGui.QLabel(container)
    window = Window()

    # create storage before adding views
    mayavi_widget_list = []

    # Put the first view at the top
    window.rowIndex = 0
    window.columnIndex = 0
    window.viewIndex = 0
    window.viewTally = 1

    # Add the first viewer (combo + mayavi widget)
    window.add_view()

    progress_bar = QtGui.QProgressBar()
    window.progress_slider = QSlider(Qt.Horizontal)
    window.progress_slider.valueChanged.connect(window.slider_movement)



    #These lines place the respective widgets into the overall layout that allows you to place items in appropriate positions
    #e.g. comboxBox will be added to the 1st row and 1st column with the line:
    #layout.addWidget(comboBox, 0, 0)
    mol_type_label_list = []
    mol_type_label_list.append(QtGui.QLabel())
    #layout.addWidget(comboBox, 0, 0)  # 0,0 = top left widget location, 0,1 = one to the right of it, etc.
    #layout.addWidget(mol_type_label_list[0], 0,1)
    #layout.addWidget(mayavi_widget_list[window.viewIndex-1], 4, 1) # Visualization of morphology
    #layout.addWidget(reset_button_list[0], 5, 1)
    layout.addWidget(progress_label, 2,0)
    layout.addWidget(progress_bar, 2, 0)
    layout.addWidget(progress_slider_label,3, 0)
    layout.addWidget(window.progress_slider, 3, 0)
    #mayavi_widget_list.append(MayaviQWidget(container))

    container.show()

    window.setCentralWidget(container)
    window.show()
    app.exec() # Start the main event loop.
