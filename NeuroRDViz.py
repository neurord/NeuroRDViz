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
        self.ug = create_morphology(simData)
        s = self.scene.mlab  # scene-specific mlab — correct because scene IS active right now
        self.surf = s.pipeline.surface(self.ug, opacity=1)
        s.pipeline.surface(s.pipeline.extract_edges(self.surf), color=(0, 0, 0))
        s.axes(self.surf, nb_labels=7)

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
This function sets up the animation for a specific viewer.
It loads molecule data, configures the EXISTING surface (created in
Visualization.update_plot, guaranteed to be in the correct scene),
and prepares the viewer for frame-by-frame animation via QTimer.
Returns True if setup succeeded, False otherwise.
'''
def anim_setup(simData, moleculeType, viewer_index):
    # --- Normalize molecule name so GUI strings match HDF5 keys ---
    if isinstance(moleculeType, (bytes, bytearray)):
        moleculeType = moleculeType.decode("utf-8")
    else:
        moleculeType = str(moleculeType)

    if viewer_index < 0 or viewer_index >= len(mayavi_widget_list):
        return False
    viewer = mayavi_widget_list[viewer_index]

    # Ensure the scene has been activated and update_plot has run
    viz = viewer.visualization
    if not hasattr(viz, 'surf') or viz.surf is None:
        print(f"Viewer {viewer_index}: scene not yet activated, cannot start animation")
        return False

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
        return False
    # --- Load voxel-wise concentration data and animation length ---

    viewer.population = get_voxel_molecule_conc(simData, moleculeType, out_location)
        # sanity check: population columns vs voxel volumes
    if viewer.population.ndim == 2:
        vv = get_voxel_volumes()
        if viewer.population.shape[1] != len(vv):
            print("Warning: population voxel count != voxel_volumes count", viewer.population.shape, len(vv))

    viewer.iterations = int(out_location[moleculeType]['samples'])

    if viewer.population is None or viewer.population.size == 0:
        return False
    # --- Compute global min/max across the ENTIRE animation for the colorbar ---

    global_min = float(np.min(viewer.population))
    global_max = float(np.max(viewer.population))
    viewer.colorbar_min, viewer.colorbar_max = global_min, global_max

    # --- Use the EXISTING surface from update_plot (already in the correct scene) ---
    # This is the key fix: we never create new pipeline objects here, so we
    # bypass mlab's global-figure routing entirely.

    surf = viz.surf
    ug = viz.ug
    viewer.surf = surf  # alias for anim_step convenience
    viewer.ug = ug

    # --- Set initial concentration scalars on the existing UG ---

    first_frame = viewer.population[0, :]
    init_scalars = np.repeat(first_frame, 8)
    ug.point_data.scalars = init_scalars
    ug.point_data.scalars.name = 'concentrations'
    ug.modified()

    # --- Configure the existing surface for concentration display ---

    lut_mgr = surf.module_manager.scalar_lut_manager
    lut_mgr.lut_mode = 'hot'
    lut_mgr.use_default_range = False
    lut_mgr.data_range = np.array([global_min, global_max])

    # Show the scalar bar (colorbar) on the existing surface's LUT manager
    try:
        lut_mgr.show_scalar_bar = True
        lut_mgr.scalar_bar.title = 'Concentration'
        viewer.colorBar = lut_mgr
    except Exception:
        viewer.colorBar = None

    # Force the scene to pick up the new scalars and render
    try:
        viz.scene.render()
    except Exception:
        pass

    viewer.setCurrentFrame(0)
    return True


def anim_step(viewer_index):
    """Advance one animation frame for the given viewer. Called by QTimer."""
    if viewer_index < 0 or viewer_index >= len(mayavi_widget_list):
        return
    v = mayavi_widget_list[viewer_index]

    if v.getCurrentFrame() >= v.iterations:
        # Animation complete — stop timer and reset
        if getattr(v, "anim_timer", None) is not None:
            v.anim_timer.stop()
        v.setCurrentFrame(0)
        return

    frame_idx = v.getCurrentFrame()
    frame_idx = max(0, min(frame_idx, v.iterations - 1))

    concentrations = v.population[frame_idx, :]
    scalars_pts = np.repeat(concentrations, 8)

    # Update the EXISTING UG's scalars directly — this UG is bound to a surface
    # that lives in this viewer's scene (created during update_plot), so the
    # VTK pipeline update will render in the CORRECT scene.
    scene = v.visualization.scene
    try:
        scene.disable_render = True
    except Exception:
        pass

    v.ug.point_data.scalars = scalars_pts
    v.ug.point_data.scalars.name = 'concentrations'
    v.ug.modified()

    # Re-lock the colorbar to global range
    try:
        v.surf.module_manager.scalar_lut_manager.data_range = np.array(
            [v.colorbar_min, v.colorbar_max])
    except Exception:
        pass

    try:
        scene.disable_render = False
    except Exception:
        pass

    # Force THIS viewer's scene to re-render
    try:
        scene.render()
    except Exception:
        pass

    v.setCurrentFrame(frame_idx + 1)

    # Update per-viewer progress bar, slider, and time label
    try:
        pct = int((v.getCurrentFrame() / max(1, v.iterations)) * 100)
        if v.progress_bar is not None:
            v.progress_bar.setValue(pct)
        if v.progress_label is not None:
            v.progress_label.setText(f"{v.getCurrentFrame()/1000:.3f}s")
        if v.progress_slider is not None:
            v.progress_slider.blockSignals(True)
            v.progress_slider.setValue(pct)
            v.progress_slider.blockSignals(False)
    except Exception:
        pass

    # Also update global progress for the currently selected viewer
    try:
        if viewer_index == window.viewIndex - 1:
            pct = int((v.getCurrentFrame() / max(1, v.iterations)) * 100)
            window.global_progress_bar.setValue(pct)
            window.global_progress_label.setText(f"{v.getCurrentFrame()/1000:.3f}s")
            window.global_progress_slider.blockSignals(True)
            window.global_progress_slider.setValue(pct)
            window.global_progress_slider.blockSignals(False)
    except Exception:
        pass


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
        self.anim_timer = None
        self.currentFrame = 0
        self.ug = None
        self.colorbar_min, self.colorbar_max = 0, 0
        self.colorBar = None
        self.surf = None
        self.population = None
        self.iterations = 0
        self.progress_bar = None
        self.progress_slider = None
        self.progress_label = None


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

        #Main Menu details for "Start All Animations" button
        startAllAction = QtGui.QAction("&Start All Animations -", self)
        startAllAction.setShortcut("Ctrl+G")
        startAllAction.setStatusTip('Start animations for all viewers simultaneously')
        startAllAction.triggered.connect(self.start_all_animations)

        #Main Menu details for "Stop All Animations" button
        stopAllAction = QtGui.QAction("S&top All Animations -", self)
        stopAllAction.setShortcut("Ctrl+T")
        stopAllAction.setStatusTip('Stop all running animations')
        stopAllAction.triggered.connect(self.stop_all_animations)

        #Main Menu details for "Help" button
        helpAction = QtGui.QAction("&Help -", self)
        helpAction.setShortcut("Ctrl+H")
        helpAction.setStatusTip('Learn More About How to Use the Visualizer')
        helpAction.triggered.connect(self.help_action) #.triggered = .clicked

        mainMenu = self.menuBar()

        #Adds Main Menu Toolbar "File" & assigns items, created above, to its dropdown.
        fileMenu = mainMenu.addMenu('&File')
        fileMenu.addAction(addAction)
        fileMenu.addAction(startAllAction)
        fileMenu.addAction(stopAllAction)
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

        toolBarStartAll = QtGui.QAction(QtGui.QIcon('startAllIcon.png'), "Start All Animations", self)
        toolBarStartAll.setStatusTip('Start animations for all viewers simultaneously (in unison)')
        toolBarStartAll.triggered.connect(self.start_all_animations)
        toolBarStopAll = QtGui.QAction(QtGui.QIcon('stopAllIcon.png'), "Stop All Animations", self)
        toolBarStopAll.setStatusTip('Stop all running animations')
        toolBarStopAll.triggered.connect(self.stop_all_animations)

        self.toolBar = self.addToolBar("ToolBar")
        self.toolBar.addAction(toolBarColorBarMinMax)
        self.toolBar.addAction(toolBarAddView)
        self.toolBar.addAction(toolBarSelectView)
        self.toolBar.addAction(toolBarStartAll)
        self.toolBar.addAction(toolBarStopAll)
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

        # Create combobox and store reference on the viewer
        comboBox = populate_comboBox(viewer_index)
        viewer = mayavi_widget_list[viewer_index]
        viewer.comboBox = comboBox

        # Create per-viewer Start and Stop buttons
        idx = viewer_index
        start_btn = QPushButton(f"Start Viewer {viewer_index + 1}")
        start_btn.clicked.connect(lambda checked, vi=idx: self.start_viewer_animation(vi))
        viewer.start_button = start_btn

        stop_btn = QPushButton(f"Stop Viewer {viewer_index + 1}")
        stop_btn.clicked.connect(lambda checked, vi=idx: self.stop_viewer_animation(vi))
        viewer.stop_button = stop_btn

        # Per-viewer progress bar, slider, and time label
        viewer.progress_label = QLabel("0.000s")
        viewer.progress_bar = QProgressBar()
        viewer.progress_bar.setValue(0)
        viewer.progress_slider = QSlider(Qt.Horizontal)
        viewer.progress_slider.setRange(0, 100)
        viewer.progress_slider.setValue(0)
        viewer.progress_slider.valueChanged.connect(
            lambda val, vi=idx: self.viewer_slider_movement(vi)
        )

        # Horizontal control bar: [combobox] [Start] [Stop]
        ctrl_widget = QWidget()
        ctrl_layout = QtWidgets.QHBoxLayout(ctrl_widget)
        ctrl_layout.setContentsMargins(0, 0, 0, 0)
        ctrl_layout.addWidget(comboBox, stretch=1)
        ctrl_layout.addWidget(start_btn)
        ctrl_layout.addWidget(stop_btn)

        # Per-viewer progress widget: [progress_bar | time_label] / [slider]
        progress_widget = QWidget()
        progress_vlayout = QtWidgets.QVBoxLayout(progress_widget)
        progress_vlayout.setContentsMargins(0, 0, 0, 0)
        progress_vlayout.setSpacing(2)
        prog_top = QWidget()
        prog_top_layout = QtWidgets.QHBoxLayout(prog_top)
        prog_top_layout.setContentsMargins(0, 0, 0, 0)
        prog_top_layout.addWidget(viewer.progress_bar, stretch=1)
        prog_top_layout.addWidget(viewer.progress_label)
        progress_vlayout.addWidget(prog_top)
        progress_vlayout.addWidget(viewer.progress_slider)

        if self.viewTally % 2 != 0:
            self.columnIndex=0
        else:
            self.columnIndex=1

        layout.addWidget(ctrl_widget, self.rowIndex, self.columnIndex)
        layout.addWidget(mayavi_widget_list[viewer_index], self.rowIndex+1, self.columnIndex)
        layout.addWidget(progress_widget, self.rowIndex+2, self.columnIndex)

        if self.viewTally % 2 == 0:
            self.rowIndex += 3

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

        # stop previous animation timer for THIS viewer (if any)
        if getattr(viewer, "anim_timer", None) is not None:
            try:
                viewer.anim_timer.stop()
            except Exception:
                pass
            viewer.anim_timer = None

        # set up the animation data and surface in the correct scene
        if not anim_setup(simData, text, viewer_index):
            return

        # create a QTimer that drives this viewer's animation independently
        timer = QtCore.QTimer()
        timer.setInterval(10)  # 10ms between frames, same as old @mlab.animate(delay=10)
        timer.timeout.connect(lambda vi=viewer_index: anim_step(vi))
        viewer.anim_timer = timer
        timer.start()

        # also set this viewer as the selected one for slider/progress bar
        self.viewIndex = viewer_index + 1

        print(f"Started animation for '{text}' in viewer {viewer_index + 1}")

    def start_viewer_animation(self, viewer_index):
        """Start animation for a specific viewer using its combobox selection."""
        if viewer_index < 0 or viewer_index >= len(mayavi_widget_list):
            return
        viewer = mayavi_widget_list[viewer_index]
        comboBox = getattr(viewer, 'comboBox', None)
        if comboBox is None:
            return
        molecule = comboBox.currentText()
        if not molecule:
            return
        self.molecule_selected_for_viewer(molecule, viewer_index)

    def start_all_animations(self):
        """Start animations for all viewers simultaneously (in unison)."""
        # Phase 1: stop all existing timers and reset frames
        for vi in range(len(mayavi_widget_list)):
            viewer = mayavi_widget_list[vi]
            if getattr(viewer, "anim_timer", None) is not None:
                try:
                    viewer.anim_timer.stop()
                except Exception:
                    pass
                viewer.anim_timer = None
            viewer.setCurrentFrame(0)

        # Phase 2: set up all animations, then start all timers together
        timers_to_start = []
        for vi in range(len(mayavi_widget_list)):
            viewer = mayavi_widget_list[vi]
            comboBox = getattr(viewer, 'comboBox', None)
            if comboBox is None:
                continue
            molecule = comboBox.currentText()
            if not molecule:
                continue
            if not anim_setup(simData, molecule, vi):
                continue
            timer = QtCore.QTimer()
            timer.setInterval(10)
            timer.timeout.connect(lambda v=vi: anim_step(v))
            viewer.anim_timer = timer
            timers_to_start.append(timer)

        for timer in timers_to_start:
            timer.start()

        if mayavi_widget_list:
            self.viewIndex = 1
        print("Started all viewer animations simultaneously")

    def stop_viewer_animation(self, viewer_index):
        """Stop animation for a specific viewer."""
        if viewer_index < 0 or viewer_index >= len(mayavi_widget_list):
            return
        viewer = mayavi_widget_list[viewer_index]
        if getattr(viewer, "anim_timer", None) is not None:
            try:
                viewer.anim_timer.stop()
            except Exception:
                pass
            viewer.anim_timer = None
        print(f"Stopped animation in viewer {viewer_index + 1}")

    def stop_all_animations(self):
        """Stop animations for all viewers."""
        for vi in range(len(mayavi_widget_list)):
            viewer = mayavi_widget_list[vi]
            if getattr(viewer, "anim_timer", None) is not None:
                try:
                    viewer.anim_timer.stop()
                except Exception:
                    pass
                viewer.anim_timer = None
        print("Stopped all viewer animations")

    def global_slider_movement(self):
        """When the global slider is dragged, move ALL viewers to that percentage."""
        if not mayavi_widget_list:
            return
        position = self.global_progress_slider.value()
        for vi in range(len(mayavi_widget_list)):
            viewer = mayavi_widget_list[vi]
            try:
                iterations = int(getattr(viewer, "iterations", 0))
            except Exception:
                continue
            if iterations <= 0:
                continue
            x = int((position / 100.0) * iterations)
            x = max(0, min(x, iterations - 1))
            viewer.setCurrentFrame(x)
            # Sync per-viewer progress widgets
            try:
                pct = int((x / max(1, iterations)) * 100)
                viewer.progress_bar.setValue(pct)
                viewer.progress_label.setText(f"{x/1000:.3f}s")
                viewer.progress_slider.blockSignals(True)
                viewer.progress_slider.setValue(pct)
                viewer.progress_slider.blockSignals(False)
            except Exception:
                pass

    def viewer_slider_movement(self, viewer_index):
        """When a per-viewer slider is dragged, move only that viewer."""
        if viewer_index < 0 or viewer_index >= len(mayavi_widget_list):
            return
        viewer = mayavi_widget_list[viewer_index]
        try:
            iterations = int(getattr(viewer, "iterations", 0))
        except Exception:
            return
        if iterations <= 0:
            return
        position = viewer.progress_slider.value()
        x = int((position / 100.0) * iterations)
        x = max(0, min(x, iterations - 1))
        viewer.setCurrentFrame(x)
        try:
            viewer.progress_label.setText(f"{x/1000:.3f}s")
            viewer.progress_bar.setValue(int((x / max(1, iterations)) * 100))
        except Exception:
            pass



    def resetAnimation(self, resetButtonNumber):
        # clamp index for safety
        idx = max(0, min(self.viewIndex - 1, len(mayavi_widget_list) - 1))
        viewer = mayavi_widget_list[idx]

        # reset viewer internal frame counter
        try:
            viewer.setCurrentFrame(0)
        except Exception:
            viewer.currentFrame = 0

        # reset per-viewer progress widgets
        try:
            if viewer.progress_bar is not None:
                viewer.progress_bar.setValue(0)
            if viewer.progress_label is not None:
                viewer.progress_label.setText("0.000s")
            if viewer.progress_slider is not None:
                viewer.progress_slider.blockSignals(True)
                viewer.progress_slider.setValue(0)
                viewer.progress_slider.blockSignals(False)
        except Exception:
            pass

        # update global progress controls
        try:
            window.global_progress_slider.blockSignals(True)
            window.global_progress_slider.setValue(0)
            window.global_progress_slider.blockSignals(False)
            window.global_progress_bar.setValue(0)
            window.global_progress_label.setText("0.000s")
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

    # Main container: VBoxLayout with viewer grid on top, global controls on bottom
    container = QtGui.QWidget()
    main_layout = QtWidgets.QVBoxLayout(container)
    main_layout.setContentsMargins(4, 4, 4, 4)
    main_layout.setSpacing(6)

    # Viewer grid (comboboxes, mayavi widgets, per-viewer progress)
    viewer_grid_widget = QtGui.QWidget()
    layout = QtGui.QGridLayout(viewer_grid_widget)
    main_layout.addWidget(viewer_grid_widget, stretch=1)

    # Dummy labels required by Window.__init__ (replaced by per-viewer labels)
    progress_label = QtGui.QLabel()
    progress_slider_label = QtGui.QLabel()

    window = Window()

    # create storage before adding views
    mayavi_widget_list = []

    # Put the first view at the top
    window.rowIndex = 0
    window.columnIndex = 0
    window.viewIndex = 0
    window.viewTally = 1

    # Add the first viewer (combo + mayavi widget + per-viewer progress)
    window.add_view()

    # ── Global controls at bottom ────────────────────────────────────────
    global_controls = QtGui.QWidget()
    gc_layout = QtWidgets.QVBoxLayout(global_controls)
    gc_layout.setContentsMargins(0, 0, 0, 0)
    gc_layout.setSpacing(4)

    # Start/Stop All buttons
    btn_row = QWidget()
    btn_row_layout = QtWidgets.QHBoxLayout(btn_row)
    btn_row_layout.setContentsMargins(0, 0, 0, 0)

    start_all_btn = QPushButton("Start All Animations")
    start_all_btn.setStyleSheet("font-weight: bold; padding: 6px;")
    start_all_btn.clicked.connect(window.start_all_animations)

    stop_all_btn = QPushButton("Stop All Animations")
    stop_all_btn.setStyleSheet("font-weight: bold; padding: 6px;")
    stop_all_btn.clicked.connect(window.stop_all_animations)

    btn_row_layout.addWidget(start_all_btn)
    btn_row_layout.addWidget(stop_all_btn)
    gc_layout.addWidget(btn_row)

    # Global progress bar + label
    window.global_progress_bar = QProgressBar()
    window.global_progress_label = QLabel("0.000s")
    prog_row = QWidget()
    prog_row_layout = QtWidgets.QHBoxLayout(prog_row)
    prog_row_layout.setContentsMargins(0, 0, 0, 0)
    prog_row_layout.addWidget(window.global_progress_bar, stretch=1)
    prog_row_layout.addWidget(window.global_progress_label)
    gc_layout.addWidget(prog_row)

    # Global scrub slider (moves ALL viewers)
    slider_label = QLabel("Global Scrub (all viewers):")
    gc_layout.addWidget(slider_label)
    window.global_progress_slider = QSlider(Qt.Horizontal)
    window.global_progress_slider.setRange(0, 100)
    window.global_progress_slider.setValue(0)
    window.global_progress_slider.valueChanged.connect(window.global_slider_movement)
    gc_layout.addWidget(window.global_progress_slider)

    # Add Viewer button at the very bottom
    add_viewer_btn = QPushButton("+ Add Viewer")
    add_viewer_btn.setStyleSheet("font-weight: bold; padding: 6px;")
    add_viewer_btn.clicked.connect(window.add_view)
    gc_layout.addWidget(add_viewer_btn)

    main_layout.addWidget(global_controls)

    container.show()
    window.setCentralWidget(container)
    window.show()
    app.exec() # Start the main event loop.
