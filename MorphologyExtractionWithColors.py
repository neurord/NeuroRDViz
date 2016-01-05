from __future__ import print_function, division

# Delete from numpy import array, arange, random
import numpy as np
from tvtk.api import tvtk
from mayavi.scripts import mayavi2
import h5py as h5
import sys
import pandas as pd
from mayavi import mlab
from mayavi.sources.vtk_data_source import VTKDataSource


#include simData arguments to main once scripts are codependent
def unstructuredGridMorpho(simData, molnum):
    """A slightly more complex example of how to generate an
    unstructured grid with different cell types.  Returns a created
    unstructured grid.
    """    

    #Extracting points, WHAT DOES view DO???  What package is it in.    
    grid = np.array(simData['trial0']['model']['grid']).view(np.recarray)
    print(grid)
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
    # FOR DEBUGGING.  DELETE AND GET RID OF PANDAS 
    df = pd.DataFrame(grid)
    df['our_volume'] = (xmax-xmin) * (zmax-zmin) * (zmax-zmin)
    print(df)
    

   
    voxels = np.arange(points.shape[0]).reshape(-1, 8)

    voxel_type = tvtk.Hexahedron().cell_type # VTK_HEXAHEDRON == 12
    ug = tvtk.UnstructuredGrid(points=points)
    ug.set_cells(voxel_type, voxels)
    
    concentrations = getMoleculeConcForEachVoxel(simData, 0, molnum)   # Have this pass less through it ********************
    #temperature = np.repeat(colors, 8)

    ug.point_data.scalars = np.repeat(concentrations, 8)
    ug.point_data.scalars.name = 'concentrations'
    # Some vectors.
    #ug1.point_data.vectors = velocity
    #ug1.point_data.vectors.name = 'velocity'
    
    return ug

## Uncomment this to save the file to a VTK XML file.
##save_xml(ug2, 'file.vtu')

# Now view the data.
@mlab.show
def view(ug):
    """ Open up a mayavi scene and display the dataset in it.
    """
    fig = mlab.figure(bgcolor=(1, 1, 1), fgcolor=(0, 0, 0),)
                      
    surf = mlab.pipeline.surface(ug, opacity=0.1)
    
    mlab.pipeline.surface(mlab.pipeline.extract_edges(surf),
                            color=(0, 0, 0), )

        
@mlab.animate(delay=10, ui=False)
def anim(ug):
    f = mlab.gcf()
    
    for i in range(len(simData['trial0']['simulation']['dend']['times'])):
        try:
            concentrations = getMoleculeConcForEachVoxel(simData, i, molnum)
            ug.point_data.scalars = np.repeat(concentrations, 8) #make 8 to static variable of "voxelpts"
            surf = mlab.pipeline.surface(ug, opacity=0.1)            
            mlab.pipeline.surface(mlab.pipeline.extract_edges(surf),
                            color=(0, 0, 0), )
            #print(concentrations)
            #f.scene.camera.azimuth(10)  #Rotates camera by 10
            f.scene.render()
        except AttributeError:
            pass
        yield
        
        #Increase Start/stop/delay window size
        #Adjust default - x10
      


def getMoleculeConcForEachVoxel(simData, ms, molnum):
    #Function takes data file and returns array concentrations 
    #Function takes array of functions and creates ug
    #Function takes morphology the h5
    #Function that displays unstructured grid 

    '''Must make the following generic for instances with A. Varying sizes. B.No Soma  etc..''' #################################<=======
    dendSnapshot = simData['trial0']['simulation']['dend']['concentrations'][ms,:,molnum] #takes [milisecond, all voxel's of it's data, molecule of interest)
    somaSnapshot = simData['trial0']['simulation']['soma']['concentrations'][ms,:,molnum] 
    wholeCellSnapshot = np.concatenate((dendSnapshot,somaSnapshot),axis=0)
    return wholeCellSnapshot


if __name__ == '__main__':
    simData = h5.File(sys.argv[1],"r")
    molname=sys.argv[2]
    molnum=0 #place holder for function that converts molname to num
    ug = unstructuredGridMorpho(simData, molnum)
    a = anim(ug) # Starts the animation.
    view(ug)