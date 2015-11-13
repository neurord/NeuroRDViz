from __future__ import print_function, division

# Delete from numpy import array, arange, random
import numpy as np
from numpy import array, arange, random
from tvtk.api import tvtk
from mayavi.scripts import mayavi2
import h5py as h5
import sys
import pandas as pd
from mayavi import mlab
from mayavi.sources.vtk_data_source import VTKDataSource

pandasyn=1

#include simData arguments to main once scripts are codependent
def single_type_ug(simData):
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
    
    concentrations = np.arange(grid.shape[0])
    concentrations = getMoleculeConcForEachVoxel(simData)   # Have this pass less through it ********************
    #temperature = np.repeat(colors, 8)

    #velocity = random.randn(points.shape[0], points.shape[1]) # Can show direction of predominant molecule movement at some point
    ug.point_data.scalars = np.repeat(concentrations, 8)
    ug.point_data.scalars.name = 'concentrations'
    # Some vectors.
    #ug1.point_data.vectors = velocity
    #ug1.point_data.vectors.name = 'velocity'
    
    return ug

## Uncomment this to save the file to a VTK XML file.
##save_xml(ug2, 'file.vtu')

# Now view the data.
def view(dataset):
    """ Open up a mayavi scene and display the dataset in it.
    """
    fig = mlab.figure(bgcolor=(1, 1, 1), fgcolor=(0, 0, 0),
                      figure=dataset.class_name[3:])
    surf = mlab.pipeline.surface(dataset, opacity=0.1)
    mlab.pipeline.surface(mlab.pipeline.extract_edges(surf),
                            color=(0, 0, 0), )

def getMoleculeConcForEachVoxel(simData):
    #Function takes data file and returns array concentrations 
    #Function takes array of functions and creates ug
    #Function takes morphology the h5
    #Function that displays unstructured grid 

    molNum = 0 #Change this to be fed in, for now 0 = glu
    '''Must make the following generic for instances with A. Varying sizes. B.No Soma  etc..''' #################################<=======
    dendSnapshot = simData['trial0']['simulation']['dend']['concentrations'][0,:,2] #takes [0'th molecule, all voxel's of it's data, 0'th milisecond)
    somaSnapshot = simData['trial0']['simulation']['soma']['concentrations'][0,:,0] 
    wholeCellSnapshot = np.concatenate((dendSnapshot,somaSnapshot),axis=0)
    return wholeCellSnapshot

@mlab.show
def main(ug):
    view(ug)


if __name__ == '__main__':
    simData = h5.File(sys.argv[1],"r")
    ug = single_type_ug(simData)
    main(ug)
