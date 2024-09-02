
from dataclasses import dataclass, asdict
import json
import numpy
import numpy as np
import matplotlib.pyplot as plt
from skimage import morphology as skmorpho
import os
from os import listdir
from SPHARM import *
from surfaces import *
from distribAnalysisFunctions import *
from shtools import *
from cytoparam import *
from shparam import *
import vedo
from vedo import Plotter, Axes, Mesh, show, Video, Points, Text2D, Volume
import napari
from scipy import interpolate as spinterp
from vedo import show
from vedo import Mesh, dataurl, Plotter
import math 
import numpy as np
from vtk  import *
import vtkmodules


from sklearn.decomposition import PCA


import vedo.vtkclasses as vtki  # a wrapper for lazy imports

import vedo
from vedo.colors import get_color
from vedo.pointcloud import Points
from vedo.utils import buildPolyData, is_sequence, mag, precision
from vedo.utils import numpy2vtk, vtk2numpy, OperationNode
from vedo.visual import MeshVisual
from skimage import morphology as skmorpho


class Cell:
    # number of cycles already in the data 
    def __init__(self, title, startingTime, timestep, numberOfCycles, xyz_um2pixel, segmentationVertices, segmentationFaces, shape_ind_description,sh_coefs_um, numTimepoints, numPlanes):
        self._title = title
        self._segmentationVertices = segmentationVertices
        self._startingTime = startingTime
        self._timestep = timestep
        self._numberOfCycles = numberOfCycles
        self._tps = [startingTime + timestep*i for i in range(numberOfCycles)]
        self._segmentationFaces = segmentationFaces
        #self._sh_coef_dict_pixel = sh_coef_dict_pixel
       # self._sh_coef_dict_um = sh_coef_dict_um
        self._shape_ind_description = shape_ind_description
        self._numTimepoints = numTimepoints
       # self._sh_coefs_pixel = sh_coefs_pixel
        self._sh_coefs_um = sh_coefs_um
        # put array with discrete tp rather than just number
        self._numPlanes = numPlanes
        self._xyz_um2pixel = xyz_um2pixel
        # postHGF hour 
        # add experimental to do stats later if needed laser poswwe, in or out cell, bin nucleated, divide and what tp 

    def get_dict(self):
        dict_obj = {
            "title": self._title,
            "segmentationVertices": self._segmentationVertices,
            "segmentationFaces": self._segmentationFaces,
            "startingTime": self._startingTime,
            "numberOfCycles": self._numberOfCycles,
            "tps": self._tps,
            "timestep": self._timestep,
          #  "sh_coef_dict_pixel" : self._sh_coef_dict_pixel,
          #  "sh_coef_dict_um" : self._sh_coef_dict_um,
            "shape_ind_description" : self._shape_ind_description,
            "numTimepoints" : self._numTimepoints,
            "numPlanes" : self._numPlanes,
            "xyz_um2pixel" : self._xyz_um2pixel,
           # "sh_coefs_pixel" : self._sh_coefs_pixel,
            "sh_coefs_um" : self._sh_coefs_um
            }

        return dict_obj



def getArray(filename, numPlanes, normalise=True):
    # read tif file and normalise every tps 
    # returns a np.array 
    timelist = []
    l = []
    i = 0
    for image in sorted(os.listdir(filename)):
        if not image.startswith('.'):
            i += 1
            print(i)
            im=np.array(Image.open(filename + '/' + image))
            
            l.append(im)
            if i % numPlanes == 0: 
                print('i am changing t')
                array = np.array(l)
                np.swapaxes(array, 0,2)
                if normalise: 
                    array = (array / array.sum()) 
                timelist.append(array)
                l = []

    return timelist


def add_cell_to_json(cell_you_want_to_add, nameOfFile):
    

    with open(nameOfFile, 'r') as f:
        books_json = f.read()

    books = json.loads(books_json)

    # check that book is not already in books
    book_not_already_in_books = True
    for livre in books:
        if cell_you_want_to_add._title == livre['title']:
            book_not_already_in_books = False
        
    if book_not_already_in_books:
        books.append(cell_you_want_to_add.get_dict())

    with open(nameOfFile, 'w') as f:
        f.write(json.dumps(books))


def create_cell(path, title, xyz_um2pixel, startingTime, timestep, numberOfCycles, num_Timepoints, num_Planes):
    
    currentCell = Cell(title, xyz_um2pixel, startingTime, timestep, numberOfCycles, segmentationVertices = [], segmentationFaces = [], sh_coefs_um = [], shape_ind_description = [], numTimepoints = num_Timepoints, numPlanes = num_Planes)

    # go through cytoparam piepline and create cell object with all the required attributes
    sh_coefsDict_pixel = []
    sh_coefs_pixel = []
    sh_coefsDict_um = []
    sh_coefs_um = []
    shapeIndRepresenation_list = []

    ###### get actin signal from the tiffs 
    os.chdir(path)
    gfp = getArray(currentCell._title, currentCell._numPlanes)

    path = '/Users/baptistevauleon/Desktop/Paluch Lab/actin_distribution_analysis/Github/data/' + currentCell._title + '.wrl'
    vertexCoord, normalVector, surfaceIndices = importwrl(path)

    currentCell._segmentationVertices = [l.tolist() for l in vertexCoord]
    currentCell._segmentationFaces = [l.tolist() for l in surfaceIndices]



    #get list of centroids 
    centroidList = []
    for surface_in_wrlfile in vertexCoord: 
        centroidList.append((surface_in_wrlfile*xyz_um2pixel).mean(axis=0, keepdims=True))

    # go through every tp 
    for tp in range(currentCell._numTimepoints):
        vertexCoord_tp = vertexCoord[tp]

        #### Get sh coefficients in um 
        centroid_tp_um = vertexCoord_tp.mean(axis=0, keepdims=True)
        centroid_tp_um = list(centroid_tp_um[0])
        vertexCoord_tp_um = shapecentroidtoorigin(vertexCoord_tp)

        grid, coeffs_um, coeffs_dict_um = getshcoefficients(vertexCoord_tp_um, lmax=16)
        sh_coefsDict_um.append(coeffs_dict_um)
        sh_coefs_um.append(coeffs_um)



        vertexCoord_tp_pixel = vertexCoord_tp*currentCell._xyz_um2pixel
        centroid_tp_pixel = vertexCoord_tp_pixel.mean(axis=0, keepdims=True)
        centroid_tp_pixel = list(centroid_tp_pixel[0])
        vertexCoord_tp_pixel = shapecentroidtoorigin(vertexCoord_tp_pixel)
        
        
        surfaceIndices_tp = surfaceIndices[tp]
        mesh = Mesh([vertexCoord_tp_pixel, surfaceIndices_tp])
    
        #### create fake nucleus mesh 
        # Define the vertices and faces that make up the mesh
        vertsNuc = vertexCoord_tp_pixel  #/1.2
        cellsNuc = surfaceIndices_tp # cells same as faces
        # Build the polygonal Mesh object from the vertices and faces
        meshNuc = Mesh([vertsNuc, cellsNuc])

        #### Get sh coefficients in pixel
        grid, coeffs_pixel, coeffs_dict_pixel = getshcoefficients(vertexCoord_tp_pixel, lmax=16)
        sh_coefsDict_pixel.append(coeffs_dict_pixel)
        sh_coefs_pixel.append(coeffs_pixel)
        gridNuc, coeffsNuc, coeffs_dictNuc = getshcoefficients(vertsNuc, lmax=16)
    
        #### Reconstruct shape from sh coefficient
        grid_rec, reconstructed_vertices = reconstructshape(coeffs_pixel, npoints= 10000)

        # Run the cellular mapping to create a parameterized intensity representation
        # for the FP image:
        gfp_representation = cellular_mapping(
            coeffs_mem=coeffs_dict_pixel,
            centroid_mem=centroid_tp_pixel,
            coeffs_nuc=coeffs_dictNuc,
            centroid_nuc=centroid_tp_pixel,
            nisos=[1000, 0],
            images_to_probe=[('gfp', gfp[tp])]
        ).data.squeeze()

        shapeIndRepresenation_list.append(gfp_representation.tolist())
  
   # currentCell._sh_coef_dict_pixel = sh_coefsDict_pixel
   # currentCell._sh_coefs_pixel = [l.tolist() for l in sh_coefs_pixel]

   # currentCell._sh_coef_dict_um = sh_coefsDict_um
    currentCell._sh_coefs_um = [l.tolist() for l in sh_coefs_um]

    currentCell._shape_ind_description =  shapeIndRepresenation_list

    return currentCell


def read_one_cell_in_json(cellId, jsonFileName):

    with open(jsonFileName, 'r') as f:
        cells_json = f.read()

    cells = json.loads(cells_json)
    for i in range(len(cells)):
        if cells[i]['title'] == cellId:
            print('cell: ' + cellId + ' was found')
            new_cell = Cell(title = cells[i]['title'], shape_ind_description = cells[i]['shape_ind_description'], sh_coefs_um = cells[i]['sh_coefs_um'], startingTime = cells[i]['startingTime'], timestep = cells[i]['timestep'], numberOfCycles = cells[i]['numberOfCycles'], xyz_um2pixel =  cells[i]['xyz_um2pixel'], numTimepoints = cells[i]['numTimepoints'], segmentationVertices = cells[i]['segmentationVertices'], numPlanes = cells[i]['numPlanes'], segmentationFaces = cells[i]['segmentationFaces'])

    return new_cell



def read_all_cells_in_json_tp(jsonFileName, wantedTp):
    # not working for the moment
    
    # imprvement idea make it bigger than tps and allow for any category selection 
    listOfcellswithTps = []
    
    with open(jsonFileName, 'r') as f:
        cells_json = f.read()

    cells = json.loads(cells_json)
    for i in range(len(cells)):
        if wantedTp in cells[i]['tps']:
            print('cell:  has the selected timepoint')
            new_cell = Cell(title = cells[i]['title'], startingTime = cells[i]['startingTime'], sh_coefs_um = cells[i]['sh_coefs_um'][cells[i]['tps'].index(wantedTp)], timestep = cells[i]['timestep'], numberOfCycles = cells[i]['numberOfCycles'], shape_ind_description = cells[i]['shape_ind_description'][cells[i]['tps'].index(wantedTp)],  xyz_um2pixel =  cells[i]['xyz_um2pixel'], numTimepoints = cells[i]['numTimepoints'], segmentationVertices = cells[i]['segmentationVertices'][cells[i]['tps'].index(wantedTp)], numPlanes = cells[i]['numPlanes'], segmentationFaces = cells[i]['segmentationFaces'][cells[i]['tps'].index(wantedTp)])
            listOfcellswithTps.append(new_cell)
    
    return listOfcellswithTps



def reconstruct_onceCell_oneTp_um(cellTitle, jsontitle, tp):
    
    new_cell = read_one_cell_in_json(cellTitle, jsontitle)
    # # Reconstruct the shape from coefficients
    grid_rec, reconstructed_vertices_tp = spharm.reconstructshape(np.array(new_cell._sh_coefs_um[tp]),10000)

    # # Reconstruct the mesh from spherical harmonics coefficients and plot it
    reconstructionplot(spharm.shapecentroidtoorigin(np.array(new_cell._segmentationVertices[tp])), new_cell._segmentationFaces[tp], reconstructed_vertices_tp, lmax=16)


def read_one_cell_specific_tp_in_json():
    return 

def project_shapeInd_description_on_shape():
    return



def plot_histogram_array():
    return 


def visualise_cell(cellId, jsonFileName, tp):
    # plot gfp signal and segmented original mesh 
    # plot reconstruction based on sh coeff
    # plot shape independant description of the actin distribution on given shape 
    
    new_cell = read_one_cell_in_json(cellId, jsonFileName)
    gfp_raw = getArray(new_cell._title, new_cell._numPlanes, normalise=False)

    v_matrix = [i*np.array(new_cell._xyz_um2pixel) for i in new_cell._segmentationVertices]
    v_matrix = [ np.array(i) for i in v_matrix]
    v_matrix_swapped = []
    for i in v_matrix:
        i[:, [2, 0]] = i[:, [0, 2]]
        v_matrix_swapped.append(i)
    faces_matrix = [ np.array(i).astype('int32') for i in new_cell._segmentationFaces]
    vertices = v_matrix[tp] #np.array([[0, 0], [0, 20], [10, 0], [10, 10]])
    faces = faces_matrix[tp] #np.array([[0, 1, 2], [1, 2, 3]])
    values = np.linspace(0, 1, len(vertices))
    surface = (vertices, faces, values)
    viewer = napari.view_surface(surface)
    layer = viewer.add_image(gfp_raw)
  #  layer.scale = (1, 1, 1)  # now it will render correctly
    napari.run()
    
    

def vtkImageData_to_numpy(vtk_image_data):
    # Get the dimensions of the vtkImageData
    dims = vtk_image_data.GetDimensions()
    
    # Get the VTK data array from the vtkImageData object
    vtk_array = vtk_image_data.GetPointData().GetScalars()
    
    # Convert the VTK array to a NumPy array
    np_array = vtkmodules.util.numpy_support.vtk_to_numpy(vtk_array)
    
    # Reshape the NumPy array to the dimensions of the vtkImageData
    np_array = np_array.reshape(dims[2], dims[1], dims[0])
    
    return np_array

def binarize(
        self,
        values=(255, 0),
        spacing=None,
        dims=None,
        origin=None,
    ):
        """
        Convert a `Mesh` into a `Volume` where
        the interior voxels value is set to `values[0]` (255 by default), while
        the exterior voxels value is set to `values[1]` (0 by default).

        Arguments:
            values : (list)
                background and foreground values.
            spacing : (list)
                voxel spacing in x, y and z.
            dims : (list)
                dimensions (nr. of voxels) of the output volume.
            origin : (list)
                position in space of the (0,0,0) voxel.

        Examples:
            - [mesh2volume.py](https://github.com/marcomusy/vedo/tree/master/examples/volumetric/mesh2volume.py)

                ![](https://vedo.embl.es/images/volumetric/mesh2volume.png)
        """
        assert len(values) == 2, "values must be a list of 2 values"
        fg_value, bg_value = values

        bounds = self.bounds()
        if spacing is None:  # compute spacing
            spacing = [0, 0, 0]
            diagonal = np.sqrt(
                  (bounds[1] - bounds[0]) ** 2
                + (bounds[3] - bounds[2]) ** 2
                + (bounds[5] - bounds[4]) ** 2
            )
            spacing[0] = spacing[1] = spacing[2] = diagonal / 1200

        if dims is None:  # compute dimensions
            dim = [0, 0, 0]
            for i in [0, 1, 2]:
                dim[i] = int(np.ceil((bounds[i*2+1] - bounds[i*2]) / spacing[i]))
        else: 
            dim = dims 
        white_img = vtki.vtkImageData()
        white_img.SetDimensions(dim)
        white_img.SetSpacing(spacing)
        white_img.SetExtent(0, dim[0]-1, 0, dim[1]-1, 0, dim[2]-1)

        if origin is None:
            origin = [0, 0, 0]
            origin[0] = bounds[0] + spacing[0]
            origin[1] = bounds[2] + spacing[1]
            origin[2] = bounds[4] + spacing[2]
        white_img.SetOrigin(origin)

        # if direction_matrix is not None:
        #     white_img.SetDirectionMatrix(direction_matrix)

        white_img.AllocateScalars(vtki.VTK_UNSIGNED_CHAR, 1)

        # fill the image with foreground voxels:
        white_img.GetPointData().GetScalars().Fill(fg_value)

        # polygonal data --> image stencil:
        pol2stenc = vtki.new("PolyDataToImageStencil")
        pol2stenc.SetInputData(self.dataset)
        pol2stenc.SetOutputOrigin(white_img.GetOrigin())
        pol2stenc.SetOutputSpacing(white_img.GetSpacing())
        pol2stenc.SetOutputWholeExtent(white_img.GetExtent())
        pol2stenc.Update()

        # cut the corresponding white image and set the background:
        imgstenc = vtki.new("ImageStencil")
        imgstenc.SetInputData(white_img)
        imgstenc.SetStencilConnection(pol2stenc.GetOutputPort())
        # imgstenc.SetReverseStencil(True)
        imgstenc.SetBackgroundValue(bg_value)
        imgstenc.Update()

        vol = vedo.Volume(imgstenc.GetOutput())
        vol.name = "BinarizedVolume"
        vol.pipeline = OperationNode(
            "binarize",
            parents=[self],
            comment=f"dims={tuple(vol.dimensions())}",
            c="#e9c46a:#0096c7",
        )
        return vol

def project_distribution_on_shape(binary_array_of_shape, shape_ind_description): 
    mem_round = nuc_round = binary_array_of_shape
    gfp_representation = shape_ind_description

    coords_round, _ = parameterize_image_coordinates(
        seg_mem=mem_round,
        seg_nuc=nuc_round,
        lmax=16,
        nisos=[1000, 0])

    gfp_morphed = morph_representation_on_shape(
        img= mem_round + nuc_round,
        param_img_coords=coords_round,
        representation=gfp_representation)


    viewer = napari.Viewer()
    layer = viewer.add_image(gfp_morphed)
    layer = viewer.add_image(binary_array_of_shape)
    layer.scale = (1, 1, 1)  # now it will render correctly
    napari.run()

    return gfp_morphed

def get_binary_array_from_vertices(v_matrix_tp, faces):

    v_matrix = np.array(v_matrix_tp) 
    v_matrix = spharm.shapecentroidtoorigin(v_matrix) 
    
    mesh_mem = Mesh([v_matrix , faces])

    # once we have the mesh, get the binarized array 
    # we choose to have an 100% nucleus and 0% cytoplasm 
    mesh_mem_bin = binarize(self = mesh_mem)
    mesh_mem_vtk = Volume.imagedata(mesh_mem_bin)
    mesh_mem_array = vtkImageData_to_numpy(mesh_mem_vtk)



    # # # # get the shape we want to project our distribution on, specifically get the sh coefficients 
# # # w = 300
# # # mem_round = np.zeros((w, w, w), dtype = np.uint8)
# # # mem_round[int(w*0.2):int(w*0.8), int(w*0.2):int(w*0.8), int(w*0.2):int(w*0.8)] = 1

# # # nuc_round = np.zeros((w, w, w), dtype = np.uint8)
# # # nuc_round[int(w*0.4):int(w*0.6), int(w*0.4):int(w*0.6), int(w*0.4):int(w*0.6)] = 1


# w = 300
# mem_round = skmorpho.ball(w // 3) # radius of our round cell
# nuc_round = mem_round

    return mesh_mem_array

def reconstruct_shapeInd_array(raw_list_shapeInd, numberOfInterpolatedPlanes):
    shapeIndDescription = np.array([raw_list_shapeInd])

    return shapeIndDescription

def save_morphed_signa_tiffs(folderName, binaryMask, shapeInddescription):
    currentPath = os.getcwd()
    newpath = currentPath + '/' + str(folderName)
    os.makedirs(newpath)
    os.chdir(newpath)

    # store morphed distribution 
    os.makedirs(newpath + '/morphed distriubtion stack')
    os.chdir(newpath + '/morphed distriubtion stack')

    i =0 
    for plane in shapeInddescription: 
        i +=1 
        d = np.array(plane)

   # im = Image.fromarray(d, mode='F') # float32
        im = Image.fromarray(d) # float32
        im.save(folderName + 'morphed distriubtion stack' + str(i) + ".tiff", "TIFF")


    # store shape binary mask 

    os.makedirs(newpath + '/binary mask stack')
    os.chdir(newpath + '/binary mask stack')

    i =0 
    for plane in binaryMask: 
        i +=1 
        d = np.array(plane * 255)

   # im = Image.fromarray(d, mode='F') # float32
        im = Image.fromarray(d) # float32
        im.save(folderName + 'binary mask stack' + str(i) + ".tiff", "TIFF")
    
path = '/Users/baptistevauleon/Desktop/Paluch Lab/actin_distribution_analysis/Github/data'
os.chdir(path)

new_cell = read_one_cell_in_json('noHGFstack4', 'your.json')

visualise_cell('noHGFstack4', 'your.json', tp=0)

# # # # # get the shape ind description we want to project on a choosen shape 
# # gfp_representation_fake = np.zeros([1001, 8194])
# # gfp_representation_fake[:,2:4] = 1 
# # gfp_representation = gfp_representation_fake
# # # shape_ind_description = np.array(new_cell._shape_ind_description[0])

plt.imshow(new_cell._shape_ind_description[0], aspect='auto')
plt.show()
# print(new_cell._shape_ind_description[0])

# # # get mesh of shape we want to project gfp_morphed on
# # v_matrix_tp = new_cell._segmentationVertices[0]
# # faces = new_cell._segmentationFaces[0]
# #binary_array_of_shape = get_binary_array_from_vertices(v_matrix_tp, faces)

#     # # # # get the shape we want to project our distribution on, specifically get the sh coefficients 
# w = 300
# # mem_round = np.zeros((w, w, w), dtype = np.uint8)
# # mem_round[int(w*0.2):int(w*0.8), int(w*0.2):int(w*0.8), int(w*0.2):int(w*0.8)] = 1
# mem_round = skmorpho.ball(w // 3)

# # nuc_round = np.zeros((w, w, w), dtype = np.uint8)
# # nuc_round[int(w*0.4):int(w*0.6), int(w*0.4):int(w*0.6), int(w*0.4):int(w*0.6)] = 1


# # w = 300
# # mem_round = skmorpho.ball(w // 3) # radius of our round cell
# # nuc_round = mem_round

# project_distribution_on_shape(mem_round, gfp_representation)