# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 14:34:04 2023

@author: ts854
"""

import numpy as np
from plyfile import PlyData#, PlyElement
from vedo import Plotter, Axes, Mesh, show, Video, Points, Text2D
import pandas as pd
from pandas import *
import SPHARM as spharm
from IPython.display import display

def importply(filename):
    '''
    filename = string specifying path to file
    '''

    # Import ply data

    plydata = PlyData.read(filename)

    # Import vertices

    vertices = plydata.elements[0].data
    v_matrix = np.empty((len(vertices),3))
    for i in range(len(vertices)):
        vertex = np.array(vertices[i].tolist())
        v_matrix[i,:] = vertex[:3]

    # OPTIONAL Import faces (needed for plotting of triangulated surface)

    faces = plydata.elements[1].data
    f_matrix = np.empty((len(faces),3))
    for i in range(len(faces)):
        f_matrix[i,:] = (faces[i])[0]

    f_matrix = f_matrix.astype(int)

    return v_matrix, f_matrix
 
def importwrl(filepath):
    '''
    **************************************************************************
    *filepath: path to the file to be imported (including file extension)    *
    *                 *
    * e.g., 'filepath = '/Users/username/mydir/2018-01-24.wrl'               *
    **************************************************************************
    *This function reads a vrml file line by line,                           *
    * parses vertex coordinate, vertex normal, and coordinate index per node.*
    *The output is three lists of numpy arrays containing the vertex         *
    *coordinates, normal vectors and coordinate indices of each surface      *
    **************************************************************************
    * Adapted from a vrml2csv converter by Yi Fan Zhao                       *
    * Original code https://github.com/yifnzhao/vrml2csv                     *
    **************************************************************************
    '''
    
    f=open(filepath)
    
    nodeCount = 0

    coord_list = []
    normalVector_list = []
    coordIndex_list = []

    while 1 :   # skip initial parameters    
      ln=f.readline().split() #python3
      # termination condition:
      #eof = 0
      #while ln == []:
      #       ln=f.readline().split()
      #       eof += 1
      #       if eof > 10:
      #           return

      # 09.01.23
      # Changed termination condition to stop at 'texture'. This means that the last coord
      # and coordIndex arrays (which do not correspond to shapes) are not included
    # before it was the word texture that needed to be added manually at the end the wrl file
      if 'bpTracksViewerInventor' in ln:
          return coord_list, normalVector_list, coordIndex_list

      if (ln !=[]) and (ln[0] == 'point'):
          nodeCount+=1
          coord = []
          print(str(nodeCount))      
          ln[4] = ln[4][:-1] 
          coord.append(ln[2:5]) #first coordinate
          while 1:

            ln = f.readline().split()
            if len(ln) > 2:
                ln[2] = ln[2][:-1] #remove comma
                coord.append(ln[0:3])
            if ln == ['}']:
                #print(coord)
                coord = np.asarray(coord, dtype = np.float64) # Converts to an array of floats
                coord_list.append(coord) # Appends array of coordinates of cell i to list
                break
                
          # get normal
  
          normalVector = []
          f.readline() #normal
          f.readline() #Normal {
          ln = f.readline().split()
          ln[4] = ln[4][:-1] #remove comma
          normalVector.append(ln[2:5])
          while 1:
          
            ln = f.readline().split()
            if len(ln)>2:
                if ln[2].endswith(','):
                     ln[2] = ln[2][:-1] 
                     normalVector.append(ln[0:3])
            if ln == ['}']:
                  normalVector = np.asarray(normalVector, dtype = np.float64)
                  normalVector_list.append(normalVector)
                  break
          # then get coordIndex
          coordIndex = []
          ln = f.readline().split() #first coordIndex 
          coordIndex.append([ln[2][:-1],ln[3][:-1],ln[4][:-1]])
          coordIndex.append([ln[6][:-1],ln[7][:-1],ln[8][:-1]])
          while 1:
            
            ln = f.readline().split()
            if len(ln) > 7:
                  coordIndex.append([ln[0][:-1],ln[1][:-1],ln[2][:-1]])
                  coordIndex.append([ln[4][:-1],ln[5][:-1],ln[6][:-1]])
            if len(ln) == 9:
                coordIndex = np.asarray(coordIndex, dtype = np.float64)
                coordIndex_list.append(coordIndex) 
                break

def vertices_um2pixel(gfp, correction):
    return gfp * np.array(correction)

def plot(v_matrix, f_matrix):

    # v_matrix = matrix with format [[x1 y1 z1] [x2 y2 z2]] specifying the vertices of a shape
    # f_matrix = matrix which specifies indices of v_matrix that should be joined to form triangles for triangulated surface

    mesh = Mesh([v_matrix, f_matrix]).c('gold')

    show(mesh, __doc__, viewup='z').close()

def sliderplot(v_matrix, f_matrix):
    meshlst = []
    for i in range(len(v_matrix)):
        meshlst.append(Mesh([v_matrix[i], f_matrix[i]]))

    def sliderfunc(widget, event):
        t = int(widget.value)
        widget.title = "Timepoint "+str(t)
        # remove the old and add the new shape
        # (no need to render as the slider makes a call to rendering)
        plt.pop().add(objs[t])


    objs = meshlst  # load a list of shapes

    plt = Plotter()
    #plt += Text2D(__doc__, pos="top-center", s=1.2)
    plt += Axes(objs[-1])
    plt += objs[0]
    plt.add_slider3d(
        sliderfunc,
        xmin=0,
        xmax=len(objs) - 1,
        pos1 = [20,80,20],
        pos2 = [50,80,20],
        show_value=False,
        s=0.04,
        rotation=45,
        title="Timepoint 0"
    )
    plt.show(zoom=1.2)
    plt.close()

def animation(v_matrix, f_matrix): # maybe file name should be a variable here??

    video = Video("shapeanimation.mp4", duration=8, backend='opencv')

    meshlst = []
    for i in range(len(v_matrix)):
    	meshlst.append(Mesh([v_matrix[i], f_matrix[i]]))
        
    for i in range(len(v_matrix)):
        mesh = meshlst[i]
        plt = show(mesh, viewup='z', camera={'pos':(0,35,30)}, interactive=0)
        video.add_frame()
        plt.clear()
    video.close()                         # merge all the recorded frames
    plt.interactive().close()
    
def reconstructmesh(vertices):
    '''
    vertices: a numpy array with (x,y,z) coordinates of a 3D surface
    reco (output): reconstructed surface from the vertices (vedo Mesh object)
    '''
    
    recopoints = Points(vertices)
    recopoints.subsample(0.005)
    recomesh = recopoints.reconstruct_surface(dims=100, radius=5).c('tomato')
    
    return recomesh, recopoints

def spharmtrajectory(lmax,v_matrix):

    df_coeffs = pd.DataFrame([])
    #coeffs_list = []
    for i in range(len(v_matrix)):

        # Get coordinates for surface i
        v_matrix_i = v_matrix[i]

        # Get spherical harmonics coefficients for surface i
        grid, coeffs, coeffs_dict = spharm.getshcoefficients(v_matrix_i,lmax)

        #grid_rec, reconstructed_vertices = spharm.reconstructshape(coeffs,10000)

        #coeffs_list.append(coeffs)
        
        # Make a dataframe of the spherical harmonics coefficients for all the surfaces
        coeffs_dict.update({'Timepoint': str(i)})
        df_dict = pd.DataFrame([coeffs_dict])
        df_coeffs = pd.concat([df_coeffs,df_dict], ignore_index=True)
        
    with pd.option_context('display.max_rows', 5, 'display.max_columns', 5):
        display(df_coeffs)
        
    return df_coeffs

def reconstructionplot(v_matrix, f_matrix, reco_vertices, lmax):

    '''
    Plot of the raw surfaces alongside the SPHARM reconstructed surfaces
    
    v_matrix: centred vertices of one raw surface (numpy array)
    f_matrix: faces of one raw surface (numpy array)
    reco_vertices: vertices of the reconstructed surface from SPHARM function
                   (numpy array)
    lmax: SPHARM order
    '''
    
    plot = Plotter(shape=(1,3))
    
    # Make raw surface into Mesh object
    mesh = Mesh([v_matrix, f_matrix])
    
    # Reconstruct a mesh from the SPHARM vertices
    recomesh, recopoints = reconstructmesh(reco_vertices)
    
    plot.at(0).show(mesh, __doc__, viewup='z').show(recopoints, __doc__, viewup='z')
    plot.at(0).add(Text2D('Original shape', pos='top-left', font="Arial", s=1))
    
    plot.at(1).show(recopoints, __doc__, viewup='z')
    plot.at(1).add(Text2D('Reconstruction', pos='top-left', font="Arial", s=1))
    plot.at(1).add(Text2D('lmax = '+str(lmax), pos='bottom-left', font="Arial", s=1))
    
    plot.at(2).show(recomesh, __doc__, viewup='z')
    plot.at(2).add(Text2D('Reconstruction', pos='top-left', font="Arial", s=1))
    plot.at(2).add(Text2D('lmax = '+str(lmax), pos='bottom-left', font="Arial", s=1))

    plot.interactive().close()
    
