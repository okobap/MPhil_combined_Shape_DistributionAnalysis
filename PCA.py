
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import pandas as pd
import os 
from cellClass import *
from vedo import show
from vedo import Mesh, dataurl, Plotter
import math 
import numpy as np
from vtk  import *
import vtkmodules
from cellClass import *
from cellClass import read_one_cell_in_json



### EXAMPLE

# # Get PCA coefficients and plot respective contributions to variance

# pca = PCA(n_components=2)
# trans = pca.fit_transform(df_coeffs.drop(columns=['Timepoint']))
# df_trans = pd.DataFrame(trans)
# df_trans.columns = ['PC1', 'PC2']
# df_trans['Timepoint'] = df_coeffs.Timepoint

# with pd.option_context('display.max_rows', 5, 'display.max_columns', 5):
#     display(df_trans)

# variance_ratio = pca.explained_variance_ratio_
# print(variance_ratio)

# # PCA plot for a single trajectory

# timepoints = [float(i) for i in df_trans.Timepoint.tolist()]

# fig = plt.figure(figsize=(8,5))
# ax = fig.subplots(1,1)
# sc = ax.scatter(df_trans.PC1, df_trans.PC2, s=50, c=timepoints, cmap = 'viridis', zorder=2) # zorder is order of layering points

# # If joining up with lines need to make sure that PC1,2 columns are ordered by increasing time

# ax.plot(df_trans.PC1, df_trans.PC2, 'k', zorder=1)
# plt.xlabel('PC1 (80%)',fontsize=25)
# plt.ylabel('PC2 (7.8%)',fontsize=25)
# plt.xticks(fontsize=25)
# plt.yticks(fontsize=25)
# cbar = plt.colorbar(sc)
# cbar.ax.tick_params(labelsize=20)
# #cbar.ax.set_ylabel("Timepoint", rotation=270,fontsize=25)
# plt.show()

### GUIDELINES ON HOW TO USE PCA MODULE 
# data should be in a dataframe 



### END OF THE GUIDELINES




def center_df_data(df_data):
    # take a dataframe as input and substracts the mean of each column to every value in this column for every columns of the dataframe
    # input 
    #       df_data: dataframe
    #
    # returns data - <data> and <data>

    # mean = df_data.loc[:,df_data.columns != "TimePoints"].mean()
    mean = df_data.mean()
    df_centered_data = df_data - mean 

    return mean, df_centered_data

def create_PC_space(n_components, df_centered_data, ):
    # create basis for PCA to visualise cells on 
    #  input: 
    #   df_centered_data: dataframe including the data the PCA will be done on, should be centered using the center_df_data
    #   n_components: number of PCs wanted 
    #   
    #
    # returns: 
    #   pca: the pca object from class PCA fitted with the data givin as input


    pca = PCA(n_components)
    pca_fitted = pca.fit(df_centered_data)
    return pca_fitted

def project_on_PC_space(pca_fitted, df_centered_data, columns_names):

    # columns_names: names of the columns in the returned dataframe, should be of size n_components, example: ['PC1', 'PC2']
    # df_data_projected: dataframe containing the projected df_data

    data_projected = pca_fitted.transform(df_centered_data)
    df_data_projected = pd.DataFrame(data_projected)
    df_data_projected.columns = columns_names

    return df_data_projected

def PC_percentages_plot():
   # plot poucentegaes of PC check main with line 81
    return

def plot_single_trajectory_on_PC_space(df_projected_data):
    # project the trajectory of a single cell on the PC space, for now needs columns named PC1 and PC2
    timepoints = [float(i) for i in df_projected_data.Timepoint.tolist()]
    fig = plt.figure(figsize=(8,5))
    ax = fig.subplots(1,1)
    sc = ax.scatter(df_projected_data.PC1, df_projected_data.PC2, s=50, c=timepoints, cmap = 'viridis', zorder=2) # zorder is order of layering points
    ax.plot(df_projected_data.PC1, df_projected_data.PC2, 'k', zorder=1)
    plt.xlabel('PC1',fontsize=25)
    plt.ylabel('PC2',fontsize=25)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    cbar = plt.colorbar(sc)
    cbar.ax.tick_params(labelsize=20)
    #cbar.ax.set_ylabel("Timepoint", rotation=270,fontsize=25)
    plt.show()

def plot_single_points_on_PC_space(df_projected_data):
    # project the trajectory of a single cell on the PC space, for now needs columns named PC1 and PC2
    timepoints = [float(i) for i in df_projected_data.Timepoint.tolist()]
    fig = plt.figure(figsize=(8,5))
    ax = fig.subplots(1,1)
    sc = ax.scatter(df_projected_data.PC1, df_projected_data.PC2, s=50, c=timepoints, cmap = 'viridis', zorder=2) # zorder is order of layering points
    #ax.plot(df_projected_data.PC1, df_projected_data.PC2, 'k', zorder=1)
    plt.xlabel('PC1',fontsize=25)
    plt.ylabel('PC2',fontsize=25)
    plt.xticks(fontsize=25)
    plt.yticks(fontsize=25)
    cbar = plt.colorbar(sc)
    cbar.ax.tick_params(labelsize=20)
    #cbar.ax.set_ylabel("Timepoint", rotation=270,fontsize=25)
    plt.show()

def get_shape_df_from_cellId(cellId, jsonFileName, align): 
    df_coeffs_cells = []
    new_cell = read_one_cell_in_json(cellId, jsonFileName)
    v_matrix = [ np.array(i) for i in new_cell._segmentationVertices]
    v_matrix = [ spharm.shapecentroidtoorigin(i) for i in v_matrix]
    recomesh, recopoints = reconstructmesh(v_matrix[0])
    show(recomesh, __doc__, axes=1).close()
    plt.show()
   # print(v_matrix)
    temp = []
    if align: 
        for single_shape in v_matrix:
            
            centered_df = pd.DataFrame(data=single_shape)
            pca_fitted = create_PC_space(n_components=2, df_centered_data = centered_df)
            PC1coordinates = pca_fitted.inverse_transform(np.array([1,0]))
        #    print(PC1coordinates)
            angle = math.atan(PC1coordinates[0]/PC1coordinates[1])

            rotated_v_matrix = [rotate([0,0,0], point, angle) for point in single_shape]
            temp.append(np.array(rotated_v_matrix))
        v_matrix = temp

            # centered_df = pd.DataFrame(data=v_matrix[0])
            # pca_fitted = create_PC_space(n_components=2, df_centered_data = centered_df)
            # PC1coordinates = pca_fitted.inverse_transform(np.array([1,0]))
            # print(PC1coordinates)
            # angle = math.atan(PC1coordinates[0]/PC1coordinates[1])
            # rotated_v_matrix = [rotate([0,0,0], point, angle) for point in v_matrix[0]]
    x = [l[0] for l in v_matrix[0]]
    y = [l[1] for l in v_matrix[0]]
    z = [l[2] for l in v_matrix[0]]
    fig = plt.figure(figsize = (10,10))
    ax = plt.axes(projection='3d')
    ax.grid()
    ax.scatter(x, y, z, c = 'r', s = 5)
    plt.show()
    recomesh, recopoints = reconstructmesh(v_matrix[0])
    show(recomesh, __doc__, axes=1).close()
    plt.show()
    df_coeffs = spharmtrajectory(25,v_matrix)
    df_withoutTp = df_coeffs.drop(columns=['Timepoint'])

    return df_coeffs, df_withoutTp, v_matrix

def get_shapeInd_df_from_cellId(cellId, jsonFileName, align): 
    new_cell = read_one_cell_in_json(cellId, jsonFileName)
    shape_ind_descriptions = []
    df_shape_ind_descriptions = []
    v_matrix = [ np.array(i) for i in new_cell._segmentationVertices]
    v_matrix = [ spharm.shapecentroidtoorigin(i) for i in v_matrix]

    for tp in range(new_cell._numTimepoints):
        shape_ind_description_tp = new_cell._shape_ind_description[tp]
        print('i am in the tp loop')
        shape_ind_description_tp_singleRaw = []
    
        if align: 
            single_shape = v_matrix[tp]
            centered_df = pd.DataFrame(data=single_shape)
            pca_fitted = create_PC_space(n_components=2, df_centered_data = centered_df)
            PC1coordinates = pca_fitted.inverse_transform(np.array([1,0]))
            print(PC1coordinates)
            angle = math.atan(PC1coordinates[0]/PC1coordinates[1])
            print(angle)
            angle = int(angle * 127 / (2*(math.pi)))
            angle = (angle+1) * 64 
            print(angle)
            indexList = [i for i in range(8194)]
            indexListPart3 = indexList[2:(8194-angle)]
            indexListPart2 = indexList[(8194-angle):]
            indexListPart1 = [0,1]
            indexList = indexListPart1 + indexListPart2 + indexListPart3
            
            gfp_representation_unrotated = np.array(shape_ind_description_tp)
            gfp_representation_rotated = gfp_representation_unrotated[:, indexList]
            
            shape_ind_description_tp = [i.tolist() for i in gfp_representation_rotated ] 
            print('i am in the align loop')
        shape_ind_descriptions.append(shape_ind_description_tp)

        for shape_ind_description_lign in shape_ind_description_tp:
            shape_ind_description_tp_singleRaw += shape_ind_description_lign

        df_shape_ind_description_tp = pd.DataFrame(shape_ind_description_tp_singleRaw).transpose()
        df_shape_ind_descriptions.append(df_shape_ind_description_tp)

    df_shape_ind_descriptions  = pd.concat(df_shape_ind_descriptions)
    df_shape_ind_descriptions.insert(0,"Timepoint",  [i for i in range(new_cell._numTimepoints)])
    df_withoutTp = df_shape_ind_descriptions.drop(columns=['Timepoint'])

    return df_shape_ind_descriptions, df_withoutTp, shape_ind_descriptions

def plot_shape_from_PCspace_coords(pca_fitted, coordinates, lmaximum, mean):

    meanPCA_raw_uncentered = pca_fitted.inverse_transform(np.array(coordinates))

    meanPCA_raw = [meanPCA_raw_uncentered[i] + list(mean)[i] for i in range(len(meanPCA_raw_uncentered))]

    meanPCA = []
    tempi = []
    tempj = []
    i = 0
    for coefficient in meanPCA_raw: 
        i += 1
        tempi.append(coefficient)
        if i % (lmaximum+1) == 0:
            tempj.append(tempi)
            tempi =[]
        if i % (lmaximum)**2 == 0: 
            meanPCA.append(tempj)
            tempj = []

# # # Reconstruct the shape from coefficients
    grid_rec, reconstructed_vertices_tp = spharm.reconstructshape(np.array(meanPCA),10000)
    # print(reconstructed_vertices_tp)
    # x = [l[0] for l in reconstructed_vertices_tp]
    # y = [l[1] for l in reconstructed_vertices_tp]
    # z = [l[2] for l in reconstructed_vertices_tp]
    # fig = plt.figure(figsize = (10,10))
    # ax = plt.axes(projection='3d')
    # ax.grid()
    # ax.scatter(x, y, z, c = 'r', s = 5)
    # plt.show()


# # Reconstruct the mesh from spherical harmonics coefficients and plot it
    recomesh, recopoints = reconstructmesh(reconstructed_vertices_tp)

    show(recomesh, __doc__, axes=1).close()
    plt.show()

    return meanPCA, recomesh

def save_distribution_fromPC_coordinates(pca_fittedShapeInd, coords, meanShapeInd, binary_array_of_shape, title):
# #project PC distribution coordinate on mean shape of the PC space 
    print('project mean distribution on mean shape of the PC space') 
    meanDistribPCA_raw_uncentered = np.array(pca_fittedShapeInd.inverse_transform(np.array(coords))) 
    print(meanDistribPCA_raw_uncentered)
    print('i a m just after the inverse transofrm')
    #meanDistribPCA_raw = [meanDistribPCA_raw_uncentered[i] + list(meanShapeInd)[i] for i in range(len(meanDistribPCA_raw_uncentered))]
    meanDistribPCA_raw = meanDistribPCA_raw_uncentered + np.array(meanShapeInd)# np.array(meanDistribPCA_raw)
    print('i am after the decentering of the mean distrib')
    meanDistribPCA = meanDistribPCA_raw.reshape(1001,8194)

    print('i am after the array reconstruction')
    meanDistrib = np.array(meanDistribPCA)


   # binary_array_of_shape = get_binary_array_from_vertices(v_matrix_tp, faces)
    gfp_morphed = project_distribution_on_shape(binary_array_of_shape, meanDistrib)

    save_morphed_signa_tiffs(title, binary_array_of_shape, gfp_morphed)



def plot_mean_shape_from_df_shcoeffs(df_shcoeffs, lmaximum):

    

    meanShape_raw = list(df_shcoeffs.mean())


    meanShape = []
    tempi = []
    tempj = []
    i = 0
    for coefficient in meanShape_raw: 
        i += 1
        tempi.append(coefficient)
        if i % (lmaximum+1) == 0:
            tempj.append(tempi)
            tempi =[]
        if i % (lmaximum)**2 == 0: 
            meanShape.append(tempj)
            tempj = []

# # # Reconstruct the shape from coefficients
    grid_rec, reconstructed_vertices_tp = spharm.reconstructshape(np.array(meanShape),10000)
    print(reconstructed_vertices_tp)
    x = [l[0] for l in reconstructed_vertices_tp]
    y = [l[1] for l in reconstructed_vertices_tp]
    z = [l[2] for l in reconstructed_vertices_tp]
    fig = plt.figure(figsize = (10,10))
    ax = plt.axes(projection='3d')
    ax.grid()
    ax.scatter(x, y, z, c = 'r', s = 5)
    plt.show()


# # Reconstruct the mesh from spherical harmonics coefficients and plot it
    recomesh, recopoints = reconstructmesh(reconstructed_vertices_tp)
    show(recomesh, __doc__, axes=1).close()
    plt.show()

    return meanShape, recomesh


def rotate(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    ox = origin[0]
    oy = origin[1]
    px = point[0]
    py = point[1]

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return [qx,qy, point[2]]

def plot_vertices(cellId, jsonFileName): 
    new_cell = read_one_cell_in_json(cellId, jsonFileName)
    v_matrix = [ np.array(i) for i in new_cell._segmentationVertices]
    v_matrix = [ spharm.shapecentroidtoorigin(i) for i in v_matrix]
    x = [l[0] for l in v_matrix[0]]
    y = [l[1] for l in v_matrix[0]]
    z = [l[2] for l in v_matrix[0]]
    fig = plt.figure(figsize = (10,10))
    ax = plt.axes(projection='3d')
    ax.grid()
    ax.set_xlim(-20,20)
    ax.set_ylim(-20,20)
    ax.scatter(x, y, z, c = 'r', s = 5)
    plt.show()






# # # main under construction 


# path = '/Users/baptistevauleon/Desktop/Paluch Lab/actin_distribution_analysis/Github/data'
# os.chdir(path)

# # # # # gather cells to create the PC space
# df_coeffs_cells = []
# # # #names_cells_of_interest = ['noHGFstack1', 'noHGFstack2', 'noHGFstack5', 'postEMTcell1' , 'postEMTcell2', 'postEMTcell3' ]
# names_cells_of_interest = [ 'postEMTcell1' ]
# for cellname in names_cells_of_interest:
#     df_coeffs_PCspace, df_withoutTp_PCspace, rotated_v_matrix_PCspace = get_shape_df_from_cellId(cellname, 'your.json', align=True)
#     df_coeffs_cells.append(df_withoutTp_PCspace)
# df_withoutTp_PCspace = pd.concat(df_coeffs_cells)

# # df_shape_ind_descriptions, df_withoutTp, gfp_representation = get_shapeInd_df_from_cellId('noHGFstack2', 'your.json', align=True)
# # df_shape, df_withoutTp, rotated_v_matrix = get_shape_df_from_cellId('noHGFstack2', 'your.json', align=True)

# #plot_mean_shape_from_df_shcoeffs(df_withoutTp_PCspace, lmaximum=25)



# # # # center the dataframe used to fit the PC space
# # mean, centered_df = center_df_data(df_withoutTp_PCspace)


# # # # create the PCs space
# # pca_fitted = create_PC_space(n_components=2, df_centered_data = centered_df)


# # # # # plot mean PCA shape 
# # meanPCA, reconstructed_mesh = plot_shape_from_PCspace_coords(pca_fitted, [0,0], 25, mean)
# # # missing plot of more than 1 set of coordinates inside the PC space


# # # # # print the variance ratios of the PCs 
# # variance_ratio = pca_fitted.explained_variance_ratio_
# # print(variance_ratio)


# # # # # gather cells to project 
# # df_coeffs_toProject, df_withoutTp_toProject, rotated_v_matrix_toProject = get_shape_df_from_cellId('20240529-sample1-cell2', 'your.json')


# # # #project data we want to visualise on the PC space 
# # df_projected_data = project_on_PC_space(pca_fitted, df_centered_data = df_withoutTp_toProject - mean, columns_names=['PC1', 'PC2'])
# # #df_projected_data = project_on_PC_space(pca_fitted, df_centered_data = centered_df, columns_names=['PC1', 'PC2'])
# # #df_projected_data.insert(0,"Timepoint", list(df_coeffs_toProject.Timepoint) )
# # df_projected_data.insert(0,"Timepoint", [0,1,2,100, 200, 300] )


# # # # plot projected data in projected space
# # #plot_single_trajectory_on_PC_space(df_projected_data)
# # plot_single_points_on_PC_space(df_projected_data)