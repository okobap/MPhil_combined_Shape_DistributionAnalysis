from cellClass import * 
from PCA import * 
from vedo import *
from pandas import *





path = '/Users/baptistevauleon/Desktop/Paluch Lab/actin_distribution_analysis/Github/data'
os.chdir(path)

# # # gather cells to create the PC space
df_coeffs_cells = []
names_cells_of_interest = ['noHGFstack1', 'noHGFstack2', 'noHGFstack3', 'noHGFstack4', 'noHGFstack5', 'noHGFstack6', 'noHGFstack7', 'noHGFstack8', 'noHGFstack9', 'noHGFstack10', 'postEMTcell1' , 'postEMTcell2', 'postEMTcell3' ,  'postEMTcell4' , 'postEMTcell5', 'postEMTcell6', 'postEMTcell7' ,  'postEMTcell8' , 'postEMTcell9', 'postEMTcell10']
#names_cells_of_interest = ['noHGFstack1', 'noHGFstack2']

# # gathers shape coefficients 
# for cellname in names_cells_of_interest:
#     df_coeffs_PCspace, df_withoutTp_PCspace, rotated_v_matrix_PCspace = get_shape_df_from_cellId(cellname, 'your.json', align=True)
#     df_coeffs_cells.append(df_withoutTp_PCspace)
# df_withoutTp_PCspace = pd.concat(df_coeffs_cells)
# print(df_withoutTp_PCspace)

# gathers shape independant descriptions
df_coeffs_cellsShapeInd = []
for cellname in names_cells_of_interest:
    df_coeffs_PCspaceShapeInd, df_withoutTp_PCspaceShapeInd, rotated_v_matrix_PCspaceShapeInd = get_shapeInd_df_from_cellId(cellname, 'your.json', align=True)
    df_coeffs_cellsShapeInd.append(df_withoutTp_PCspaceShapeInd)
df_withoutTp_PCspaceShapeInd = pd.concat(df_coeffs_cellsShapeInd)
print(df_withoutTp_PCspaceShapeInd)


# # center the dataframe used to fit the PC space
# mean, centered_df = center_df_data(df_withoutTp_PCspace)
meanShapeInd, centered_dfShapeInd = center_df_data(df_withoutTp_PCspaceShapeInd)

# # create the PCs space
# pca_fitted = create_PC_space(n_components=2, df_centered_data = centered_df)
pca_fittedShapeInd = create_PC_space(n_components=2, df_centered_data = centered_dfShapeInd)


# # #plot different shapes 

# # #plot the mean shape of the whole PC space 
# print('mean shape of PC space')
# meanPCA, reconstructed_mesh = plot_shape_from_PCspace_coords(pca_fitted, [0,0], 25, mean)

# # # #plot the mean shape of a subgroup of cell sin the PC space: epithelial cells 
# sub_df_witouthTp_PCspace = df_withoutTp_PCspace.iloc[:10]
# # print(sub_df_witouthTp_PCspace)
# print('mean epithelial shape')
# meanEpithelialShape, meanEpithelialMesh = plot_mean_shape_from_df_shcoeffs(sub_df_witouthTp_PCspace, lmaximum = 25)



# # # #plot the mean shape of a subgroup of cell sin the PC space: mesenchymal cells 
# sub_df_witouthTp_PCspace = df_withoutTp_PCspace.iloc[10:]
# # print(sub_df_witouthTp_PCspace)
# print('mean mesenchymal shape')
# meanMesenchymalShape, meanMesenchymalMesh =plot_mean_shape_from_df_shcoeffs(sub_df_witouthTp_PCspace, lmaximum = 25)

# # #plot the shapes corresponding to different coordinates in the PC space

# # meanPCA, recomesh0 = plot_shape_from_PCspace_coords(pca_fitted, [2,0.5], 25, mean)
# # meanPCA, recomesh1 = plot_shape_from_PCspace_coords(pca_fitted, [0,0.5], 25, mean) 
# # meanPCA, recomesh2 = plot_shape_from_PCspace_coords(pca_fitted, [-2,0.5], 25, mean)   
# # meanPCA, recomesh3 = plot_shape_from_PCspace_coords(pca_fitted, [2,0], 25, mean) 
# # meanPCA, recomesh4 = plot_shape_from_PCspace_coords(pca_fitted, [0,0], 25, mean) 
# # meanPCA, recomesh5 = plot_shape_from_PCspace_coords(pca_fitted, [-2,0], 25, mean) 
# # meanPCA, recomesh6 = plot_shape_from_PCspace_coords(pca_fitted, [2,-0.5], 25, mean) 
# # meanPCA, recomesh7 = plot_shape_from_PCspace_coords(pca_fitted, [0,-0.5], 25, mean) 
# # meanPCA, recomesh8 = plot_shape_from_PCspace_coords(pca_fitted, [-2,-0.5], 25, mean) 

# # recomesh = [recomesh0, recomesh1, recomesh2, recomesh3, recomesh4, recomesh5, recomesh6, recomesh7, recomesh8]

# # plot = Plotter(shape=(3,3))
# # plot.at(0).show(recomesh[0], __doc__, viewup='z')
# # plot.at(0).add(Text2D('PC1 = -2, PC2 = 0.5', pos='top-left', font="Arial", s=1))
# # plot.at(1).show(recomesh[1], __doc__, viewup='z')
# # plot.at(1).add(Text2D('PC1 = 0, PC2 = 0.5', pos='top-left', font="Arial", s=1))
# # plot.at(2).show(recomesh[2], __doc__, viewup='z')
# # plot.at(2).add(Text2D('PC1 = 2, PC2 = 0.5' , pos='top-left', font="Arial", s=1))
# # plot.at(3).show(recomesh[3], __doc__, viewup='z')
# # plot.at(3).add(Text2D('PC1 = -2, PC2 = 0', pos='top-left', font="Arial", s=1))
# # plot.at(4).show(recomesh[4], __doc__, viewup='z')
# # plot.at(4).add(Text2D('PC1 = 0, PC2 = 0', pos='top-left', font="Arial", s=1))
# # plot.at(5).show(recomesh[5], __doc__, viewup='z')
# # plot.at(5).add(Text2D('PC1 = 2, PC2 = 0', pos='top-left', font="Arial", s=1))
# # plot.at(6).show(recomesh[6], __doc__, viewup='z')
# # plot.at(6).add(Text2D('PC1 = 2, PC2 = -0.5', pos='top-left', font="Arial", s=1))
# # plot.at(7).show(recomesh[7], __doc__, viewup='z')
# # plot.at(7).add(Text2D('PC1 = 0, PC2 = -0.5', pos='top-left', font="Arial", s=1))
# # plot.at(8).show(recomesh[8], __doc__, viewup='z')
# # plot.at(8).add(Text2D('PC1 = 2, PC2 = -0.5', pos='top-left', font="Arial", s=1))
# # plot.interactive().close()

# # # # Project shape ind description on a shape 

# # #project mean distribution on mean shape of the PC space 
# reconstructed_mesh = Mesh(reconstructed_mesh)
# v_matrix_tp = np.array(reconstructed_mesh.vertices)
# faces = np.array(reconstructed_mesh.cells)
# binary_array_of_shape = get_binary_array_from_vertices(v_matrix_tp, faces)

# save_distribution_fromPC_coordinates(pca_fittedShapeInd, coords=[0,0], meanShapeInd, binary_array_of_shape, "mean distribution on mean shape")
# save_distribution_fromPC_coordinates(pca_fittedShapeInd, coords=[-1,0], meanShapeInd, v_matrix_tp, faces, "-PC1 on mean shape")
# save_distribution_fromPC_coordinates(pca_fittedShapeInd, coords=[1,0], meanShapeInd, v_matrix_tp, faces, "+PC1 on mean shape")
# save_distribution_fromPC_coordinates(pca_fittedShapeInd, coords=[0,-1], meanShapeInd, v_matrix_tp, faces, "-PC2 on mean shape")
# save_distribution_fromPC_coordinates(pca_fittedShapeInd, coords=[0,1], meanShapeInd, v_matrix_tp, faces, "+PC2 on mean shape")


# # #project mean epithelial distribution on mean shape of PC space 
# print('project mean epithelial distribution on mean shape of PC space ')
# meanEpiDistrib_raw = np.array(df_withoutTp_PCspaceShapeInd.iloc[:10].mean())
# meanDistribPCA = meanEpiDistrib_raw.reshape(1001,8194)
# meanDistrib = np.array(meanDistribPCA)

# gfp_morphed = project_distribution_on_shape(binary_array_of_shape, meanDistrib)

# save_morphed_signa_tiffs('mean epithelial distribution on mean shape', binary_array_of_shape, gfp_morphed)

# # #project mean epithelial distribution on mean epithelial shape 
# print('project mean epithelial distribution on mean epithelial shape ')

# meanEpithelialMesh = Mesh(meanEpithelialMesh)
# v_matrix_tp = np.array(meanEpithelialMesh.vertices)
# faces = np.array(meanEpithelialMesh.cells)
# binary_array_of_shape = get_binary_array_from_vertices(v_matrix_tp, faces)
# gfp_morphed = project_distribution_on_shape(binary_array_of_shape, meanDistrib)
# save_morphed_signa_tiffs('mean epithelial distribution on mean epithelial shape', binary_array_of_shape, gfp_morphed)




# # #project mean mesenchymal distribution on mean shape of PC space 
# print('project mean mesenchymal distribution on mean shape of PC space')
# meanMesDistrib_raw = np.array(df_withoutTp_PCspaceShapeInd.iloc[10:].mean())
# meanMesDistrib = meanMesDistrib_raw.reshape(1001,8194)
# meanDistrib = np.array(meanMesDistrib)


# reconstructed_mesh = Mesh(reconstructed_mesh)
# v_matrix_tp = np.array(reconstructed_mesh.vertices)
# faces = np.array(reconstructed_mesh.cells)
# binary_array_of_shape = get_binary_array_from_vertices(v_matrix_tp, faces)
# gfp_morphed = project_distribution_on_shape(binary_array_of_shape, meanDistrib)
# save_morphed_signa_tiffs('mean mesenchymal distribution on mean shape', binary_array_of_shape, gfp_morphed)



# # #project mean mesnechymal distribution on mean mesnechymal shape
# print('project mean mesnechymal distribution on mean mesnechymal shape')
# meanMesMesh = Mesh(meanMesenchymalMesh)
# v_matrix_tp = np.array(meanMesMesh.vertices)
# faces = np.array(meanMesMesh.cells)
# binary_array_of_shape = get_binary_array_from_vertices(v_matrix_tp, faces)
# gfp_morphed = project_distribution_on_shape(binary_array_of_shape, meanDistrib)

# save_morphed_signa_tiffs('mean mesenchymal distribution on mean mesenchymal shape', binary_array_of_shape, gfp_morphed)


# # #  visualise_cell('noHGFstack2', 'your.json', tp=0)

# # gfp_representation = np.array(gfp_representations[0])
# # shape_ind_description = gfp_representation
# # plt.imshow(gfp_representation, aspect='auto')
# # plt.show()

# # get mesh of shape we want to project gfp_morphed on
# # v_matrix_tp = rotated_v_matrix[0]
# # new_cell = read_one_cell_in_json('noHGFstack2', 'your.json')
# # faces = np.array(new_cell._segmentationFaces[0])

# # binary_array_of_shape = get_binary_array_from_vertices(v_matrix_tp, faces)

# # project_distribution_on_shape(binary_array_of_shape, shape_ind_description)




# # # # # # print the variance ratios of the PCs 
# # print('variance ratios of PCs shape')
# # variance_ratio = pca_fitted.explained_variance_ratio_
# # print(variance_ratio)

# # print('variance ratios of PCs distribution')
# # variance_ratioShapeInd = pca_fittedShapeInd.explained_variance_ratio_
# # print(variance_ratioShapeInd)


# # # # # # # # gather cells to project 
# # # # # df_coeffs_toProject, df_withoutTp_toProject, rotated_v_matrix_toProject = get_shape_df_from_cellId('20240529-sample1-cell2', 'your.json')


# # # #project data we want to visualise on the PC space 
# # #df_projected_data = project_on_PC_space(pca_fitted, df_centered_data = df_withoutTp_toProject - mean, columns_names=['PC1', 'PC2'])
# # df_projected_data = project_on_PC_space(pca_fitted, df_centered_data = centered_df, columns_names=['PC1', 'PC2'])
# # df_projected_data['PC1'] = df_projected_data['PC1'].apply(lambda x: x*-1)
# # #df_projected_data.insert(0,"Timepoint", list(df_coeffs_toProject.Timepoint) )
# df_projected_data.insert(0,"Timepoint", [0,0,0,0,0,0,0,0,0,0,21, 21, 21, 21, 21, 21, 21, 21, 21, 21] )
# #df_projected_data.insert(0,"Timepoint", [0,100] )

df_projected_dataShapeInd = project_on_PC_space(pca_fittedShapeInd, df_centered_data = centered_dfShapeInd, columns_names=['PC1', 'PC2'])

df_projected_dataShapeInd.insert(0,"Timepoint", [0,0,0,0,0,0,0,0,0,0,21, 21, 21, 21, 21, 21, 21, 21, 21, 21] )
# #df_projected_data.insert(0,"Timepoint", [0,100] )

# # # plot projected data in projected space
# #plot_single_trajectory_on_PC_space(df_projected_data)
# # plot_single_points_on_PC_space(df_projected_data)
plot_single_points_on_PC_space(df_projected_dataShapeInd)
 
# timepoints = [float(i) for i in df_projected_data.Timepoint.tolist()]
# fig = plt.figure(figsize=(8,5))
# ax = fig.subplots(1,1)
# sc = ax.scatter(df_projected_data.PC1, df_projected_data.PC2, s=50, c=timepoints, cmap = 'viridis', zorder=2) # zorder is order of layering points
# #ax.plot(df_projected_data.PC1, df_projected_data.PC2, 'k', zorder=1)
# plt.xlabel('PC1 shape - 63%',fontsize=25)
# plt.ylabel('PC2 shape - 12%',fontsize=25)
# plt.xticks(fontsize=25)
# plt.yticks(fontsize=25)
# cbar = plt.colorbar(sc)
# cbar.ax.tick_params(labelsize=20)
# #cbar.ax.set_ylabel("Timepoint", rotation=270,fontsize=25)
# plt.show()

# timepoints = [float(i) for i in df_projected_dataShapeInd.Timepoint.tolist()]
# fig = plt.figure(figsize=(8,5))
# ax = fig.subplots(1,1)
# sc = ax.scatter(df_projected_dataShapeInd.PC1, df_projected_dataShapeInd.PC2, s=50, c=timepoints, cmap = 'viridis', zorder=2) # zorder is order of layering points
# #ax.plot(df_projected_data.PC1, df_projected_data.PC2, 'k', zorder=1)
# plt.xlabel('PC1 distribution - 32%',fontsize=25)
# plt.ylabel('PC2 distribution - 15%',fontsize=25)
# plt.xticks(fontsize=25)
# plt.yticks(fontsize=25)
# cbar = plt.colorbar(sc)
# cbar.ax.tick_params(labelsize=20)
# #cbar.ax.set_ylabel("Timepoint", rotation=270,fontsize=25)
# plt.show()

# timepoints = [float(i) for i in df_projected_data.Timepoint.tolist()]
# fig = plt.figure(figsize=(8,5))
# ax = fig.subplots(1,1)
# sc = ax.scatter(df_projected_data.PC1, df_projected_dataShapeInd.PC1, s=50, c=timepoints, cmap = 'viridis', zorder=2) # zorder is order of layering points
# #ax.plot(df_projected_data.PC1, df_projected_data.PC2, 'k', zorder=1)
# plt.xlabel('PC1 shape',fontsize=25)
# plt.ylabel('PC1 distribution',fontsize=25)
# plt.xticks(fontsize=25)
# plt.yticks(fontsize=25)
# cbar = plt.colorbar(sc)
# cbar.ax.tick_params(labelsize=20)
# #cbar.ax.set_ylabel("Timepoint", rotation=270,fontsize=25)
# plt.show()
