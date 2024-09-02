# Import required packages
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate as spinterp
import pyshtools
import plotly.express as px
#import plotly.figure_factory as ff
from vedo import Mesh, show

def shapecentroidtoorigin(v_matrix):
    
    # v_matrix = matrix with format [[x1 y1 z1] [x2 y2 z2]] specifying the vertices of a shape
    centroid = v_matrix.mean(axis=0, keepdims=True)
    # Translate to origin
    v_matrix_new = v_matrix - centroid

    return v_matrix_new


def getshcoefficients(v_matrix, lmax):

    # v_matrix = matrix with format [[x1 y1 z1] [x2 y2 z2]] specifying the vertices of a shape
    # lmax = maximum l coefficient

    x = v_matrix[:,0]
    y = v_matrix[:,1]
    z = v_matrix[:,2]

    # Cartesian to spherical coordinates convertion

    rad = np.sqrt(x**2 + y**2 + z**2)
    lat = np.arccos(np.divide(z, rad, out=np.zeros_like(rad), where=(rad != 0)))
    lon = np.pi + np.arctan2(y, x)

    # Put lon and lat coordinates into the same numpy array

    points = np.concatenate(
        [np.array(lon).reshape(-1, 1), np.array(lat).reshape(-1, 1)], axis=1)

    # Make a lon lat equally spaced grid

    grid_lon, grid_lat = np.meshgrid(
        np.linspace(start=0, stop=2 * np.pi, num=256, endpoint=True),
        np.linspace(start=0, stop=1 * np.pi, num=128, endpoint=True),)

    # Interpolate the (lon,lat,r) data into a grid
    grid = spinterp.griddata(points, rad, (grid_lon, grid_lat), method="nearest")

    # Get spherical harmonics coefficients from the grid

    coeffs = pyshtools.expand.SHExpandDH(grid, sampling=2, lmax_calc=lmax)
    
    # Make a dict of coefficients and which l,m they correspond to
    # This is to make a dataframe of the spherical harmonics coefficients for each shape
    
    lvalues = np.repeat(np.arange(lmax + 1).reshape(-1, 1), lmax + 1, axis=1)

    keys = []
    for suffix in ["C", "S"]:
        for (l, m) in zip(lvalues.flatten(), lvalues.T.flatten()):
            keys.append(f"shcoeffs_L{l}M{m}{suffix}")

    coeffs_dict = dict(zip(keys, coeffs.flatten()))

    return grid, coeffs, coeffs_dict

def reconstructshape(coeffs, npoints):

    # coeffs = spherical harmonics coefficients (format is output from pyshtools.expand.SHExpandDH function)
    # npoints = number of vertices wanted for the reconstructed shape

    # Reconstruct grid from coefficients

    grid_rec = pyshtools.expand.MakeGridDH(coeffs, sampling=2)

    # Reconstruct shape from grid

    rcentroid = (0,0,0)

    res_lat = grid_rec.shape[0]
    res_lon = grid_rec.shape[1]

    # Creates an interpolator
    rlon = np.linspace(start=0, stop=2 * np.pi, num=res_lon, endpoint=True)
    rlat = np.linspace(start=0, stop=1 * np.pi, num=res_lat, endpoint=True)

    # Input grid from coefficients (probably quite coarse grained) is interpolated onto a finer grid

    fgrid = spinterp.RectBivariateSpline(rlon, rlat, grid_rec.T)

    # Create x,y,z coordinates based on the Fibonacci Lattice

    golden_ratio = 0.5 * (1 + np.sqrt(5))
    idxs = np.arange(0, npoints, dtype=np.float32)
    fib_theta = np.arccos(2 * ((idxs + 0.5) / npoints) - 1)
    fib_phi = (2 * np.pi * (idxs / golden_ratio)) % (2 * np.pi) - np.pi

    fib_lat = fib_theta
    fib_lon = fib_phi + np.pi

    fib_grid = fgrid.ev(fib_lon, fib_lat) # Gives values of spline (finer grid) at equal intervals on a sphere

    # Assign to sphere

    fib_x = rcentroid[0] + fib_grid * np.sin(fib_theta) * np.cos(fib_phi)
    fib_y = rcentroid[1] + fib_grid * np.sin(fib_theta) * np.sin(fib_phi)
    fib_z = rcentroid[2] + fib_grid * np.cos(fib_theta)

    # Concatenate fib_x,y,z so that they can be input into the Delaunay triangulation function

    reconstructed_vertices = np.column_stack((fib_x,fib_y,fib_z))

    return grid_rec, reconstructed_vertices

def plotshapeverticesradius(v_matrix, title):

    # v_matrix = matrix with format [[x1 y1 z1] [x2 y2 z2]] specifying the vertices of a shape

    x = v_matrix[:,0]
    y = v_matrix[:,1]
    z = v_matrix[:,2]

    rad = np.sqrt(x**2 + y**2 + z**2)

    fig = px.scatter_3d(x=x, y=y, z=z, color=rad, title=title)
    fig.update_layout(coloraxis_colorbar_title_text = 'Radius')
    fig.show()

def plotgrid(grid, title):

    # grid = interpolated longitude latitude grid, dimensions 2n*n where each point has a radius value attached to it
    # title = plot title (string)

    plt.imshow(grid, origin='lower', interpolation='none',extent=[0,2*np.pi,0,np.pi])
    plt.title(title)
    plt.xlabel('Longitude (phi)')
    plt.ylabel('Latitude (theta)')
    cbar = plt.colorbar(shrink=0.5)
    cbar.ax.set_title('    Radius')
    plt.show()