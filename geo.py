import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def getModel(resolution=50):
  x = np.linspace(-10, 10, resolution)  # X coordinates from -10 to 10
  y = np.linspace(-10, 10, resolution)  # Y coordinates from -10 to 10
  z = np.linspace(-10, 10, resolution)  # Z coordinates from -10 to 10
  X, Y, Z = np.meshgrid(x, y, z, indexing='ij')  # Create a 3D meshgrid for X, Y, and Z coordinates

  # Combine flattened arrays into a 2D numpy array where each row is an (x, y, z) coordinate
  xyz = np.column_stack((X.flatten(), Y.flatten(), Z.flatten()))
  return xyz, X, Y, Z

def rotate(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / np.sqrt(np.dot(axis, axis))
    a = np.cos(theta / 2.0)
    b, c, d = -axis * np.sin(theta / 2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

class Layer:
    def __init__(self, base, width, value):
        self.base = base
        self.width = width
        self.value = value

    def run(self, xyz, data):
        # Extract z-coordinates
        z_coords = xyz[:, 2]

        # Create a mask where z is within the specified range and data is None
        mask = (self.base <= z_coords) & (z_coords <= self.base + self.width) & (np.isnan(data))

        # Apply the mask and update data where condition is met
        data[mask] = self.value

        # Return the unchanged xyz and the potentially modified data
        return xyz, data

class Tilt:
    def __init__(self, strike, dip):
        self.strike = np.radians(strike)  # Convert degrees to radians
        self.dip = np.radians(dip)  # Convert degrees to radians

    def run(self, xyz, data):
        # Calculate rotation axis from strike (rotation around z-axis)
        axis = rotate([0, 0, 1], -self.strike) @ [0, 1, 0]

        # Calculate rotation matrix from dip (tilt)
        R = rotate(axis, -self.dip)

        # Apply rotation to xyz points
        return xyz @ R.T, data  # Assuming xyz is an Nx3 numpy array

class Fold:
    def __init__(self, strike = 0, dip = 90, rake = 0, period = 50, amplitude = 10, shape = 0, offset=0, point=[0, 0, 0]):
        self.strike = np.radians(strike)
        self.dip = np.radians(dip)
        self.rake = np.radians(rake)
        self.period = period
        self.amplitude = amplitude
        self.shape = shape
        self.offset = offset
        self.point = np.array(point)

    def run(self, xyz, data):
        # Calculate rotation matrices
        M1 = rotate([0, 0, 1], -(self.rake + np.pi / 2))
        M2 = rotate([0, 1, 0], -(self.dip))
        M3 = rotate([0, 0, 1], -(self.strike))

        # Normal vector of the fold
        N = M3 @ M2 @ M1 @ [0., 1.0, 0.0]
        M1 = rotate([0, 0, 1], -(self.rake))
        V = M3 @ M2 @ M1 @ [0., 1.0, 0.0]
        U = np.cross(N[:3], V[:3])

        new_xyz = np.empty_like(xyz)
        for i, point in enumerate(xyz):
            v0 = point - self.point
            fU = np.dot(v0, U) / np.linalg.norm(U) - self.offset * self.period
            inside = 2 * np.pi * fU / self.period
            off = self.amplitude * (np.cos(inside) + self.shape * np.cos(3 * inside))
            T = N * off
            new_xyz[i] = point + T[:3]

        return new_xyz, data

class Dike:
    def __init__(self, strike, dip, width, point, data_value):
        self.strike = np.radians(strike)
        self.dip = np.radians(dip)
        self.width = width
        self.point = np.array(point)
        self.data_value = data_value

    def run(self, xyz, data):
        # Calculate rotation matrices
        M1 = rotate([0, 0, 1], -self.strike)  # Rotation around z-axis for strike
        M2 = rotate([0, 1, 0], self.dip)      # Rotation around y-axis for dip

        # Normal vector of the dike plane
        N = M1 @ M2 @ [0.0, 0.0, 1.0]
        N /= np.linalg.norm(N)  # Normalize the normal vector

        # Calculate distances from the dike plane
        dists = np.dot(xyz - self.point, N)

        # Update data based on whether points are within the width of the dike
        data[np.abs(dists) <= self.width / 2.0] = self.data_value

        return xyz, data


def plotCrossSection(data, X, Y, Z):
    # Reshape the data back to the 3D grid
    data_reshaped = data.reshape(X.shape)
    # # Plotting a cross-section, for example at z index 10
    y_index = 10  # Change this index to view different vertical cross-sectional slices
    plt.figure(figsize=(10, 8))
    plt.pcolormesh(X[:, y_index, :], Z[:, y_index, :], data_reshaped[:, y_index, :], shading='auto')
    plt.colorbar()  # Show color scale
    plt.title(f'Vertical Cross-section at Y index {y_index}')
    plt.xlabel('X coordinate')
    plt.ylabel('Z coordinate')
    plt.show()

def plot3D(data, X, Y, Z):
    # Reshape the data back to the 3D grid
    data_reshaped = data.reshape(X.shape)
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Choose indices for the cross sections
    y_index = 10  # for the X-Z plane cross-section
    x_index = 10  # for the Y-Z plane cross-section

    # Plot the X-Z cross-section as a texture on a Y constant plane
    XZ_x, XZ_z = X[:, y_index, :], Z[:, y_index, :]
    XZ_slice = data_reshaped[:, y_index, :]
    ax.plot_surface(XZ_x, np.full_like(XZ_x, Y[0, y_index, 0]), XZ_z, rstride=1, cstride=1, facecolors=plt.cm.viridis(XZ_slice / np.nanmax(XZ_slice)))

    # Plot the Y-Z cross-section as a texture on an X constant plane
    YZ_y, YZ_z = Y[x_index, :, :], Z[x_index, :, :]
    YZ_slice = data_reshaped[x_index, :, :]
    ax.plot_surface(np.full_like(YZ_y, X[x_index, 0, 0]), YZ_y, YZ_z, rstride=1, cstride=1, facecolors=plt.cm.viridis(YZ_slice / np.nanmax(YZ_slice)))

    # Set labels
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_zlabel('Z Coordinate')

    # Color bar
    mappable = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(vmin=data_reshaped.min(), vmax=data_reshaped.max()))
    mappable.set_array([])
    cbar = plt.colorbar(mappable, ax=ax, orientation='vertical')
    cbar.set_label('Data value')

    # Set the aspect ratio of the plot
    ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio

    plt.title('3D Plot with X-Z and Y-Z Cross Section Images')
    plt.show()


def ModelHistory(xyz, history):
    # Clone the xyz array to avoid modifying the original input
    cloned_xyz = np.array(xyz, copy=True)
    # Initialize data for each coordinate point as NaNs, in a numpy array
    data = np.full(cloned_xyz.shape[0], np.nan)

    # Apply each transformation in the history in reverse order
    for transformation in reversed(history):
        cloned_xyz, data = transformation.run(cloned_xyz, data)

    return data
