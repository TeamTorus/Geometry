# Import necessary libraries
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pyvista as pv
from nurbs_gen import nurbs_gen


# symbolic variables
t, s = sp.symbols('t s', real=True)

# Domain and resolution
s_domain = [0, 1]  # Domain for the curve parameter s
t_domain = [0, 2]  # Domain for the shape parameter t
s_resolution = 100  # Resolution for discretizing the curve parameter s
t_resolution = 100  # Resolution for discretizing the shape parameter t

# Modifiable parameters
hub_radius = 5  # Radius of the cylindrical hub
hub_length = 20  # Length of the cylindrical hub
num_blades = 3  # Number of blades

# Airfoil Params
m = 0.02
p = 0.4
thickness = 0.5

# Centerline Params
loc_ctrl_point2 = [3, 3, 10]
loc_ctrl_point3 = [7, 6, 8]
blade_vector = [8, 8]   # offset between the two endpoints

# Angle of Attack
a_AoA = 0
b_AoA = 0
c_AoA = 0
d_AoA = 2 * np.pi
e_AoA = 0

# Scaling Params
a_scX = 0
b_scX = 0
c_scX = 0
d_scX = 1
e_scX = 2

a_scY = 0
b_scY = 0
c_scY = 0
d_scY = 1
e_scY = 2

# define the 2D Shape
# https://en.wikipedia.org/wiki/NACA_airfoil
x_2D = t * sp.Heaviside(1 - t) + (2 - t) * sp.Heaviside(t - 1)

yt = 5 * thickness * (0.2969 * sp.sqrt(t) - 0.1260 * t - 0.3516 * t**2 + 0.2843 * t**3 - 0.1036 * t**4)
yc = (m / p**2) * (2 * p * t - t**2) * sp.Heaviside(p - t) + (m / (1 - p)**2) * ((1 - 2 * p) + 2 * p * t - t**2) * sp.Heaviside(t - p)
yu = yt + yc

# domain shift
yl_1 = yc - yt
yl = yl_1.subs(t, (2 - t))

# use activations for upper/lower
y_2D = yu * sp.Heaviside(1 - t) + yl * sp.Heaviside(t - 1)

# define the 3D Curve

blade_hub_radius = hub_radius - hub_radius / 8

# Pick first point
ctrl_point1 = [blade_hub_radius, 0, hub_length / 2 - 1]

# For last control point, offset by blade_vector
# blade_vector[0] represents a movement along the circumference of the hub
# blade_vector[1] represents a movement along the z-axis
disp_theta = blade_vector[0] / (2 * np.pi * blade_hub_radius) * 2 * np.pi
ctrl_point4 = [blade_hub_radius * np.cos(disp_theta), blade_hub_radius * np.sin(disp_theta), -1 * blade_vector[1]]

# Define the control points for the curve (converting to cylindrical and translating)
disp_theta2 = loc_ctrl_point2[0] / (2 * np.pi * blade_hub_radius) * 2 * np.pi
loc2_radius = blade_hub_radius + loc_ctrl_point2[2]
ctrl_point2 = [loc2_radius * np.cos(disp_theta2), loc2_radius * np.sin(disp_theta2), -1 * loc_ctrl_point2[1]]

disp_theta3 = loc_ctrl_point3[0] / (2 * np.pi * blade_hub_radius) * 2 * np.pi
loc3_radius = blade_hub_radius + loc_ctrl_point3[2]
ctrl_point3 = [loc3_radius * np.cos(disp_theta3), loc3_radius * np.sin(disp_theta3), -1 * loc_ctrl_point3[1]]

control_points = np.array([ctrl_point1, ctrl_point2, ctrl_point3, ctrl_point4])  # [x, y, z]


# could add weights as a parameter, but for now just use uniform weights
weights = [1, 1, 1, 1]
x_curve, y_curve, z_curve = nurbs_gen(s, control_points, weights, to_plot=False)

# Define the Frenet-Serret Frame

# Derivatives to get tangent vector (T)
dx_ds = sp.diff(x_curve, s)
dy_ds = sp.diff(y_curve, s)
dz_ds = sp.diff(z_curve, s)
T = sp.Matrix([dx_ds, dy_ds, dz_ds])
T_norm = sp.sqrt(T.dot(T))
T = T / T_norm

# Derivative of tangent to get normal vector (N)
dT_ds = sp.diff(T, s)
N_norm = sp.sqrt(dT_ds.dot(dT_ds))
N = dT_ds / N_norm

# Binormal vector (B) is the cross product of T and N
B = T.cross(N)

# rotation angle function of s
AoA = a_AoA * s**4 + b_AoA * s**3 + c_AoA * s**2 + d_AoA * s + e_AoA

# counterclockwise in the local x-y plane
rotation_matrix = sp.Matrix([[sp.cos(AoA), -sp.sin(AoA)], [sp.sin(AoA), sp.cos(AoA)]])

# apply rotation
XY_rotated = rotation_matrix * sp.Matrix([x_2D, y_2D]) 
X_rotated = XY_rotated[0]
Y_rotated = XY_rotated[1]

# scaling/stretching functions (parametric by s) ~ approximates airfoil rotation transformations
scale_x = a_scX * s**4 + b_scX * s**3 + c_scX * s**2 + d_scX * s + e_scX  # Parametric scaling for x
scale_y = a_scY * s**4 + b_scY * s**3 + c_scY * s**2 + d_scY * s + e_scY  # Parametric scaling for y

# Apply scaling in the local frame (after rotation)
X_rotated_scaled = X_rotated * scale_x
Y_rotated_scaled = Y_rotated * scale_y

# Transform the 2D shape Along the Curve

# Define the parametric position along the curve
C = sp.Matrix([x_curve, y_curve, z_curve])

# Express the final 3D coordinates of the airfoil after extrusion
X_final = C[0] + X_rotated_scaled * N[0] + Y_rotated_scaled * B[0]
Y_final = C[1] + X_rotated_scaled * N[1] + Y_rotated_scaled * B[1]
Z_final = C[2] + X_rotated_scaled * N[2] + Y_rotated_scaled * B[2]

print("X_final:")
print(X_final)


# Generate Cylinder Hub (parametric)
theta_hub = np.linspace(0, 2 * np.pi, s_resolution)
z_hub = np.linspace(-hub_length / 2, hub_length / 2, t_resolution)
TH, ZH = np.meshgrid(theta_hub, z_hub)

X_hub = hub_radius * np.cos(TH)
Y_hub = hub_radius * np.sin(TH)
Z_hub = ZH

# convert to cylindrical coordinates
R = sp.sqrt(X_final**2 + Y_final**2)
Theta = sp.atan2(Y_final, X_final)

# create arrays for storing all blades
X_prop = np.zeros((t_resolution, s_resolution * num_blades))
Y_prop = np.zeros((t_resolution, s_resolution * num_blades))
Z_prop = np.zeros((t_resolution, s_resolution * num_blades))

# generate multiple blades around the hub
for i in range(num_blades):
    # rotate the blade by 2pi/num_blades
    Loc_Theta = Theta + i * 2 * np.pi / num_blades

    # Cartesian coordinates
    X_rotated = R * sp.cos(Loc_Theta)
    Y_rotated = R * sp.sin(Loc_Theta)

    # convert to evaluatable lambdas
    X_func = sp.lambdify((s, t), X_rotated, modules=['numpy', {'Heaviside': np.heaviside}])
    Y_func = sp.lambdify((s, t), Y_rotated, modules=['numpy', {'Heaviside': np.heaviside}])
    Z_func = sp.lambdify((s, t), Z_final, modules=['numpy', {'Heaviside': np.heaviside}])

    # discretize the parameters now, it makes it easier
    # TODO: get symbolic expression for all 3 blades, do a domain rescale or shift with heavisides to get the 3
    s_vals = np.linspace(s_domain[0], s_domain[1], s_resolution)
    t_vals = np.linspace(t_domain[0], t_domain[1], t_resolution)
    s_mesh, t_mesh = np.meshgrid(s_vals, t_vals)

    # evaluate it on the meshgrid
    X_vals = X_func(s_mesh, t_mesh)
    Y_vals = Y_func(s_mesh, t_mesh)
    Z_vals = Z_func(s_mesh, t_mesh)

    # stack them horizontally
    X_prop[:, i * s_resolution:(i + 1) * s_resolution] = X_vals
    Y_prop[:, i * s_resolution:(i + 1) * s_resolution] = Y_vals
    Z_prop[:, i * s_resolution:(i + 1) * s_resolution] = Z_vals

# Combine the hub and blades into a single array
X_tot = np.hstack((X_prop, X_hub))
Y_tot = np.hstack((Y_prop, Y_hub))
Z_tot = np.hstack((Z_prop, Z_hub))

# Visualize the assembled propeller in 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X_tot, Y_tot, Z_tot, cmap='viridis', edgecolor='none')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('Extruded Toroidal Propeller with Multiple Blades')
plt.show()

#------------------------ PT 2 ------------------------
def create_mesh_from_grids(X, Y, Z):
    # Flatten the grid arrays
    points = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))

    mesh = pv.PolyData(points)
    # generate faces from the grid
    nx, ny = X.shape
    faces = []

    for i in range(nx - 1):
        for j in range(ny - 1):
            idx0 = i * ny + j
            idx1 = idx0 + 1
            idx2 = idx0 + ny + 1
            idx3 = idx0 + ny
            faces.extend([4, idx0, idx1, idx2, idx3])

    # faces list to NumPy array
    faces = np.array(faces)
    mesh.faces = faces

    # triangulate grid mesh (it don't like quads)
    mesh = mesh.triangulate()
    return mesh

mesh_blades = create_mesh_from_grids(X_prop, Y_prop, Z_prop)
mesh_hub = create_mesh_from_grids(X_hub, Y_hub, Z_hub)

print("Blades mesh is manifold:", mesh_blades.is_manifold)
print("Hub mesh is manifold:", mesh_hub.is_manifold)

mesh_hub_m = mesh_hub.fill_holes(100).clean()
mesh_blades_m = mesh_blades.fill_holes(100).clean()

final_mesh = mesh_hub.boolean_union(mesh_blades, tolerance=1e-5).clean()
final_mesh = final_mesh.fill_holes(1, inplace=True)  # fill small mesh holes idk y they here


final_mesh.plot_normals(mag=0.25, show_edges=True)

# export the final mesh to STL
final_mesh.save('toroidal_propeller.stl')


