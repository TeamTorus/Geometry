# Import necessary libraries
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import pyvista as pv
from scipy.integrate import quad
from scipy.interpolate import interp1d
from nurbs_gen import nurbs_gen

# symbolic variables
t, s = sp.symbols('t s', real=True)

# Domain and resolution
s_domain = [0, 1]  # Domain for the curve parameter s
t_domain = [0, 2]  # Domain for the shape parameter t
s_resolution = 100  # Resolution for discretizing the curve parameter s
t_resolution = 100  # Resolution for discretizing the shape parameter t
hub_resolution = 50 # Resolution for discretizing the hub (needs to equal s_resolution for matplotlib viz)

normalize_blade_mesh = False        # Normalize the blade mesh to have uniform arc length wrt t
apply_thickness_normal = False      # Apply airfoil thickness normal to camber line
close_cylinder = True               # Close the cylinder mesh for the hub with top and bottom faces

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

# ---------------------------------------- PT 1: Define the 2D shape ----------------------------------------
# https://en.wikipedia.org/wiki/NACA_airfoil

yt = 5 * thickness * (0.2969 * sp.sqrt(t) - 0.1260 * t - 0.3516 * t**2 + 0.2843 * t**3 - 0.1036 * t**4)
yc = sp.Piecewise(((m / p**2) * (2 * p * t - t**2), t <= p), ((m / (1 - p)**2) * ((1 - 2 * p) + 2 * p * t - t**2), t > p))

if apply_thickness_normal:
    # define camber normals for proper thickness distribution
    dyc_dt = sp.Piecewise((2*m/p**2 * (p - t), t <= p), (2*m/(1 - p)**2 * (p - t), t > p))
    theta_c = sp.atan(dyc_dt)
else:
    theta_c = 0

# apply camber + thickness
yu = yc + yt * sp.cos(theta_c)
yl_1 = yc - yt * sp.cos(theta_c)

xu = t - yt * sp.sin(theta_c)
xl = (2 - t) + yt * sp.sin(theta_c)             # 2-t for the t domain shift

# domain shift for two halves of the airfoil 
yl = yl_1.subs(t, (2 - t))

# use activations for upper/lower
y_2D = yu * sp.Heaviside(1 - t) + yl * sp.Heaviside(t - 1)
x_2D = xu * sp.Heaviside(1 - t) + xl * sp.Heaviside(t - 1)

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

# -------------------------------------- PT 2: Normalize parameter domains  --------------------------------------

if normalize_blade_mesh:
    # get gradient components
    dyu_dt = sp.diff(yu, t)
    dyl_dt = sp.diff(yl, t)

    dxu_dt = sp.diff(xu, t)
    dxl_dt = sp.diff(xl, t)

    # arc length of the airfoil
    ds_u = sp.sqrt(dxu_dt**2 + dyu_dt**2)
    ds_l = sp.sqrt(dxl_dt**2 + dyl_dt**2)

    # NACA airfoil isn't analytically integrable, so use numerical integration
    arc_length_upper_func = sp.lambdify(t, ds_u, 'numpy')
    arc_length_lower_func = sp.lambdify(t, ds_l, 'numpy')

    # define fine-grained t values for numerical integration
    t_fine_upper = np.linspace(t_domain[0], t_domain[1]/2, 500)
    t_fine_lower = np.linspace(t_domain[1]/2, t_domain[1], 500)

    # integrate for the upper curve
    s_upper_vals = [quad(arc_length_upper_func, 0, t_)[0] for t_ in t_fine_upper]

    # start from s_upper_vals[-1] to maintain continuity
    s_lower_vals = [quad(arc_length_lower_func, 1, t_)[0] + s_upper_vals[-1] for t_ in t_fine_lower]

    t_fine = np.concatenate((t_fine_upper, t_fine_lower))
    s_fine = np.concatenate((s_upper_vals, s_lower_vals))

    # Create an interpolation function to estimate t given s
    t_of_s = interp1d(s_fine, t_fine, kind='linear', fill_value='extrapolate')

    # Now t_of_s(s) can be used to get t for any s
    t_vals = t_of_s(np.linspace(0, s_fine[-1], t_resolution))
    t_vals = np.clip(t_vals, t_domain[0], t_domain[1])
else:
    t_vals = np.linspace(t_domain[0], t_domain[1], t_resolution)

s_vals = np.linspace(s_domain[0], s_domain[1], s_resolution)



# ------------------------------------------ PT 3: Create the 3D curve  ------------------------------------------

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

# -------------------------------------- PT 4: Frenet-Serret frame along curve  --------------------------------------

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

# Transform the 2D shape Along the Curve

# Define the parametric position along the curve
C = sp.Matrix([x_curve, y_curve, z_curve])

# Express the final 3D coordinates of the airfoil after extrusion
X_final = C[0] + X_rotated_scaled * N[0] + Y_rotated_scaled * B[0]
Y_final = C[1] + X_rotated_scaled * N[1] + Y_rotated_scaled * B[1]
Z_final = C[2] + X_rotated_scaled * N[2] + Y_rotated_scaled * B[2]

print("X_final:")
print(X_final)

# -------------------------------------- PT 5: Assemble the 3D Propeller  --------------------------------------

def generate_cylinder_mesh(radius, length, resolution):
    # Create a grid of points for the cylinder
    theta_vals = np.linspace(0, 2 * np.pi, resolution)
    z_vals = np.linspace(-length / 2, length / 2, resolution)

    # Create the meshgrid
    theta_mesh, z_mesh = np.meshgrid(theta_vals, z_vals)

    # Convert to Cartesian coordinates
    X = radius * np.cos(theta_mesh)
    Y = radius * np.sin(theta_mesh)
    Z = z_mesh

    return X, Y, Z

# Generate the hub mesh
X_hub, Y_hub, Z_hub = generate_cylinder_mesh(hub_radius, hub_length, hub_resolution)

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
    s_mesh, t_mesh = np.meshgrid(s_vals, t_vals)

    # evaluate it on the meshgrid
    X_vals = X_func(s_mesh, t_mesh)
    Y_vals = Y_func(s_mesh, t_mesh)
    Z_vals = Z_func(s_mesh, t_mesh)

    # stack them horizontally
    X_prop[:, i * s_resolution:(i + 1) * s_resolution] = X_vals
    Y_prop[:, i * s_resolution:(i + 1) * s_resolution] = Y_vals
    Z_prop[:, i * s_resolution:(i + 1) * s_resolution] = Z_vals


# -------------------------------- PT 6: Visualize the propeller in matplotlib  --------------------------------
# # Combine the hub and blades into a single array
# X_tot = np.hstack((X_prop, X_hub))
# Y_tot = np.hstack((Y_prop, Y_hub))
# Z_tot = np.hstack((Z_prop, Z_hub))

# # Visualize the assembled propeller in 3D
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.plot_surface(X_tot, Y_tot, Z_tot, cmap='viridis', edgecolor='none')

# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')
# ax.set_title('Extruded Toroidal Propeller with Multiple Blades')
# plt.show()

# --------------------------------------- PT 7: Boolean it with PyVista  ---------------------------------------
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

def close_cylinder_mesh(X, Y, Z):
    # Flatten the grid arrays
    points = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))

    nx, ny = X.shape
    faces = []

    for i in range(nx - 1):
        for j in range(ny):
            idx0 = i * ny + j
            idx1 = i * ny + (j + 1) % ny
            idx2 = (i + 1) * ny + (j + 1) % ny
            idx3 = (i + 1) * ny + j
            faces.extend([4, idx0, idx1, idx2, idx3])

    # center points for top and bottom faces
    top_center = np.array([[0, 0, Z.max()]])
    bottom_center = np.array([[0, 0, Z.min()]])
    points = np.vstack((points, top_center, bottom_center))
    idx_top_center = len(points) - 2
    idx_bottom_center = len(points) - 1

    # Indices of top and bottom circle points
    idx_top = (nx - 1) * ny + np.arange(ny)
    idx_bottom = 0 * ny + np.arange(ny)

    # Create faces for the top cap
    for j in range(ny):
        idx0 = idx_top[j]
        idx1 = idx_top[(j + 1) % ny]
        faces.extend([3, idx0, idx1, idx_top_center])

    # Create faces for the bottom cap
    for j in range(ny):
        idx0 = idx_bottom[j]
        idx1 = idx_bottom[(j + 1) % ny]
        faces.extend([3, idx_bottom_center, idx1, idx0])

    faces = np.array(faces)

    mesh = pv.PolyData(points)
    mesh.faces = faces

    mesh = mesh.triangulate()
    return mesh

mesh_blades = create_mesh_from_grids(X_prop, Y_prop, Z_prop)
mesh_hub = close_cylinder_mesh(X_hub, Y_hub, Z_hub)

print("Blades mesh is manifold:", mesh_blades.is_manifold)
print("Hub mesh is manifold:", mesh_hub.is_manifold)

mesh_hub_m = mesh_hub.fill_holes(100).clean()
mesh_blades_m = mesh_blades.fill_holes(100).clean()


final_mesh = mesh_hub_m.boolean_union(mesh_blades, tolerance=1e-4)
final_mesh = final_mesh.fill_holes(3, inplace=True)  # fill small mesh holes idk y they here


final_mesh.plot_normals(mag=0.25, show_edges=True)
# final_mesh.plot(show_edges=True)

# export the final mesh to STL
final_mesh.save('toroidal_propeller.stl')


