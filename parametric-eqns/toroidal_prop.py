# Import necessary libraries
import numpy as np
import sympy as sp
import pyvista as pv
from scipy.spatial import Delaunay

# Define symbolic variables
t, s = sp.symbols('t s', real=True)

# Domain and resolution
s_domain = [0, 1]
t_domain = [0, 2]
s_resolution = 50
t_resolution = 50

# Modifiable parameters
hub_radius = 5
hub_length = 20
num_blades = 3

# Airfoil Params
m = 0.02
p = 0.4
thickness = 0.5

# Centerline Params
loc_ctrl_point2 = [3, 3, 10]
loc_ctrl_point3 = [7, 6, 8]
blade_vector = [8, 8]

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

# Define the 2D Airfoil Shape
x_2D = t * sp.Heaviside(1 - t) + (2 - t) * sp.Heaviside(t - 1)

yt = 5 * thickness * (
    0.2969 * sp.sqrt(t) - 0.1260 * t - 0.3516 * t**2 +
    0.2843 * t**3 - 0.1036 * t**4
)
yc = (
    (m / p**2) * (2 * p * t - t**2) * sp.Heaviside(p - t) +
    (m / (1 - p)**2) * ((1 - 2 * p) + 2 * p * t - t**2) * sp.Heaviside(t - p)
)
yu = yt + yc

# Domain shift for lower surface
yl_1 = yc - yt
yl = yl_1.subs(t, (2 - t))

# Combine upper and lower surfaces
y_2D = yu * sp.Heaviside(1 - t) + yl * sp.Heaviside(t - 1)

# Define the 3D Curve (Blade Centerline)
blade_hub_radius = hub_radius - hub_radius / 8

# First control point
ctrl_point1 = [blade_hub_radius, 0, hub_length / 2 - 1]

# Last control point
disp_theta = blade_vector[0] / (2 * np.pi * blade_hub_radius) * 2 * np.pi
ctrl_point4 = [
    blade_hub_radius * np.cos(disp_theta),
    blade_hub_radius * np.sin(disp_theta),
    -blade_vector[1]
]

# Intermediate control points
disp_theta2 = loc_ctrl_point2[0] / (2 * np.pi * blade_hub_radius) * 2 * np.pi
loc2_radius = blade_hub_radius + loc_ctrl_point2[2]
ctrl_point2 = [
    loc2_radius * np.cos(disp_theta2),
    loc2_radius * np.sin(disp_theta2),
    -loc_ctrl_point2[1]
]

disp_theta3 = loc_ctrl_point3[0] / (2 * np.pi * blade_hub_radius) * 2 * np.pi
loc3_radius = blade_hub_radius + loc_ctrl_point3[2]
ctrl_point3 = [
    loc3_radius * np.cos(disp_theta3),
    loc3_radius * np.sin(disp_theta3),
    -loc_ctrl_point3[1]
]

control_points = np.array([
    ctrl_point1, ctrl_point2, ctrl_point3, ctrl_point4
])

# NURBS generation function
def nurbs_gen(s, control_points, weights):
    degree = 3
    knot_vector = np.array([0, 0, 0, 0, 1, 1, 1, 1])
    n = len(control_points) - 1

    def nurbs_basis(i, p, u, knot_vector):
        if p == 0:
            return sp.Piecewise(
                (1, (u >= knot_vector[i]) & (u < knot_vector[i+1])),
                (0, True)
            )
        else:
            denom1 = knot_vector[i+p] - knot_vector[i]
            denom2 = knot_vector[i+p+1] - knot_vector[i+1]
            term1 = 0
            term2 = 0
            if denom1 != 0:
                term1 = ((u - knot_vector[i]) / denom1) * \
                        nurbs_basis(i, p-1, u, knot_vector)
            if denom2 != 0:
                term2 = ((knot_vector[i+p+1] - u) / denom2) * \
                        nurbs_basis(i+1, p-1, u, knot_vector)
            return term1 + term2

    x_curve = 0
    y_curve = 0
    z_curve = 0
    sum_weights = 0

    for i in range(n + 1):
        N_i = nurbs_basis(i, degree, s, knot_vector)
        wN_i = N_i * weights[i]
        x_curve += wN_i * control_points[i, 0]
        y_curve += wN_i * control_points[i, 1]
        z_curve += wN_i * control_points[i, 2]
        sum_weights += wN_i

    x_curve = x_curve / sum_weights
    y_curve = y_curve / sum_weights
    z_curve = z_curve / sum_weights

    return x_curve, y_curve, z_curve

# Generate the NURBS curve
weights = [1, 1, 1, 1]
x_curve, y_curve, z_curve = nurbs_gen(s, control_points, weights)

# Define the Frenet-Serret Frame
dx_ds = sp.diff(x_curve, s)
dy_ds = sp.diff(y_curve, s)
dz_ds = sp.diff(z_curve, s)
T = sp.Matrix([dx_ds, dy_ds, dz_ds])
T_norm = sp.sqrt(T.dot(T))
T = T / T_norm

dT_ds = sp.diff(T, s)
N_norm = sp.sqrt(dT_ds.dot(dT_ds))
N = dT_ds / N_norm

B = T.cross(N)

# Define Rotation Angle based on curve parameter s
AoA = a_AoA * s**4 + b_AoA * s**3 + c_AoA * s**2 + d_AoA * s + e_AoA

# Define the rotation matrix
rotation_matrix = sp.Matrix([
    [sp.cos(AoA), -sp.sin(AoA)],
    [sp.sin(AoA), sp.cos(AoA)]
])

# Rotate and scale the 2D shape
XY_rotated = rotation_matrix * sp.Matrix([x_2D, y_2D])
X_rotated = XY_rotated[0]
Y_rotated = XY_rotated[1]

scale_x = a_scX * s**4 + b_scX * s**3 + c_scX * s**2 + d_scX * s + e_scX
scale_y = a_scY * s**4 + b_scY * s**3 + c_scY * s**2 + d_scY * s + e_scY

X_rotated_scaled = X_rotated * scale_x
Y_rotated_scaled = Y_rotated * scale_y

# Transform the 2D shape along the curve
C = sp.Matrix([x_curve, y_curve, z_curve])

X_final = C[0] + X_rotated_scaled * N[0] + Y_rotated_scaled * B[0]
Y_final = C[1] + X_rotated_scaled * N[1] + Y_rotated_scaled * B[1]
Z_final = C[2] + X_rotated_scaled * N[2] + Y_rotated_scaled * B[2]


# Define a mapping for the Heaviside function
heaviside = lambda x: np.heaviside(x, 0.5)

# Lambdify the expressions
X_func = sp.lambdify((s, t), X_final, modules=['numpy', {'Heaviside': heaviside}])
Y_func = sp.lambdify((s, t), Y_final, modules=['numpy', {'Heaviside': heaviside}])
Z_func = sp.lambdify((s, t), Z_final, modules=['numpy', {'Heaviside': heaviside}])

# Create meshgrid for s and t
s_vals = np.linspace(s_domain[0], s_domain[1], s_resolution)
t_vals = np.linspace(t_domain[0], t_domain[1], t_resolution)
s_mesh, t_mesh = np.meshgrid(s_vals, t_vals)

# Evaluate the functions on the meshgrid
X_vals_base = X_func(s_mesh, t_mesh)
Y_vals_base = Y_func(s_mesh, t_mesh)
Z_vals_base = Z_func(s_mesh, t_mesh)

# remove NaN values completely from the list
X_vals_base = X_vals_base[~np.isnan(X_vals_base)]
Y_vals_base = Y_vals_base[~np.isnan(Y_vals_base)]
Z_vals_base = Z_vals_base[~np.isnan(Z_vals_base)]

# Adjust blade geometry to ensure connection with hub
# Extend the blade slightly into the hub
hub_extension = 0.05  # Adjust as needed
Z_vals_base += hub_extension * hub_length


# Prepare vertices and faces for the blade mesh
def create_mesh_from_surface(X_vals, Y_vals, Z_vals):
    # Flatten the arrays
    points = np.column_stack((X_vals.ravel(), Y_vals.ravel(), Z_vals.ravel()))

    # Generate faces using Delaunay triangulation
    tri = Delaunay(points[:, :2])
    faces = np.hstack(
        (np.full((tri.simplices.shape[0], 1), 3), tri.simplices)
    )
    return points, faces

# Initialize final mesh
final_mesh = None

# Create the hub mesh using PyVista
hub_length_z = hub_length  # Total length of the hub along Z-axis
hub = pv.Cylinder(
    center=(0, 0, 0), direction=(0, 0, 1), radius=hub_radius,
    height=hub_length_z, resolution=100
)
# Shift the hub to match the coordinate system
hub.translate((0, 0, -hub_length_z / 2))
hub_mesh = hub.triangulate()

# Generate blades and perform boolean union
for i in range(num_blades):
    theta_rot = i * (2 * np.pi / num_blades)
    # Rotate the blade coordinates
    X_vals_rot = (
        X_vals_base * np.cos(theta_rot) - Y_vals_base * np.sin(theta_rot)
    )
    Y_vals_rot = (
        X_vals_base * np.sin(theta_rot) + Y_vals_base * np.cos(theta_rot)
    )
    Z_vals_rot = Z_vals_base

    # Create blade mesh
    points, faces = create_mesh_from_surface(
        X_vals_rot, Y_vals_rot, Z_vals_rot
    )
    blade_mesh = pv.PolyData(points, faces)

    # Perform boolean union with the hub
    if final_mesh is None:
        # First blade union with hub
        final_mesh = hub_mesh.boolean_union(blade_mesh)
    else:
        # Union the next blade
        final_mesh = final_mesh.boolean_union(blade_mesh)

# Visualization
plotter = pv.Plotter()
plotter.add_mesh(final_mesh, color='silver', show_edges=False)
plotter.show()

# Export the final mesh to STL
final_mesh.save('propeller.stl')
