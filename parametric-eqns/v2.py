# Import necessary libraries
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Define symbolic variables
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

# Define the 2D Shape
# https://en.wikipedia.org/wiki/NACA_airfoil
x_2D = t * sp.Heaviside(1 - t) + (2 - t) * sp.Heaviside(t - 1)

yt = 5 * thickness * (0.2969 * sp.sqrt(t) - 0.1260 * t - 0.3516 * t**2 + 0.2843 * t**3 - 0.1036 * t**4)
yc = (m / p**2) * (2 * p * t - t**2) * sp.Heaviside(p - t) + (m / (1 - p)**2) * ((1 - 2 * p) + 2 * p * t - t**2) * sp.Heaviside(t - p)
yu = yt + yc

# Domain shift
yl_1 = yc - yt
yl = yl_1.subs(t, (2 - t))

# Use activations for upper/lower
y_2D = yu * sp.Heaviside(1 - t) + yl * sp.Heaviside(t - 1)

# Define the 3D Curve

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

# NURBS generation function
def nurbs_gen(s, control_points, weights, to_plot=False):
    # Degree of the NURBS curve (cubic)
    degree = 3

    # Knot vector for NURBS curve (for degree 3 with clamped ends)
    knot_vector = np.array([0, 0, 0, 0, 1, 1, 1, 1])

    # Number of control points
    n = len(control_points) - 1

    # Define the NURBS basis function recursively
    def nurbs_basis(i, p, u, knot_vector):
        if p == 0:
            return sp.Piecewise((1, (u >= knot_vector[i]) & (u < knot_vector[i+1])), (0, True))
        else:
            denom1 = knot_vector[i+p] - knot_vector[i]
            denom2 = knot_vector[i+p+1] - knot_vector[i+1]
            term1 = 0
            term2 = 0
            if denom1 != 0:
                term1 = ((u - knot_vector[i]) / denom1) * nurbs_basis(i, p-1, u, knot_vector)
            if denom2 != 0:
                term2 = ((knot_vector[i+p+1] - u) / denom2) * nurbs_basis(i+1, p-1, u, knot_vector)
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

    x_curve = sp.simplify(x_curve)
    y_curve = sp.simplify(y_curve)
    z_curve = sp.simplify(z_curve)

    if to_plot:
        s_vals_spline = np.linspace(0, 1, 100)
        x_spline = [float(x_curve.subs(s, s_val)) for s_val in s_vals_spline]
        y_spline = [float(y_curve.subs(s, s_val)) for s_val in s_vals_spline]
        z_spline = [float(z_curve.subs(s, s_val)) for s_val in s_vals_spline]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(x_spline, y_spline, z_spline, 'b-', linewidth=2)
        ax.plot(control_points[:, 0], control_points[:, 1], control_points[:, 2], 'ro')
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('3D NURBS Spline')
        plt.show()

    return x_curve, y_curve, z_curve

# Generate the NURBS curve
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

# Define Rotation Angle based on curve parameter s
AoA = a_AoA * s**4 + b_AoA * s**3 + c_AoA * s**2 + d_AoA * s + e_AoA

# Define the rotation matrix (counterclockwise in the local x-y plane)
rotation_matrix = sp.Matrix([[sp.cos(AoA), -sp.sin(AoA)], [sp.sin(AoA), sp.cos(AoA)]])

# Rotate the 2D shape symbolically
XY_rotated = rotation_matrix * sp.Matrix([x_2D, y_2D])  # Apply rotation
X_rotated = XY_rotated[0]
Y_rotated = XY_rotated[1]

# Define symbolic scaling/stretching functions (parametric by s)
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

# Define a mapping for the Heaviside function
heaviside = lambda x: np.heaviside(x, 0.5)

# Lambdify the expressions with custom mapping
X_func = sp.lambdify((s, t), X_final, modules=['numpy', {'Heaviside': heaviside}])
Y_func = sp.lambdify((s, t), Y_final, modules=['numpy', {'Heaviside': heaviside}])
Z_func = sp.lambdify((s, t), Z_final, modules=['numpy', {'Heaviside': heaviside}])

# Create meshgrid for s and t
s_vals = np.linspace(s_domain[0], s_domain[1], s_resolution)
t_vals = np.linspace(t_domain[0], t_domain[1], t_resolution)
s_mesh, t_mesh = np.meshgrid(s_vals, t_vals)

print("s_mesh:")
print(s_mesh)

# Evaluate the functions on the meshgrid
X_vals = X_func(s_mesh, t_mesh)
Y_vals = Y_func(s_mesh, t_mesh)
Z_vals = Z_func(s_mesh, t_mesh)

# Visualize the 3D blade
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X_vals, Y_vals, Z_vals, cmap='viridis', edgecolor='none')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('3D Blade Visualization')
plt.show()
