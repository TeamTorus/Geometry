import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

from nurbs_gen import nurbs_gen, nurbs_basis, plot_nurbs_curve

# Define symbols
s, t = sp.symbols('s t')

# Define domain and resolution parameters
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
loc_ctrl_point2 = np.array([3, 3, 10])
loc_ctrl_point3 = np.array([7, 6, 8])
blade_vector = np.array([8, 8])  # Offset between the two endpoints

# Angle of Attack Params
a_AoA = 0
b_AoA = 0
c_AoA = 0
d_AoA = 2 * sp.pi
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

# Define the 2D Shape (NACA airfoil)
x_2D = t * sp.Heaviside(1 - t) + (2 - t) * sp.Heaviside(t - 1)

yt = 5 * thickness * (0.2969 * sp.sqrt(t) - 0.1260 * t - 0.3516 * t**2 + 0.2843 * t**3 - 0.1036 * t**4)
yc = m / p**2 * (2 * p * t - t**2) * sp.Heaviside(p - t) + m / (1 - p)**2 * ((1 - 2*p) + 2*p*t - t**2) * sp.Heaviside(t - p)
yu = yt + yc
yl_1 = yc - yt
yl = yl_1.subs(t, (2 - t))

# Activations for upper/lower surface
y_2D = yu * sp.Heaviside(1 - t) + yl * sp.Heaviside(t - 1)

# Define the 3D Curve
blade_hub_radius = hub_radius - hub_radius / 8

# First point of the control points
ctrl_point1 = np.array([blade_hub_radius, 0, hub_length / 2 - 1])

# Last control point offset by blade_vector
disp_theta = blade_vector[0] / (2 * np.pi * blade_hub_radius) * 2 * np.pi
ctrl_point4 = np.array([blade_hub_radius * np.cos(disp_theta), blade_hub_radius * np.sin(disp_theta), -1 * blade_vector[1]])

# Second control point
disp_theta2 = loc_ctrl_point2[0] / (2 * np.pi * blade_hub_radius) * 2 * np.pi
loc2_radius = blade_hub_radius + loc_ctrl_point2[2]
ctrl_point2 = np.array([loc2_radius * np.cos(disp_theta2), loc2_radius * np.sin(disp_theta2), -1 * loc_ctrl_point2[1]])

# Third control point
disp_theta3 = loc_ctrl_point3[0] / (2 * np.pi * blade_hub_radius) * 2 * np.pi
loc3_radius = blade_hub_radius + loc_ctrl_point3[2]
ctrl_point3 = np.array([loc3_radius * np.cos(disp_theta3), loc3_radius * np.sin(disp_theta3), -1 * loc_ctrl_point3[1]])

# Control points
control_points = np.array([ctrl_point1, ctrl_point2, ctrl_point3, ctrl_point4])

# Use the NURBS generator to calculate the 3D curve
x_curve, y_curve, z_curve = nurbs_gen(s, control_points, [1, 1, 1, 1], False)

plot_nurbs_curve(x_curve, y_curve, z_curve, control_points)

# Define the Frenet-Serret Frame
dx_ds = sp.diff(x_curve, s)
dy_ds = sp.diff(y_curve, s)
dz_ds = sp.diff(z_curve, s)
T = sp.Matrix([dx_ds, dy_ds, dz_ds])  # Tangent vector
T = T / T.norm()

dT_ds = sp.diff(T, s)
N = dT_ds / dT_ds.norm()  # Normal vector

B = T.cross(N)  # Binormal vector

# Define the rotation angle based on curve parameter s
AoA = a_AoA * s**4 + b_AoA * s**3 + c_AoA * s**2 + d_AoA * s + e_AoA

# Rotation matrix (counterclockwise in the local x-y plane)
rotation_matrix = sp.Matrix([[sp.cos(AoA), -sp.sin(AoA)], [sp.sin(AoA), sp.cos(AoA)]])

# Rotate the 2D shape symbolically
XY_rotated = rotation_matrix * sp.Matrix([x_2D, y_2D])
X_rotated = XY_rotated[0]
Y_rotated = XY_rotated[1]

# Define scaling/stretching functions
scale_x = a_scX * s**4 + b_scX * s**3 + c_scX * s**2 + d_scX * s + e_scX  # Parametric scaling for x
scale_y = a_scY * s**4 + b_scY * s**3 + c_scY * s**2 + d_scY * s + e_scY  # Parametric scaling for y

# Apply scaling in the local frame (after rotation)
X_rotated_scaled = X_rotated * scale_x
Y_rotated_scaled = Y_rotated * scale_y

# Define the parametric position along the curve
C = sp.Matrix([x_curve, y_curve, z_curve])

# Express the final 3D coordinates of the airfoil after extrusion
X_final = C[0] + X_rotated_scaled * N[0] + Y_rotated_scaled * B[0]
Y_final = C[1] + X_rotated_scaled * N[1] + Y_rotated_scaled * B[1]
Z_final = C[2] + X_rotated_scaled * N[2] + Y_rotated_scaled * B[2]

# Discretizing and plotting the final surface

# Create a meshgrid for s and t
s_vals = np.linspace(s_domain[0], s_domain[1], s_resolution)
t_vals = np.linspace(t_domain[0], t_domain[1], t_resolution)
S, T = np.meshgrid(s_vals, t_vals)

# Lambdify the symbolic expressions for numerical evaluation
X_final_func = sp.lambdify((s, t), X_final, "numpy")
Y_final_func = sp.lambdify((s, t), Y_final, "numpy")
Z_final_func = sp.lambdify((s, t), Z_final, "numpy")

# Evaluate the functions on the grid
X_vals = X_final_func(S, T)
Y_vals = Y_final_func(S, T)
Z_vals = Z_final_func(S, T)

# Plot the surface
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X_vals, Y_vals, Z_vals, cmap='viridis')

# Set labels and title
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title('3D Propeller Blade Surface')

plt.show()
