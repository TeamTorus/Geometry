import sympy as sp
import numpy as np

# Define symbolic variables
t, s = sp.symbols('t s')

# Define parameters
hub_radius = 5  # Radius of the cylindrical hub
hub_length = 20  # Length of the cylindrical hub
num_blades = 3  # Number of blades

# Airfoil parameters
m = 0.02
p = 0.4
thickness = 0.5

# Define control points and other parameters
loc_ctrl_point2 = [3, 3, 10]
loc_ctrl_point3 = [7, 6, 8]
blade_vector = [8, 8]

# Define the NACA airfoil shape equations
x_2D = t * sp.Heaviside(1 - t) + (2 - t) * sp.Heaviside(t - 1)

yt = 5 * thickness * (0.2969 * sp.sqrt(t) - 0.1260 * t - 0.3516 * t**2 + 0.2843 * t**3 - 0.1036 * t**4)
yc = m / p**2 * (2 * p * t - t**2) * sp.Heaviside(p - t) + \
     m / (1 - p)**2 * ((1 - 2*p) + 2*p*t - t**2) * sp.Heaviside(t - p)
yu = yt + yc

# Domain shift for lower surface
yl_1 = yc - yt
yl = yl_1.subs(t, (2 - t))

# Use heaviside activations for upper and lower surface
y_2D = yu * sp.Heaviside(1 - t) + yl * sp.Heaviside(t - 1)


# Define the NURBS basis function recursively
def nurbs_basis(i, p, u, knot_vector):
    if p == 0:
        return sp.Heaviside(u - knot_vector[i]) * sp.Heaviside(knot_vector[i + 1] - u)
    else:
        if knot_vector[i + p] == knot_vector[i]:
            c1 = 0
        else:
            c1 = (u - knot_vector[i]) / (knot_vector[i + p] - knot_vector[i]) * nurbs_basis(i, p - 1, u, knot_vector)

        if knot_vector[i + p + 1] == knot_vector[i + 1]:
            c2 = 0
        else:
            c2 = (knot_vector[i + p + 1] - u) / (knot_vector[i + p + 1] - knot_vector[i + 1]) * nurbs_basis(i + 1, p - 1, u, knot_vector)

        return c1 + c2

# Define control points and weights
control_points = np.array([
    [hub_radius - hub_radius / 8, 0, hub_length / 2 - 1],  # control point 1
    # Add other control points here (converted as in MATLAB)
])

weights = [1, 1, 1, 1]
knot_vector = [0, 0, 0, 0, 1, 1, 1, 1]
degree = 3

# Calculate NURBS curve
x_curve = 0
y_curve = 0
z_curve = 0

for i in range(len(control_points)):
    N_i = nurbs_basis(i, degree, s, knot_vector)
    x_curve += N_i * weights[i] * control_points[i, 0]
    y_curve += N_i * weights[i] * control_points[i, 1]
    z_curve += N_i * weights[i] * control_points[i, 2]

# Weighted sum of basis functions
sum_weights = sum(nurbs_basis(i, degree, s, knot_vector) * weights[i] for i in range(len(control_points)))

# Final NURBS curve equations
x_curve /= sum_weights
y_curve /= sum_weights
z_curve /= sum_weights

x_curve = sp.simplify(x_curve)
y_curve = sp.simplify(y_curve)
z_curve = sp.simplify(z_curve)

# Derivatives for tangent vector (T)
dx_ds = sp.diff(x_curve, s)
dy_ds = sp.diff(y_curve, s)
dz_ds = sp.diff(z_curve, s)

T = sp.Matrix([dx_ds, dy_ds, dz_ds])
T = T / T.norm()

# Derivative of tangent for normal vector (N)
dT_ds = sp.diff(T, s)
N = dT_ds / dT_ds.norm()

# Binormal vector (B) as the cross product of T and N
B = T.cross(N)

# Define symbolic angle of attack
a_AoA, b_AoA, c_AoA, d_AoA, e_AoA = 0, 0, 0, 2 * np.pi, 0
AoA = a_AoA * s**4 + b_AoA * s**3 + c_AoA * s**2 + d_AoA * s + e_AoA

# Define rotation matrix
rotation_matrix = sp.Matrix([
    [sp.cos(AoA), -sp.sin(AoA)],
    [sp.sin(AoA), sp.cos(AoA)]
])

# Rotate 2D airfoil shape
XY_rotated = rotation_matrix * sp.Matrix([x_2D, y_2D])
X_rotated = XY_rotated[0]
Y_rotated = XY_rotated[1]

# Define scaling functions
a_scX, b_scX, c_scX, d_scX, e_scX = 0, 0, 0, 1, 2
a_scY, b_scY, c_scY, d_scY, e_scY = 0, 0, 0, 1, 2

scale_x = a_scX * s**4 + b_scX * s**3 + c_scX * s**2 + d_scX * s + e_scX
scale_y = a_scY * s**4 + b_scY * s**3 + c_scY * s**2 + d_scY * s + e_scY

# Apply scaling
X_rotated_scaled = X_rotated * scale_x
Y_rotated_scaled = Y_rotated * scale_y

# Final parametric position along the curve
C = sp.Matrix([x_curve, y_curve, z_curve])

# Final 3D coordinates after extrusion
X_final = C[0] + X_rotated_scaled * N[0] + Y_rotated_scaled * B[0]
Y_final = C[1] + X_rotated_scaled * N[1] + Y_rotated_scaled * B[1]
Z_final = C[2] + X_rotated_scaled * N[2] + Y_rotated_scaled * B[2]


# plot the 3D curve
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Create a 3D plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Plot the 3D curve for the blade, using the parametric equations for X, Y, Z, in terms of s and t
s_vals = np.linspace(0, 1, 100)
t_vals = np.linspace(0, 2, 200)

X_vals = np.array([[X_final.subs({s: s_val, t: t_val}) for s_val in s_vals] for t_val in t_vals], dtype=np.float64)
Y_vals = np.array([[Y_final.subs({s: s_val, t: t_val}) for s_val in s_vals] for t_val in t_vals], dtype=np.float64)
Z_vals = np.array([[Z_final.subs({s: s_val, t: t_val}) for s_val in s_vals] for t_val in t_vals], dtype=np.float64)

ax.plot_surface(X_vals, Y_vals, Z_vals, cmap='viridis')

# Set the aspect ratio of the 3D plot to be equal
ax.set_box_aspect([1, 1, 1])

# Set the axis labels
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

# Show the plot
plt.show()