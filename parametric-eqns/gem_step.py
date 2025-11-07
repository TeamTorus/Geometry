# Import necessary libraries
import numpy as np
import sympy as sp
import time

# --- CadQuery Imports ---
# Use the high-level CadQuery library as requested
import cadquery as cq

# --- Your existing imports ---
from nurbs_gen import nurbs_gen
# (Other imports from your original file are no longer needed)

# symbolic variables
t, s = sp.symbols('t s', real=True)

# ---------------------------------------- PARAMETERS ----------------------------------------
# Domain and resolution
s_domain = [0, 1]  # Domain for the curve parameter s
t_domain = [0, 2]  # Domain for the shape parameter t

# --- CAD Generation Parameters ---
s_resolution_cad = 25  # Number of "ribs" along the blade's length
t_resolution_cad = 60  # Number of points to define each airfoil "rib"

# Modifiable parameters
hub_radius = 5  # Radius of the cylindrical hub
hub_length = 20  # Length of the cylindrical hub
num_blades = 3  # Number of blades

# Airfoil Params
m = 0.04
p = 0.4
thickness = .75

# Centerline Params
loc_ctrl_point2 = [2, -5, 25]
loc_ctrl_point3 = [5, 0.75, 30]
blade_vector = [12, 1.5]  # offset between the two endpoints

# Angle of Attack
a_AoA = 0
b_AoA = 0
c_AoA = 0
d_AoA = 0.2 * np.pi
e_AoA = np.pi

# Scaling Params
a_scX = 1
b_scX = 0
c_scX = -2
d_scX = 2.5
e_scX = 3

a_scY = 0
b_scY = 0
c_scY = -1
d_scY = 0
e_scY = 2

apply_thickness_normal = False

print("Parameters set. Starting symbolic definition...")

# ---------------------------------------- PT 1: Define the 2D shape (Unchanged) ----------------------------------------
# ... (This entire section is identical to your code) ...

yt = 5 * thickness * (0.2969 * sp.sqrt(t) - 0.1260 * t - 0.3516 * t**2 + 0.2843 * t**3 - 0.1036 * t**4)
yc = sp.Piecewise(((m / p**2) * (2 * p * t - t**2), t <= p), ((m / (1 - p)**2) * ((1 - 2 * p) + 2 * p * t - t**2), t > p))

if apply_thickness_normal:
    dyc_dt = sp.Piecewise((2*m/p**2 * (p - t), t <= p), (2*m/(1 - p)**2 * (p - t), t > p))
    theta_c = sp.atan(dyc_dt)
else:
    theta_c = 0

yu = yc + yt * sp.cos(theta_c)
yl_1 = yc - yt * sp.cos(theta_c)
xu = t - yt * sp.sin(theta_c)
# Define the lower-surface x-coordinate function
xl_1 = t + yt * sp.sin(theta_c) 
# Now, re-parameterize it 
xl = xl_1.subs(t, (2 - t))
yl = yl_1.subs(t, (2 - t))
y_2D = yu * sp.Heaviside(1 - t) + yl * sp.Heaviside(t - 1)
x_2D = xu * sp.Heaviside(1 - t) + xl * sp.Heaviside(t - 1) - 0.5
AoA = a_AoA * s**4 + b_AoA * s**3 + c_AoA * s**2 + d_AoA * s + e_AoA
rotation_matrix = sp.Matrix([[sp.cos(AoA), -sp.sin(AoA)], [sp.sin(AoA), sp.cos(AoA)]])
XY_rotated = rotation_matrix * sp.Matrix([x_2D, y_2D]) 
X_rotated = XY_rotated[0]
Y_rotated = XY_rotated[1]
scale_x = a_scX * s**4 + b_scX * s**3 + c_scX * s**2 + d_scX * s + e_scX
scale_y = a_scY * s**4 + b_scY * s**3 + c_scY * s**2 + d_scY * s + e_scY
scale_x = 0.2 + (scale_x - 0.2) * sp.Heaviside(0.99 - s) * sp.Heaviside(s - 0.01)
scale_y = 0.2 + (scale_y - 0.2) * sp.Heaviside(0.99 - s) * sp.Heaviside(s - 0.01)
X_rotated_scaled = X_rotated * scale_x
Y_rotated_scaled = Y_rotated * scale_y

# -------------------------------------- PT 2: Set Parameter Domains (Simplified) --------------------------------------

t_vals_cad = np.linspace(t_domain[0], t_domain[1], t_resolution_cad, endpoint=False) # <-- ADD 'endpoint=False'
s_vals_cad = np.linspace(s_domain[0], s_domain[1], s_resolution_cad)

# ------------------------------------------ PT 3: Create the 3D curve (Unchanged) ------------------------------------------
# ... (This entire section is identical to your code) ...

scale0 = float(max(scale_x.subs(s, 0), scale_y.subs(s, 0))) 
scale1 = float(max(scale_x.subs(s, 1), scale_y.subs(s, 1))) 
scale = max(scale0, scale1)
inset_ratio = 4/8
blade_hub_radius = inset_ratio * hub_radius
ctrl_point1 = [blade_hub_radius, 0, hub_length / 2 - 1]
disp_theta = blade_vector[0] / (2 * np.pi * blade_hub_radius) * 2 * np.pi
ctrl_point4 = [blade_hub_radius * np.cos(disp_theta), blade_hub_radius * np.sin(disp_theta), -1 * blade_vector[1]]
disp_theta2 = loc_ctrl_point2[0] / (2 * np.pi * blade_hub_radius) * 2 * np.pi
loc2_radius = blade_hub_radius + loc_ctrl_point2[2]
ctrl_point2 = [loc2_radius * np.cos(disp_theta2), loc2_radius * np.sin(disp_theta2), -1 * loc_ctrl_point2[1]]
disp_theta3 = loc_ctrl_point3[0] / (2 * np.pi * blade_hub_radius) * 2 * np.pi
loc3_radius = blade_hub_radius + loc_ctrl_point3[2]
ctrl_point3 = [loc3_radius * np.cos(disp_theta3), loc3_radius * np.sin(disp_theta3), -1 * loc_ctrl_point3[1]]
control_points = np.array([ctrl_point1, ctrl_point2, ctrl_point3, ctrl_point4])
weights = [1, 1, 1, 1]
x_curve, y_curve, z_curve = nurbs_gen(s, control_points, weights, to_plot=False)

print("Symbolic centerline defined.")

# -------------------------------------- PT 4: Frenet-Serret frame along curve (Unchanged) --------------------------------------
# ... (This entire section is identical to your code) ...

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
C = sp.Matrix([x_curve, y_curve, z_curve])
X_final = C[0] + X_rotated_scaled * N[0] + Y_rotated_scaled * B[0]
Y_final = C[1] + X_rotated_scaled * N[1] + Y_rotated_scaled * B[1]
Z_final = C[2] + X_rotated_scaled * N[2] + Y_rotated_scaled * B[2]

print("Symbolic surface equations defined. Lambdifying...")

# -------------------------------------- PT 5: Lambdify Symbolic Expressions (Unchanged) --------------------------------------

C_func = sp.lambdify(s, C, 'numpy')
N_func = sp.lambdify(s, N, 'numpy')
B_func = sp.lambdify(s, B, 'numpy')
X_airfoil_func = sp.lambdify((s, t), X_rotated_scaled, modules=['numpy', {'Heaviside': np.heaviside}])
Y_airfoil_func = sp.lambdify((s, t), Y_rotated_scaled, modules=['numpy', {'Heaviside': np.heaviside}])

print("Lambdifying complete. Starting B-Rep generation with CadQuery...")

# -------------------------------------- PT 6: Build B-Rep Blade Solid (NEW: CadQuery) --------------------------------------
# -------------------------------------- PT 6: Build B-Rep Blade Solid (NEW: CadQuery) --------------------------------------
start_time = time.time()

# This list will hold all the 3D "rib" wires
rib_wires = []

for s_val in s_vals_cad:
    # 1. Get centerline position and Frenet frame vectors for this 's'
    C_vec = C_func(s_val).flatten().astype(float)
    N_vec = N_func(s_val).flatten().astype(float)
    B_vec = B_func(s_val).flatten().astype(float)

    # 2. Generate the 2D airfoil points in the local (N, B) plane
    x_local = X_airfoil_func(s_val, t_vals_cad)
    y_local = Y_airfoil_func(s_val, t_vals_cad)

    # 3. Create a list of 3D CadQuery Vectors
    points_3d = []
    for (x_l, y_l) in zip(x_local, y_local):
        # Final 3D point = C + x_local * N + y_local * B
        pt_3d = C_vec + x_l * N_vec + y_l * B_vec
        points_3d.append(cq.Vector(pt_3d[0], pt_3d[1], pt_3d[2]))

    # 4. NOW we manually close the loop (re-add this line)
    # This appends Point_A to the end of [Point_A, Point_B, ..., Point_Z]
    points_3d.append(points_3d[0])

    # 5. Create a smooth B-Spline edge from the 3D points
    spline_edge = cq.Edge.makeSpline(points_3d)
    
    # 6. Use a Workplane as a context to "promote" the edge to a wire
    wire = cq.Workplane().add(spline_edge).wires().val()

    # 7. Add the wire to our list for lofting
    rib_wires.append(wire)

# 8. Create the lofted solid from all the rib wires
blade_solid = cq.Solid.makeLoft(rib_wires)

print(f"Blade solid created in {time.time() - start_time:.2f} seconds.")

# --------------------------------------- PT 7: Build Hub and Assemble Propeller (NEW: CadQuery) ---------------------------------------

# 1. Create the hub solid
# We'll create it centered at the origin, matching your math
hub_solid = cq.Solid.makeCylinder(
    hub_radius,
    hub_length,
    pnt=cq.Vector(0, 0, -hub_length / 2), # Start point
    dir=cq.Vector(0, 0, 1)                # Extrusion direction
)

# 2. Assemble the propeller using Boolean Fuses (Unions)
# Start with the hub as the base
propeller_solid = hub_solid

print(f"Hub created. Assembling {num_blades} blades...")

for i in range(num_blades):
    # 3. Define the rotation angle
    angle_deg = i * (360.0 / num_blades)
    
    # 4. Create a rotated copy of the blade
    blade_copy = blade_solid.rotate(
        (0, 0, 0),       # rotation center
        (0, 0, 1),       # rotation axis
        angle_deg        # angle in degrees
    )
    
    # 5. Fuse the rotated blade to the main propeller solid
    propeller_solid = propeller_solid.fuse(blade_copy)
    
    print(f"Blade {i+1}/{num_blades} fused.")

print("Assembly complete.")

# --------------------------------------- PT 8: Export to STEP File (NEW: CadQuery) ---------------------------------------

output_filename = "toroidal_propeller.step"

# Use CadQuery's built-in exporter
try:
    cq.exporters.export(
        propeller_solid,
        output_filename,
        "STEP"
    )
    print(f"Successfully exported STEP file to: {output_filename}")

except Exception as e:
    print(f"Failed to export STEP file: {e}")