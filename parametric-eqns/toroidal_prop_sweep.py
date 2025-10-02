# Import necessary libraries
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import pyvista as pv
import os
import pickle
import itertools
from scipy.integrate import quad
from scipy.interpolate import interp1d

from nurbs_gen import nurbs_gen
from screw_prop_gen import generate_equivalent_screw_propeller
from union_with_bpy import union_stl_files

# symbolic variables
t, s = sp.symbols('t s', real=True)

# Create output directory
output_dir = "sweep_outputs"
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Load existing metadata or create new
metadata_file = os.path.join(output_dir, "metadata.pkl")
if os.path.exists(metadata_file):
    with open(metadata_file, 'rb') as f:
        metadata = pickle.load(f)
else:
    metadata = {}

# Fixed parameters
s_domain = [0, 1]
t_domain = [0, 2]
s_resolution = 100
t_resolution = 100
hub_resolution = 50

normalize_blade_mesh = False
apply_thickness_normal = False
close_cylinder = True
plot_matplotlib = False  # Disabled for batch processing

hub_radius = 5
hub_length = 20
num_blades = 3

# Parameter sweeps - adjust ranges as needed
m_values = [0.02, 0.04, 0.06]  # NACA camber
p_values = [0.3, 0.4, 0.5, 0.6]  # NACA camber position
thickness_values = [0.3, 0.5, 0.8]  # Airfoil thickness

# Control point variations (relative to base values)
loc_ctrl_point2_base = [3, 3, 10]
loc_ctrl_point3_base = [7, 6, 15]
blade_vector_base = [8, 8]

# Variation factors for control points
ctrl_point_variations = [0.75, 1.0, 1.25, 1.5]  # Scale factors
blade_vector_variations = [0.75, 1.0, 1.5]

# Angle of attack coefficients
c_AoA_values = [0, 0.5, 1.0]  # Quadratic term
d_AoA_values = [0.5 * np.pi, 1.0 * np.pi, 1.5 * np.pi]  # Linear term
e_AoA_values = [np.pi]  # Constant term

# Fixed AoA coefficients
a_AoA = 0
b_AoA = 0

# Fixed scaling parameters
scaling_params = {
    'a_scX': 1, 'b_scX': 0, 'c_scX': 0, 'd_scX': 1, 'e_scX': 2,
    'a_scY': 0, 'b_scY': 0, 'c_scY': 0, 'd_scY': 1, 'e_scY': 1
}

def generate_propeller(param_dict, filename):
    """Generate a single propeller with given parameters"""
    
    # Extract parameters
    m = param_dict['m']
    p = param_dict['p']
    thickness = param_dict['thickness']
    loc_ctrl_point2 = param_dict['loc_ctrl_point2']
    loc_ctrl_point3 = param_dict['loc_ctrl_point3']
    blade_vector = param_dict['blade_vector']
    c_AoA = param_dict['c_AoA']
    d_AoA = param_dict['d_AoA']
    e_AoA = param_dict['e_AoA']
    
    try:
        # PT 1: Define the 2D shape
        yt = 5 * thickness * (0.2969 * sp.sqrt(t) - 0.1260 * t - 0.3516 * t**2 + 0.2843 * t**3 - 0.1036 * t**4)
        yc = sp.Piecewise(((m / p**2) * (2 * p * t - t**2), t <= p), ((m / (1 - p)**2) * ((1 - 2 * p) + 2 * p * t - t**2), t > p))

        theta_c = 0  # Not applying thickness normal

        # Apply camber + thickness
        yu = yc + yt * sp.cos(theta_c)
        yl_1 = yc - yt * sp.cos(theta_c)
        xu = t - yt * sp.sin(theta_c)
        xl = (2 - t) + yt * sp.sin(theta_c)
        yl = yl_1.subs(t, (2 - t))

        # Use activations for upper/lower
        y_2D = yu * sp.Heaviside(1 - t) + yl * sp.Heaviside(t - 1)
        x_2D = xu * sp.Heaviside(1 - t) + xl * sp.Heaviside(t - 1) - 0.5

        # Angle of attack
        AoA = a_AoA * s**4 + b_AoA * s**3 + c_AoA * s**2 + d_AoA * s + e_AoA
        rotation_matrix = sp.Matrix([[sp.cos(AoA), -sp.sin(AoA)], [sp.sin(AoA), sp.cos(AoA)]])
        XY_rotated = rotation_matrix * sp.Matrix([x_2D, y_2D])
        X_rotated = XY_rotated[0]
        Y_rotated = XY_rotated[1]

        # Scaling
        scale_x = scaling_params['a_scX'] * s**4 + scaling_params['b_scX'] * s**3 + scaling_params['c_scX'] * s**2 + scaling_params['d_scX'] * s + scaling_params['e_scX']
        scale_y = scaling_params['a_scY'] * s**4 + scaling_params['b_scY'] * s**3 + scaling_params['c_scY'] * s**2 + scaling_params['d_scY'] * s + scaling_params['e_scY']

        X_rotated_scaled = X_rotated * scale_x
        Y_rotated_scaled = Y_rotated * scale_y

        # PT 2: Parameter discretization (no normalization for speed)
        t_vals = np.linspace(t_domain[0], t_domain[1], t_resolution)
        s_vals = np.linspace(s_domain[0], s_domain[1], s_resolution)

        # PT 3: Create the 3D curve
        scale0 = max(scale_x.subs(s, 0), scale_y.subs(s, 0))
        scale1 = max(scale_x.subs(s, 1), scale_y.subs(s, 1))
        scale = max(scale0, scale1)
        inset_ratio = 1 - min(scale * 1/4, 1/2)
        blade_hub_radius = inset_ratio * hub_radius

        # Control points
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

        # PT 4: Frenet-Serret frame
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

        # PT 5: Generate meshes
        def generate_cylinder_mesh(radius, length, resolution):
            theta_vals = np.linspace(0, 2 * np.pi, resolution)
            z_vals = np.linspace(-length / 2, length / 2, resolution)
            theta_mesh, z_mesh = np.meshgrid(theta_vals, z_vals)
            X = radius * np.cos(theta_mesh)
            Y = radius * np.sin(theta_mesh)
            Z = z_mesh
            return X, Y, Z

        X_hub, Y_hub, Z_hub = generate_cylinder_mesh(hub_radius, hub_length, hub_resolution)

        R = sp.sqrt(X_final**2 + Y_final**2)
        Theta = sp.atan2(Y_final, X_final)

        X_prop = np.zeros((t_resolution, s_resolution * num_blades))
        Y_prop = np.zeros((t_resolution, s_resolution * num_blades))
        Z_prop = np.zeros((t_resolution, s_resolution * num_blades))

        for i in range(num_blades):
            Loc_Theta = Theta + i * 2 * np.pi / num_blades
            X_rotated_blade = R * sp.cos(Loc_Theta)
            Y_rotated_blade = R * sp.sin(Loc_Theta)

            X_func = sp.lambdify((s, t), X_rotated_blade, modules=['numpy', {'Heaviside': np.heaviside}])
            Y_func = sp.lambdify((s, t), Y_rotated_blade, modules=['numpy', {'Heaviside': np.heaviside}])
            Z_func = sp.lambdify((s, t), Z_final, modules=['numpy', {'Heaviside': np.heaviside}])

            s_mesh, t_mesh = np.meshgrid(s_vals, t_vals)
            X_vals = X_func(s_mesh, t_mesh)
            Y_vals = Y_func(s_mesh, t_mesh)
            Z_vals = Z_func(s_mesh, t_mesh)

            X_prop[:, i * s_resolution:(i + 1) * s_resolution] = X_vals
            Y_prop[:, i * s_resolution:(i + 1) * s_resolution] = Y_vals
            Z_prop[:, i * s_resolution:(i + 1) * s_resolution] = Z_vals

        # Create PyVista meshes
        def create_mesh_from_grids(X, Y, Z):
            points = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
            mesh = pv.PolyData(points)
            nx, ny = X.shape
            faces = []
            for i in range(nx - 1):
                for j in range(ny - 1):
                    idx0 = i * ny + j
                    idx1 = idx0 + 1
                    idx2 = idx0 + ny + 1
                    idx3 = idx0 + ny
                    faces.extend([4, idx0, idx1, idx2, idx3])
            faces = np.array(faces)
            mesh.faces = faces
            mesh = mesh.triangulate()
            return mesh

        def close_cylinder_mesh(X, Y, Z):
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

            top_center = np.array([[0, 0, Z.max()]])
            bottom_center = np.array([[0, 0, Z.min()]])
            points = np.vstack((points, top_center, bottom_center))
            idx_top_center = len(points) - 2
            idx_bottom_center = len(points) - 1

            idx_top = (nx - 1) * ny + np.arange(ny)
            idx_bottom = 0 * ny + np.arange(ny)

            for j in range(ny):
                idx0 = idx_top[j]
                idx1 = idx_top[(j + 1) % ny]
                faces.extend([3, idx0, idx1, idx_top_center])

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

        # Clean meshes
        mesh_hub_m = mesh_hub.fill_holes(100).clean()
        mesh_blades_m = mesh_blades.fill_holes(100).clean()

        # Save temporary files
        temp_hub = f'temp_hub_{filename[:-4]}.stl'
        temp_blades = f'temp_blades_{filename[:-4]}.stl'
        mesh_hub_m.save(temp_hub)
        mesh_blades_m.save(temp_blades)

        # Boolean union
        output_path = os.path.join(output_dir, filename)
        success = union_stl_files(temp_hub, temp_blades, output_path)

        # Clean up temporary files
        if os.path.exists(temp_hub):
            os.remove(temp_hub)
        if os.path.exists(temp_blades):
            os.remove(temp_blades)

        return success

    except Exception as e:
        print(f"Error generating {filename}: {e}")
        return False

# Main sweep loop
print("Starting parameter sweep...")
count = 0
failed_count = 0

# Create parameter combinations
for m in m_values:
    for p in p_values:
        for thickness in thickness_values:
            for ctrl_var in ctrl_point_variations:
                for blade_var in blade_vector_variations:
                    for c_AoA in c_AoA_values:
                        for d_AoA in d_AoA_values:
                            for e_AoA in e_AoA_values:
                                
                                # Apply variations to control points
                                loc_ctrl_point2 = [x * ctrl_var for x in loc_ctrl_point2_base]
                                loc_ctrl_point3 = [x * ctrl_var for x in loc_ctrl_point3_base]
                                blade_vector = [x * blade_var for x in blade_vector_base]
                                
                                # Create parameter dictionary
                                param_dict = {
                                    's_domain': s_domain, 't_domain': t_domain, 
                                    's_resolution': s_resolution, 't_resolution': t_resolution,
                                    'normalize_blade_mesh': normalize_blade_mesh, 
                                    'apply_thickness_normal': apply_thickness_normal,
                                    'close_cylinder': close_cylinder, 'plot_matplotlib': plot_matplotlib,
                                    'hub_radius': hub_radius, 'hub_length': hub_length, 
                                    'num_blades': num_blades,
                                    'm': m, 'p': p, 'thickness': thickness,
                                    'loc_ctrl_point2': loc_ctrl_point2, 
                                    'loc_ctrl_point3': loc_ctrl_point3, 
                                    'blade_vector': blade_vector,
                                    'a_AoA': a_AoA, 'b_AoA': b_AoA, 'c_AoA': c_AoA, 
                                    'd_AoA': d_AoA, 'e_AoA': e_AoA,
                                    **scaling_params
                                }
                                
                                # Generate unique filename
                                filename = f"prop_{count:03d}_m{m:.2f}_p{p:.1f}_t{thickness:.1f}_cv{ctrl_var:.1f}_bv{blade_var:.1f}_c{c_AoA:.1f}_d{d_AoA:.2f}_e{e_AoA:.2f}.stl"
                                
                                print(f"Generating {count+1}: {filename}")
                                
                                success = generate_propeller(param_dict, filename)
                                
                                if success:
                                    # Update metadata
                                    metadata[filename] = param_dict
                                    
                                    # Save metadata after each successful generation
                                    with open(metadata_file, 'wb') as f:
                                        pickle.dump(metadata, f)
                                    
                                    print(f"  ✓ Successfully generated {filename}")
                                else:
                                    failed_count += 1
                                    print(f"  ✗ Failed to generate {filename}")
                                
                                count += 1
                                
                                # Break if we've generated enough (optional limit)
                                if count >= 300:  # Safety limit
                                    break
                            if count >= 300:
                                break
                        if count >= 300:
                            break
                    if count >= 300:
                        break
                if count >= 300:
                    break
            if count >= 300:
                break
        if count >= 300:
            break
    if count >= 300:
        break

print(f"\nSweep completed!")
print(f"Total attempted: {count}")
print(f"Successfully generated: {count - failed_count}")
print(f"Failed: {failed_count}")
print(f"Output directory: {output_dir}")
print(f"Metadata saved to: {metadata_file}")
