import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.mplot3d import Axes3D

# --- Plot Formatting for LaTeX ---
plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'legend.fontsize': 10,
    'figure.dpi': 150,
    'font.family': 'serif' # Matches LaTeX styling well
})

# --- Core Airfoil Function ---
def get_airfoil(t, thickness=0.25):
    """
    Generates NACA 4-digit symmetric airfoil based on parameter t in [0, 2].
    t in [0, 1] is the upper surface.
    t in (1, 2] is the lower surface.
    """
    x = np.where(t <= 1.0, t, 2.0 - t)
    yt = 5 * thickness * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1036 * x**4)
    y = np.where(t <= 1.0, yt, -yt)
    x = x - 0.5
    return x, y

def apply_transform(x, y, scale_x=1.0, scale_y=1.0, alpha_deg=0.0):
    """Applies scaling and rotation (angle of attack) to the 2D airfoil."""
    alpha = np.radians(alpha_deg)
    x_rot = x * np.cos(alpha) - y * np.sin(alpha)
    y_rot = x * np.sin(alpha) + y * np.cos(alpha)
    x_scaled = x_rot * scale_x
    y_scaled = y_rot * scale_y
    return x_scaled, y_scaled

# ==========================================
# FIGURE 1: Parameter t Colormap
# ==========================================
def plot_figure_1():
    fig, ax = plt.subplots(figsize=(6, 4))
    t_vals = np.linspace(0, 2, 500)
    x, y = get_airfoil(t_vals)
    
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
    norm = plt.Normalize(t_vals.min(), t_vals.max())
    lc = LineCollection(segments, cmap='plasma', norm=norm, linewidth=3)
    lc.set_array(t_vals)
    line = ax.add_collection(lc)
    
    cbar = fig.colorbar(line, ax=ax)
    cbar.set_label('Parameter $t$')
    
    ax.set_xlim(-0.6, 0.6)
    ax.set_ylim(-0.2, 0.2)
    ax.set_aspect('equal')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.grid(True, linestyle=':')
    fig.tight_layout()

# ==========================================
# FIGURES 2a, 2b, 2c: Airfoil Transformations
# ==========================================
def plot_figure_2():
    t_vals = np.linspace(0, 2, 200)
    x_base, y_base = get_airfoil(t_vals)
    
    # 2a: X-Scaling
    fig_a, ax_a = plt.subplots(figsize=(5, 4))
    scales_x = [0.5, 1.0, 1.5]
    for sx in scales_x:
        x_t, y_t = apply_transform(x_base, y_base, scale_x=sx)
        ax_a.plot(x_t, y_t, label=f'$S_x = {sx}$')
    ax_a.set_aspect('equal')
    ax_a.set_xlim(-0.8, 0.8)
    ax_a.set_ylim(-0.4, 0.4)
    ax_a.grid(True, linestyle=':')
    ax_a.legend()
    ax_a.set_xlabel('X')
    ax_a.set_ylabel('Y')
    fig_a.tight_layout()
    
    # 2b: Y-Scaling
    fig_b, ax_b = plt.subplots(figsize=(5, 4))
    scales_y = [0.5, 1.0, 2.0]
    for sy in scales_y:
        x_t, y_t = apply_transform(x_base, y_base, scale_y=sy)
        ax_b.plot(x_t, y_t, label=f'$S_y = {sy}$')
    ax_b.set_aspect('equal')
    ax_b.set_xlim(-0.8, 0.8)
    ax_b.set_ylim(-0.4, 0.4)
    ax_b.grid(True, linestyle=':')
    ax_b.legend()
    ax_b.set_xlabel('X')
    ax_b.set_ylabel('Y')
    fig_b.tight_layout()
    
    # 2c: Angle of Attack
    fig_c, ax_c = plt.subplots(figsize=(5, 4))
    alphas = [-15, 0, 15, 30]
    for a in alphas:
        x_t, y_t = apply_transform(x_base, y_base, alpha_deg=a)
        ax_c.plot(x_t, y_t, label=f'$\\alpha = {a}^\\circ$')
    ax_c.set_aspect('equal')
    ax_c.set_xlim(-0.8, 0.8)
    ax_c.set_ylim(-0.4, 0.4)
    ax_c.grid(True, linestyle=':')
    ax_c.legend()
    ax_c.set_xlabel('X')
    ax_c.set_ylabel('Y')
    fig_c.tight_layout()

# ==========================================
# Centerline & Frenet-Serret Helper
# ==========================================
def get_centerline(s_vals):
    P = np.array([
        [0.0, 0.0, 0.0],
        [0.6, 0.6, -0.5],
        [0.8, 1.5, -1.0],
        [0.0, 2.0, -1.5]
    ])
    C = np.zeros((len(s_vals), 3))
    for i, s in enumerate(s_vals):
        C[i] = ((1-s)**3 * P[0] + 
                3*(1-s)**2 * s * P[1] + 
                3*(1-s) * s**2 * P[2] + 
                s**3 * P[3])
    return C, P

def get_frenet_frame(C):
    dC = np.gradient(C, axis=0)
    T = dC / np.linalg.norm(dC, axis=1)[:, np.newaxis]
    dT = np.gradient(T, axis=0)
    N = dT / np.linalg.norm(dT, axis=1)[:, np.newaxis]
    B = np.cross(T, N)
    return T, N, B

# ==========================================
# FIGURE 3: 3D Centerline & Perpendicular Airfoils
# ==========================================
def plot_figure_3():
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection='3d')
    
    s_vals = np.linspace(0, 1, 100)
    C, P = get_centerline(s_vals)
    T, N, B = get_frenet_frame(C)
    
    ax.plot(C[:,0], C[:,1], C[:,2], 'k-', linewidth=2, label='Centerline $C(s)$')
    ax.scatter(P[:,0], P[:,1], P[:,2], c='red', s=50, label='Control Points')
    ax.plot(P[:,0], P[:,1], P[:,2], 'r--', alpha=0.5)
    
    t_vals = np.linspace(0, 2, 50)
    x_base, y_base = get_airfoil(t_vals)
    
    intervals = np.linspace(0, 99, 10).astype(int)
    for idx in intervals:
        x_s = x_base * 0.4
        y_s = y_base * 0.4
        X_final = C[idx, 0] + x_s * N[idx, 0] + y_s * B[idx, 0]
        Y_final = C[idx, 1] + x_s * N[idx, 1] + y_s * B[idx, 1]
        Z_final = C[idx, 2] + x_s * N[idx, 2] + y_s * B[idx, 2]
        ax.plot(X_final, Y_final, Z_final, color='blue', alpha=0.8)

    ax.set_xlabel('Global X')
    ax.set_ylabel('Global Y')
    ax.set_zlabel('Global Z')
    ax.legend()
    ax.set_box_aspect([1, 1, 1])
    fig.tight_layout()

    # plot just centerline
    fig = plt.figure(figsize=(7, 6))
    ax = fig.add_subplot(111, projection='3d')
    
    s_vals = np.linspace(0, 1, 100)
    C, P = get_centerline(s_vals)
    T, N, B = get_frenet_frame(C)
    
    ax.plot(C[:,0], C[:,1], C[:,2], 'k-', linewidth=2, label='Centerline $C(s)$')
    ax.scatter(P[:,0], P[:,1], P[:,2], c='red', s=50, label='Control Points')
    ax.plot(P[:,0], P[:,1], P[:,2], 'r--', alpha=0.5)
    
    ax.set_xlabel('Global X')
    ax.set_ylabel('Global Y')
    ax.set_zlabel('Global Z')
    ax.legend()
    ax.set_box_aspect([1, 1, 1])
    fig.tight_layout()

# ==========================================
# FIGURES 4a & 4b: Polynomials & Resultant Blade
# ==========================================
def plot_figure_4():
    s_vals = np.linspace(0, 1, 100)
    
    # 4a: Polynomials
    fig_a, ax_a = plt.subplots(figsize=(6, 5))
    alpha_s = 2 *(120 * s_vals**4 - 150 * s_vals**3 + 2 * s_vals**2 + 60 * s_vals)
    Sx_s = 15*s_vals**4 - 30 * s_vals**3 + 15 * s_vals**2 + 1.0
    Sy_s = 10 * s_vals**2 - 10 *s_vals + 3
    
    ax_a.plot(s_vals, alpha_s, label='$\\alpha(s)$ [Degrees]', color='red')
    ax_a.plot(s_vals, Sx_s, label='$S_x(s)$ [Scale]', color='green')
    ax_a.plot(s_vals, Sy_s, label='$S_y(s)$ [Scale]', color='blue')
    
    ax_a.set_xlabel('Parameter $s$')
    ax_a.set_ylabel('Magnitude')
    ax_a.grid(True, linestyle=':')
    ax_a.legend()
    fig_a.tight_layout()
    
    # 4b: Resultant Blade
    fig_b = plt.figure(figsize=(6, 6))
    ax_b = fig_b.add_subplot(111, projection='3d')
    C, _ = get_centerline(s_vals)
    T, N, B = get_frenet_frame(C)
    
    t_vals = np.linspace(0, 2, 40)
    x_base, y_base = get_airfoil(t_vals)
    
    S_mesh, T_mesh = np.meshgrid(s_vals, t_vals, indexing='ij')
    X_surf = np.zeros_like(S_mesh)
    Y_surf = np.zeros_like(S_mesh)
    Z_surf = np.zeros_like(S_mesh)
    
    # Pre-calculate full mesh for surface plotting
    for i in range(len(s_vals)):
        x_t, y_t = apply_transform(x_base, y_base, scale_x=Sx_s[i]*0.3, scale_y=Sy_s[i]*0.3, alpha_deg=alpha_s[i])
        X_surf[i, :] = C[i, 0] + x_t * N[i, 0] + y_t * B[i, 0]
        Y_surf[i, :] = C[i, 1] + x_t * N[i, 1] + y_t * B[i, 1]
        Z_surf[i, :] = C[i, 2] + x_t * N[i, 2] + y_t * B[i, 2]
    
    # Plot the colored surface
    ax_b.plot_surface(X_surf, Y_surf, Z_surf, cmap='viridis', edgecolor='none', alpha=0.5)
    
    # New Loop: Add "ribs" every 0.1 of s (every 10 indices)
    rib_intervals = np.linspace(0, 99, 10).astype(int)
    for idx in rib_intervals:
        # Re-calculate or reuse the already computed transformed airfoil at this s
        x_t, y_t = apply_transform(x_base, y_base, scale_x=Sx_s[idx]*0.3, scale_y=Sy_s[idx]*0.3, alpha_deg=alpha_s[idx])
        
        X_rib = C[idx, 0] + x_t * N[idx, 0] + y_t * B[idx, 0]
        Y_rib = C[idx, 1] + x_t * N[idx, 1] + y_t * B[idx, 1]
        Z_rib = C[idx, 2] + x_t * N[idx, 2] + y_t * B[idx, 2]
        
        # Plot rib as a thin black line to contrast with surface cmap
        ax_b.plot(X_rib, Y_rib, Z_rib, color='black', linewidth=2, alpha=1)

    ax_b.set_xlabel('Global X')
    ax_b.set_ylabel('Global Y')
    ax_b.set_zlabel('Global Z')
    ax_b.set_box_aspect([1, 1, 1])
    fig_b.tight_layout()

# ==========================================
# FIGURE 5: Local Non-Euclidean Frame on Hub
# ==========================================
def plot_coordinate_system():
    # Make figure slightly wider to accommodate legend
    fig = plt.figure(figsize=(9, 6))
    ax = fig.add_subplot(111, projection='3d')
    
    # Adjust subplot position to leave wide margin on right for legend
    plt.subplots_adjust(right=0.75)
    
    # Cylinder parameters
    radius = 2.0
    length = 6.0
    
    # Generate cylinder mesh
    theta = np.linspace(0, 2 * np.pi, 50)
    z = np.linspace(-length/2, length/2, 50)
    Theta, Z = np.meshgrid(theta, z)
    X = radius * np.cos(Theta)
    Y = radius * np.sin(Theta)
    
    # Plot central hub (cylinder)
    ax.plot_surface(X, Y, Z, color='lightgray', alpha=0.3, edgecolor='none')
    
    
    # Define origin of local frame L on the surface of the cylinder
    origin_theta = np.radians(-45) # offset slightly for visual separation
    origin_z = length/2 - 1.0
    origin = np.array([radius * np.cos(origin_theta), radius * np.sin(origin_theta), origin_z])
    ax.scatter(*origin, color='black', s=50, label='Local Origin $P_1$')

    # --- New Visualization of Curvilinear Axes ---
    
    # 1. +x_L axis curve (wraps horizontally around circumference - tangential)
    # We define this as a change in theta at a constant R and Z.
    theta_arc = np.linspace(origin_theta, origin_theta + np.radians(45), 50)
    X_arc = radius * np.cos(theta_arc)
    Y_arc = radius * np.sin(theta_arc)
    Z_arc = np.ones_like(theta_arc) * origin_z
    
    # Label curve as curvilinear axis
    ax.plot(X_arc, Y_arc, Z_arc, color='red', linewidth=3, label='Circumferential Axis ($+x_L$)')
    
    # Arrowhead at end of curvilinear axis
    ax.quiver(X_arc[-2], Y_arc[-2], Z_arc[-2], 
              X_arc[-1]-X_arc[-2], Y_arc[-1]-Y_arc[-2], 0, 
              color='red', length=0.4, normalize=True, arrow_length_ratio=1)

    # 2. +y_L axis curve (vertical line down the cylinder lateral face - tangent to -Z)
    z_line = np.linspace(origin_z, origin_z - 2.5, 50)
    X_line = np.ones_like(z_line) * origin[0]
    Y_line = np.ones_like(z_line) * origin[1]
    
    ax.plot(X_line, Y_line, z_line, color='green', linewidth=3, label='Downward Axis ($+y_L$)')
    
    # Arrowhead
    ax.quiver(origin[0], origin[1], z_line[-2], 
              0, 0, z_line[-1]-z_line[-2], 
              color='green', length=0.4, normalize=True, arrow_length_ratio=1)

    # 3. +z_L axis path (straight radial normal line)
    # This remains Euclidean as it simply points out radially normal to the face.
    r_path = np.linspace(radius, radius + 2.0, 50)
    X_rad = r_path * np.cos(origin_theta)
    Y_rad = r_path * np.sin(origin_theta)
    Z_rad = np.ones_like(r_path) * origin_z
    
    ax.plot(X_rad, Y_rad, Z_rad, color='blue', linewidth=3, label='Radial Normal Axis ($+z_L$)')
    
    # Arrowhead
    ax.quiver(X_rad[-2], Y_rad[-2], Z_rad[-2], 
              X_rad[-1]-X_rad[-2], Y_rad[-1]-Y_rad[-2], 0, 
              color='blue', length=0.4, normalize=True, arrow_length_ratio=1)

    # Add text labels along the curvilinear paths
    ax.text(X_arc[-1], Y_arc[-1], Z_arc[-1], r'$+x_L$', color='red', fontsize=12, fontweight='bold')
    ax.text(origin[0], origin[1], z_line[-1]-0.2, r'$+y_L$', color='green', fontsize=12, fontweight='bold')
    ax.text(X_rad[-1], Y_rad[-1], Z_rad[-1], r'$+z_L$', color='blue', fontsize=12, fontweight='bold')
    
    # Reference Hub Axis
    ax.plot([0, 0], [0, 0], [-length/2, length/2], 'k-', alpha=0.5, linewidth=3, label='Hub Axis (Global Z)')
    
    # Standard global arrows at origin for comparison
    # ax.quiver(radius, 0, origin_z, 0, 1, 0, color='r', alpha=0.3, linestyle='--', length=1, label='Global +Y')
    
    ax.set_xlabel('Global X')
    ax.set_ylabel('Global Y')
    ax.set_zlabel('Global Z')
    
    # Moved Legend to the far right margin
    ax.legend(loc=(1.05, 0.5), fontsize=9)
    
    # Equal aspect ratio
    ax.set_box_aspect([1, 1, length/(2*radius)])
    # fig.tight_layout() # tight_layout messes up the legend placement we just did manually

# --- Execution ---
if __name__ == '__main__':
    plot_figure_1()
    plot_figure_2()
    plot_figure_3()
    plot_figure_4()
    plot_coordinate_system()
    plt.show()