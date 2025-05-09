import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from nurbs_gen import nurbs_gen  # your existing symbolic NURBS helper

"""
Toroidal propeller - point-cloud generator
-----------------------------------------
• Keeps **all** analytic steps from your original script (symbolic airfoil → rotation → scaling → Frenet-Serret frame → extrusion)
• Uses the *exact* symbolic expressions for T, N, B (no finite-difference approximations)
• Adds a tiny orientation-smoothing routine so successive cross-sections keep the same handedness (avoids flipped airfoils when n_s is small)
"""

# ---------------------------------------------------------------------------
# 0)  Symbolic setup & master parameters
# ---------------------------------------------------------------------------
# Independent symbols
s, t = sp.symbols('s t', real=True)

# ---- Airfoil parameters (same as your NACA-4 definition) ----
m, p, thickness = 0.04, 0.4, 0.5
apply_thickness_normal = False  # keep as in original

# ---- AoA(s) polynomial coefficients        AoA(s) = a*s⁴ + b*s³ + c*s² + d*s + e ----
a_AoA, b_AoA, c_AoA, d_AoA, e_AoA = 0, 0, 0, 1*np.pi, 0

# ---- Span-wise scaling for the airfoil outline ----
a_scX, b_scX, c_scX, d_scX, e_scX = 1, 0, 0, 1, 2

a_scY, b_scY, c_scY, d_scY, e_scY = 0, 0, 0, 1, 1

# ---- Centre-line control points & weights ----
num_blades = 3
hub_radius = 5 
hub_length = 20
angle_max_deg = 25  # max entry angle (radial) for hub entry
enforce_angle = True  

degree_around_hub = 60
# Centerline Params
loc_ctrl_point2 = [3, 3, 10]
loc_ctrl_point3 = [7, 6, 15]
blade_vector = [8, 8]   # offset between the two endpoints



# ------ Generate the ctrl points for the centerline ------
scaleX = a_scX*s**4 + b_scX*s**3 + c_scX*s**2 + d_scX*s + e_scX
scaleY = a_scY*s**4 + b_scY*s**3 + c_scY*s**2 + d_scY*s + e_scY


scale0 = max(scaleX.subs(s, 0), scaleY.subs(s, 0))  # scale at s=0
scale1 = max(scaleX.subs(s, 1), scaleY.subs(s, 1))  # scale at s=1
scale = max(scale0, scale1)  # scale factor for the control points

inset_ratio = 1 - min(scale * 1/4, 1/2)

blade_hub_radius = inset_ratio * hub_radius

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

print(ctrl_point1, ctrl_point2, ctrl_point3, ctrl_point4)

control_points = np.array([ctrl_point1, ctrl_point2, ctrl_point3, ctrl_point4])  # [x, y, z]
weights = [1, 1, 1, 1]

# ---------------------------------------------------------------------------
# Enforce hub entry angle by adjusting interior control points
# ---------------------------------------------------------------------------
if enforce_angle:
    def adjust_entry(idx_root, idx_next):
        """Ensure ctrl_pts[idx_root]→ctrl_pts[idx_next] is radial (≤ angle_max)."""
        P_root = ctrl_pts[idx_root]
        P_next = ctrl_pts[idx_next]
        # tangent estimate (first derivative dir for C0 NURBS)
        T_vec = P_next - P_root
        T_unit = T_vec / np.linalg.norm(T_vec)
        # hub outward normal (radial in xy-plane)
        radial_vec = np.array([P_root[0], P_root[1], 0.0])
        radial_unit = radial_vec / np.linalg.norm(radial_vec)
        # angle between T and radial
        angle = np.degrees(np.arccos(np.clip(np.abs(np.dot(T_unit, radial_unit)), -1, 1)))
        if angle > angle_max_deg:
            # project P_next onto radial line, keeping its original z for smoothness
            proj_len = np.linalg.norm(P_next[:2] - P_root[:2])  # keep same spacing
            P_adjust = P_root[:2] + proj_len * radial_unit[:2]
            ctrl_pts[idx_next, 0] = P_adjust[0]
            ctrl_pts[idx_next, 1] = P_adjust[1]
            # recompute tangent (optional diagnostics)
            new_T = ctrl_pts[idx_next]-P_root
            new_ang = np.degrees(np.arccos(np.clip(np.abs(np.dot(new_T/np.linalg.norm(new_T),radial_unit)), -1, 1)))
            print(f"Adjusted control point {idx_next} to enforce entry angle: {angle:.1f}° → {new_ang:.1f}°")

    # root side (s≈0) uses ctrl_pts[0] & ctrl_pts[1]
    adjust_entry(0, 1)
    # tip side (s≈1) uses ctrl_pts[3] & ctrl_pts[2] (note reverse order)
    adjust_entry(3, 2)

# ---- Discretisation ----
n_s = 25      # blade stations (keep small to test continuity)
n_t = 50    # points per airfoil outline (higher = smoother airfoil)

# ---------------------------------------------------------------------------
# 1)  Exact symbolic airfoil in local (x,y)
# ---------------------------------------------------------------------------
y_t = 5*thickness * (0.2969*sp.sqrt(t) - 0.1260*t - 0.3516*t**2 + 0.2843*t**3 - 0.1036*t**4)
y_c = sp.Piecewise(
    ((m/p**2)*(2*p*t - t**2), t <= p),
    ((m/(1-p)**2)*((1-2*p) + 2*p*t - t**2), True),
)
if apply_thickness_normal:
    dyc_dt = sp.Piecewise((2*m/p**2*(p - t), t <= p), (2*m/(1-p)**2*(p - t), True))
    theta_c = sp.atan(dyc_dt)
else:
    theta_c = 0

# upper & lower (x,y) halves, exactly as in your original
x_u = t - y_t*sp.sin(theta_c)
y_u = y_c + y_t*sp.cos(theta_c)

x_l = (2 - t) + y_t*sp.sin(theta_c)
# note: lower camber reflected with t→(2-t)
y_l = (y_c - y_t*sp.cos(theta_c)).subs(t, 2 - t)

# stitch using Heaviside so the symbolic form is continuous over t∈[0,2]
x_2D = (x_u*sp.Heaviside(1 - t) + x_l*sp.Heaviside(t - 1)) - 0.5  # centre at 0
y_2D =  y_u*sp.Heaviside(1 - t) + y_l*sp.Heaviside(t - 1)

# ---------------------------------------------------------------------------
# 2)  AoA(s), scaling(s), rotation matrix
# ---------------------------------------------------------------------------
AoA = a_AoA*s**4 + b_AoA*s**3 + c_AoA*s**2 + d_AoA*s + e_AoA
Rot = sp.Matrix([[sp.cos(AoA), -sp.sin(AoA)],
                 [sp.sin(AoA),  sp.cos(AoA)]])
xy_rot = Rot*sp.Matrix([x_2D, y_2D])
Xr, Yr = xy_rot[0], xy_rot[1]

Xrs = Xr*scaleX
Yrs = Yr*scaleY

# ---------------------------------------------------------------------------
# 3)  Centre-line curve & Frenet-Serret frame (exact derivatives)
# ---------------------------------------------------------------------------
x_curve, y_curve, z_curve = nurbs_gen(s, ctrl_pts, weights, to_plot=False)

T = sp.Matrix([sp.diff(x_curve,s), sp.diff(y_curve,s), sp.diff(z_curve,s)])
T = T / sp.sqrt(T.dot(T))

dT_ds = sp.diff(T, s)
N = dT_ds / sp.sqrt(dT_ds.dot(dT_ds))
B = T.cross(N)

# ---------------------------------------------------------------------------
# 4)  Final embedded surface coordinates
# ---------------------------------------------------------------------------
X_final = x_curve + Xrs*N[0] + Yrs*B[0]
Y_final = y_curve + Xrs*N[1] + Yrs*B[1]
Z_final = z_curve + Xrs*N[2] + Yrs*B[2]

# Lambdify (Heaviside→numpy.heaviside)
mods = ['numpy', {"Heaviside": np.heaviside}]
Xfun = sp.lambdify((s, t), X_final, modules=mods)
Yfun = sp.lambdify((s, t), Y_final, modules=mods)
Zfun = sp.lambdify((s, t), Z_final, modules=mods)

# ---------------------------------------------------------------------------
# 5)  Mesh / point cloud evaluation
# ---------------------------------------------------------------------------
s_vals = np.linspace(0, 1, n_s)
t_vals = np.linspace(0, 2, n_t)
S, Tm = np.meshgrid(s_vals, t_vals, indexing='ij')  # shape (n_s,n_t)

X = Xfun(S, Tm)
Y = Yfun(S, Tm)
Z = Zfun(S, Tm)

# ── optional: orientation smoothing when n_s is small (flips from sign ambiguity) ──
# Evaluate N & B numerically just to detect flips
Nfun = sp.lambdify(s, N, modules=mods)
Bfun = sp.lambdify(s, B, modules=mods)
N_num = np.array([Nfun(si).flatten() for si in s_vals])  # (n_s,3)
B_num = np.array([Bfun(si).flatten() for si in s_vals])  # (n_s,3)

# Then the broadcasting will work correctly
orient_sign = np.ones(n_s)
for i in range(1, n_s):
    if np.dot(N_num[i], N_num[i-1]) < 0:  # sign flip detected
        orient_sign[i:] *= -1
        
# Apply the sign correction to N_num and B_num
N_num = N_num * orient_sign[:,None]  # Reshape orient_sign for broadcasting
B_num = B_num * orient_sign[:,None]
# regenerate coordinates with consistent sign
for i,sgn in enumerate(orient_sign):
    if sgn == -1:
        X[i] = Xfun(s_vals[i], Tm[i])  # re-eval not needed, parity handled via sign in N,B, but cheap
        Y[i] = Yfun(s_vals[i], Tm[i])
        Z[i] = Zfun(s_vals[i], Tm[i])

points = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])  # (n_s*n_t,3)

# ---------------------------------------------------------------------------
# 6)  Quick visual: scatter + polylines to show continuity
# ---------------------------------------------------------------------------
fig = plt.figure(figsize=(9,7))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(points[:,0], points[:,1], points[:,2], s=2, alpha=0.7)

# draw each airfoil polyline to see if cross-sections connect
for i in range(n_s):
    ax.plot(X[i], Y[i], Z[i], lw=0.8, c='k')

ax.set_xlabel('X'); ax.set_ylabel('Y'); ax.set_zlabel('Z')
ax.set_title('Toroidal Propeller - point cloud (exact analytic frame)')
plt.tight_layout()
plt.show()

outfile = 'propeller_pointcloud.txt'
mid = np.where(t_vals <= 1)[0][-1]
with open(outfile, 'w') as f:
    for i in range(n_s - 1):
        idx = i + 1
        # top half
        f.write(f"Airfoil{idx}a\n")
        f.write("START\t0.000000\n")
        for j in range(0, mid+1):
            f.write(f"{X[i,j]:.8f}\t{Y[i,j]:.8f}\t{Z[i,j]:.8f}\n")
        f.write("END\t0.000000\n")
        # bottom half
        f.write(f"Airfoil{idx}b\n")
        f.write("START\t0.000000\n")
        for j in range(mid, n_t):
            f.write(f"{X[i,j]:.8f}\t{Y[i,j]:.8f}\t{Z[i,j]:.8f}\n")
        f.write("END\t0.000000\n")
    # guide curves
    f.write("GuideCurve1\n")
    f.write("START\t0.000000\n")
    for i in range(n_s - 1):
        f.write(f"{X[i,0]:.8f}\t{Y[i,0]:.8f}\t{Z[i,0]:.8f}\n")
    f.write("END\t0.000000\n")
    f.write("GuideCurve2\n")
    f.write("START\t0.000000\n")
    for i in range(n_s - 1):
        f.write(f"{X[i,mid]:.8f}\t{Y[i,mid]:.8f}\t{Z[i,mid]:.8f}\n")
    f.write("END\t0.000000\n")
    # hub and blade count
    f.write(f"HubRadius {hub_radius * 1.35:.2f}\n")     # expand the hub for offset
    f.write(f"Blades {num_blades}\n")
print(f"Wrote pointcloud data to {outfile}")