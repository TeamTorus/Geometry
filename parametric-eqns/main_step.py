# Requires: pip install cadquery  (or)  pip install cadquery-ocp
# (CadQuery brings OCP/OpenCascade; we'll use CadQuery-level APIs to build NURBS & export STEP.)

import numpy as np
import sympy as sp
from scipy.integrate import quad
from scipy.interpolate import interp1d
import cadquery as cq

from nurbs_gen import nurbs_gen  # unchanged math for cubic NURBS centerline

# ----------------------------- Params (kept the same) -----------------------------
t, s = sp.symbols('t s', real=True)

s_domain = [0, 1]
t_domain = [0, 2]
s_resolution = 100
t_resolution = 100

normalize_blade_mesh = False
apply_thickness_normal = False
close_cylinder = True   # for the hub solid
plot_matplotlib = False # no plotting here

hub_radius = 5
hub_length = 20
num_blades = 3

# Airfoil (NACA-4) params
m = 0.04
p = 0.4
thickness = 0.75

# Centerline shaping params (same semantics as your code)
loc_ctrl_point2 = [2, -5, 25]
loc_ctrl_point3 = [5, 0.75, 30]
blade_vector   = [12, 1.5]  # arc displacement and z-offset between endpoints

# Angle of Attack α(s) = a s^4 + b s^3 + c s^2 + d s + e
a_AoA = 0
b_AoA = 0
c_AoA = 0
d_AoA = 0.2 * np.pi
e_AoA = np.pi

# Scaling laws (same 4th-degree polynomials)
a_scX, b_scX, c_scX, d_scX, e_scX = 1, 0, -2, 2.5, 3
a_scY, b_scY, c_scY, d_scY, e_scY = 0, 0, -1, 0, 2

# ----------------------------- Airfoil (2D) – same equations -----------------------------
yt = 5 * thickness * (0.2969 * sp.sqrt(t) - 0.1260 * t - 0.3516 * t**2 + 0.2843 * t**3 - 0.1036 * t**4)
yc = sp.Piecewise(((m / p**2) * (2 * p * t - t**2), t <= p),
                  ((m / (1 - p)**2) * ((1 - 2 * p) + 2 * p * t - t**2), t > p))

if apply_thickness_normal:
    dyc_dt = sp.Piecewise((2*m/p**2 * (p - t), t <= p),
                          (2*m/(1 - p)**2 * (p - t), t > p))
    theta_c = sp.atan(dyc_dt)
else:
    theta_c = 0

yu  = yc + yt * sp.cos(theta_c)
yl1 = yc - yt * sp.cos(theta_c)

xu = t - yt * sp.sin(theta_c)
xl = (2 - t) + yt * sp.sin(theta_c)          # 2 - t shift

yl = yl1.subs(t, (2 - t))                    # lower half on (1,2]

# Heaviside assembly (single parameter t ∈ [0,2])
y_2D = yu * sp.Heaviside(1 - t) + yl * sp.Heaviside(t - 1)
x_2D = (xu * sp.Heaviside(1 - t) + xl * sp.Heaviside(t - 1)) - 0.5

# AoA(s), rotation, scaling (same as your math)
AoA = a_AoA * s**4 + b_AoA * s**3 + c_AoA * s**2 + d_AoA * s + e_AoA
R = sp.Matrix([[sp.cos(AoA), -sp.sin(AoA)],
               [sp.sin(AoA),  sp.cos(AoA)]])

XY_rot = R * sp.Matrix([x_2D, y_2D])
X_rot, Y_rot = XY_rot[0], XY_rot[1]

scale_x = a_scX * s**4 + b_scX * s**3 + c_scX * s**2 + d_scX * s + e_scX
scale_y = a_scY * s**4 + b_scY * s**3 + c_scY * s**2 + d_scY * s + e_scY

# Clamp scaling to avoid extremes near ends (same heavy-side window)
scale_x = 0.2 + (scale_x - 0.2) * sp.Heaviside(0.99 - s) * sp.Heaviside(s - 0.01)
scale_y = 0.2 + (scale_y - 0.2) * sp.Heaviside(0.99 - s) * sp.Heaviside(s - 0.01)

X_rs = X_rot * scale_x
Y_rs = Y_rot * scale_y

# Optional arc-length normalization for t sampling (kept same behavior)
if normalize_blade_mesh:
    dyu_dt = sp.diff(yu, t)
    dyl_dt = sp.diff(yl, t)
    dxu_dt = sp.diff(xu, t)
    dxl_dt = sp.diff(xl, t)

    ds_u = sp.sqrt(dxu_dt**2 + dyu_dt**2)
    ds_l = sp.sqrt(dxl_dt**2 + dyl_dt**2)

    arc_len_u = sp.lambdify(t, ds_u, 'numpy')
    arc_len_l = sp.lambdify(t, ds_l, 'numpy')

    t_fine_u = np.linspace(t_domain[0], t_domain[1]/2, 500)
    t_fine_l = np.linspace(t_domain[1]/2, t_domain[1], 500)

    s_u_vals = [quad(arc_len_u, 0, tv)[0] for tv in t_fine_u]
    s_l_vals = [quad(arc_len_l, 1, tv)[0] + s_u_vals[-1] for tv in t_fine_l]

    t_fine = np.concatenate((t_fine_u, t_fine_l))
    s_fine = np.concatenate((s_u_vals, s_l_vals))
    inv = interp1d(s_fine, t_fine, kind='linear', fill_value='extrapolate')

    t_vals = inv(np.linspace(0, s_fine[-1], t_resolution))
    t_vals = np.clip(t_vals, t_domain[0], t_domain[1])
else:
    t_vals = np.linspace(t_domain[0], t_domain[1], t_resolution)

eps = 1e-3
s_vals = np.linspace(s_domain[0] + eps, s_domain[1] - eps, s_resolution)

# ----------------------------- NURBS centerline (same) -----------------------------
# Compute hub inset radius based on your logic
scale0 = float(max(scale_x.subs(s, 0), scale_y.subs(s, 0)))
scale1 = float(max(scale_x.subs(s, 1), scale_y.subs(s, 1)))
_ = max(scale0, scale1)  # not directly used but kept for parity

inset_ratio = 4/8
blade_hub_radius = inset_ratio * hub_radius

# Control points as before (cylindrical with angular displacements)
ctrl_point1 = [blade_hub_radius, 0, hub_length / 2 - 1]

disp_theta = blade_vector[0] / (2 * np.pi * blade_hub_radius) * 2 * np.pi
ctrl_point4 = [blade_hub_radius * np.cos(disp_theta),
               blade_hub_radius * np.sin(disp_theta),
               -1 * blade_vector[1]]

disp_theta2 = loc_ctrl_point2[0] / (2 * np.pi * blade_hub_radius) * 2 * np.pi
loc2_radius = blade_hub_radius + loc_ctrl_point2[2]
ctrl_point2 = [loc2_radius * np.cos(disp_theta2),
               loc2_radius * np.sin(disp_theta2),
               -1 * loc_ctrl_point2[1]]

disp_theta3 = loc_ctrl_point3[0] / (2 * np.pi * blade_hub_radius) * 2 * np.pi
loc3_radius = blade_hub_radius + loc_ctrl_point3[2]
ctrl_point3 = [loc3_radius * np.cos(disp_theta3),
               loc3_radius * np.sin(disp_theta3),
               -1 * loc_ctrl_point3[1]]

control_points = np.array([ctrl_point1, ctrl_point2, ctrl_point3, ctrl_point4], dtype=float)
weights = [1, 1, 1, 1]

x_curve, y_curve, z_curve = nurbs_gen(s, control_points, weights, to_plot=False)

# ----------------------------- Frenet–Serret frame (same math) -----------------------------
dx_ds = sp.diff(x_curve, s)
dy_ds = sp.diff(y_curve, s)
dz_ds = sp.diff(z_curve, s)

T_vec = sp.Matrix([dx_ds, dy_ds, dz_ds])
T_vec = T_vec / sp.sqrt(T_vec.dot(T_vec))

dT_ds = sp.diff(T_vec, s)
N_vec = dT_ds / sp.sqrt(dT_ds.dot(dT_ds))
B_vec = T_vec.cross(N_vec)

C = sp.Matrix([x_curve, y_curve, z_curve])

# Final 3D surface param equations:
X_final = C[0] + X_rs * N_vec[0] + Y_rs * B_vec[0]
Y_final = C[1] + X_rs * N_vec[1] + Y_rs * B_vec[1]
Z_final = C[2] + X_rs * N_vec[2] + Y_rs * B_vec[2]

# Lambdify for fast numeric eval (Heaviside provided by numpy)
X_func = sp.lambdify((s, t), X_final, modules=['numpy', {'Heaviside': np.heaviside}])
Y_func = sp.lambdify((s, t), Y_final, modules=['numpy', {'Heaviside': np.heaviside}])
Z_func = sp.lambdify((s, t), Z_final, modules=['numpy', {'Heaviside': np.heaviside}])

# Centerline only (for Hub reference & sanity checks)
Cx = sp.lambdify(s, C[0], 'numpy')
Cy = sp.lambdify(s, C[1], 'numpy')
Cz = sp.lambdify(s, C[2], 'numpy')

# ----------------------------- Build CAD sections and loft via CadQuery -----------------------------
def make_section_wire_at_s(sval: float,
                           tol_point=1e-7,
                           tol_close=1e-6) -> cq.Wire:
    # evaluate section curve
    x = X_func(sval, t_vals)
    y = Y_func(sval, t_vals)
    z = Z_func(sval, t_vals)

    # 1) Drop NaN/Inf rows
    arr = np.column_stack([x, y, z]).astype(float)
    msk = np.isfinite(arr).all(axis=1)
    arr = arr[msk]
    if arr.shape[0] < 4:
        raise ValueError(f"Not enough finite points at s={sval:.6f}")

    # 2) Deduplicate consecutive points (kill zero-length spans)
    cleaned = [arr[0]]
    for p in arr[1:]:
        if np.linalg.norm(p - cleaned[-1]) > tol_point:
            cleaned.append(p)
    if len(cleaned) < 4:
        raise ValueError(f"Too few unique points after cleaning at s={sval:.6f}")

    # 3) Build CQ vectors
    pts = [cq.Vector(*p) for p in cleaned]

    # Distance between ends
    close_dist = (pts[-1].sub(pts[0])).Length

    # 4) If ends are essentially the same point, drop the last to avoid zero-length close edges
    if close_dist <= tol_close:
        pts = pts[:-1]
        close_dist = (pts[-1].sub(pts[0])).Length  # recompute for paranoia
        if len(pts) < 4:
            raise ValueError(f"Loop collapsed after end merge at s={sval:.6f}")

    # 5) Try periodic spline first
    try:
        edge = cq.Edge.makeSpline(pts, periodic=True)
        wire = cq.Wire.assembleEdges([edge])
        if wire.isClosed():
            return wire
        # Not closed? explicitly close only if it’s non-degenerate
        if (pts[-1].sub(pts[0])).Length > tol_close:
            closing = cq.Edge.makeLine(pts[-1], pts[0])
            wire = cq.Wire.assembleEdges([edge, closing])
        return wire
    except Exception:
        # 6) Fallback: non-periodic spline + explicit close IF non-degenerate
        try:
            edge = cq.Edge.makeSpline(pts, periodic=False)
            if (pts[-1].sub(pts[0])).Length > tol_close:
                closing = cq.Edge.makeLine(pts[-1], pts[0])
                wire = cq.Wire.assembleEdges([edge, closing])
            else:
                # The loop is already closed geometrically; a single edge may be fine
                wire = cq.Wire.assembleEdges([edge])
            if wire.isClosed():
                return wire
        except Exception:
            pass

        # 7) Last-resort: polyline wire (always robust), then convert to BSpline
        poly = cq.Wire.makePolygon(pts, closed=True)
        try:
            bs = poly.toBSpline()   # smooths the polyline into a BSpline wire
            return bs
        except Exception:
            # If conversion fails, return the polygonal wire (still valid for loft)
            return poly

# Build one blade by lofting many section wires along s
section_wires_one_blade = [make_section_wire_at_s(float(sv)) for sv in s_vals]
blade_shell = cq.Workplane("XY").add(section_wires_one_blade[0])
for w in section_wires_one_blade[1:]:
    blade_shell = blade_shell.add(w)

# Loft through wires to get a smooth NURBS surface/solid
# Workplane.loft() returns a Solid when wires are closed profiles.
blade_solid_one = blade_shell.loft(ruled=False, close=False)

# ----------------------------- Duplicate blades around hub (same rotation logic) -----------------------------
# We rotate the *entire blade solid* by 2*pi/num_blades increments around Z.
blades = blade_solid_one
for i in range(1, num_blades):
    ang = i * 360.0 / num_blades
    blades = blades.union(blade_solid_one.rotate((0, 0, 0), (0, 0, 1), ang))

# ----------------------------- Make hub cylinder and boolean union -----------------------------
# CadQuery cylinder is along +Z, centered at origin by default if we place pnt at (-h/2) for Z-centering.
hub = cq.Solid.makeCylinder(hub_radius, hub_length, pnt=cq.Vector(0, 0, -hub_length/2.0))

# If you want the inlet/outlet end caps, cylinder is already a solid.
assembly = hub.union(blades)

# ----------------------------- Export STEP -----------------------------
cq.exporters.export(assembly, "toroidal_propeller.step")
print("Wrote toroidal_propeller.step")
