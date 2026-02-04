# toroidal_optimize.py
#
# End-to-end parametric toroidal propeller generator + optimizer.
# - Uses your symbolic / CadQuery generation pipeline
# - Optimizes a subset of parameters to match a target STEP/STL
#
# Requirements:
#   pip install cadquery trimesh scipy sympy numpy
#   (plus your existing nurbs_gen module)
#

import os
import io
import time
import tempfile
import signal
import multiprocessing as mp
import threading
import datetime

import numpy as np
import sympy as sp
import cadquery as cq
import trimesh

from scipy.optimize import differential_evolution
from nurbs_gen import nurbs_gen


# Global stop event seen by workers
_STOP = mp.Event()

# Global evaluation counter (shared across processes)
EVAL_COUNTER = mp.Value('i', 0)
START_TIME   = time.time()
LOG_EVERY    = 1          # print every evaluation
SUMMARY_EVERY = 25        # print a summary every N evaluations

def _fmt_seconds(sec):
    return f"{int(sec//3600):02d}:{int((sec%3600)//60):02d}:{int(sec%60):02d}"


# =====================================================================================
# 1. Parametric propeller generator
# =====================================================================================

def build_propeller(params,
                    s_resolution_cad=25,
                    t_resolution_cad=60,
                    verbose=False):
    """
    Build a CadQuery solid for the toroidal propeller using given parameters.

    params: dict with keys like:
        - global_scale
        - hub_radius, hub_length, num_blades
        - m, p, thickness
        - loc_ctrl_point2 (len 3)
        - loc_ctrl_point3 (len 3)
        - blade_vector (len 2)
        - a_AoA, b_AoA, c_AoA, d_AoA, e_AoA
        - a_scX..e_scX, a_scY..e_scY
        - apply_thickness_normal (bool)

    Returns:
        propeller_solid (unscaled) : cq.Solid
    """

    # -------------------- symbolic variables --------------------
    t, s = sp.symbols('t s', real=True)

    # -------------------- parameter domains ---------------------
    s_domain = [0, 1]
    t_domain = [0, 2]

    # -------------------- unpack modifiable parameters ----------
    global_scale = params.get("global_scale", 7.5)  # NOTE: not used here; used in optimizer
    hub_radius   = params.get("hub_radius",   5.0)
    hub_length   = params.get("hub_length",  20.0)
    num_blades   = params.get("num_blades",   3)

    # Airfoil params
    m         = params.get("m", 0.04)
    p         = params.get("p", 0.4)
    thickness = params.get("thickness", 0.75)

    # Centerline params
    loc_ctrl_point2 = params.get("loc_ctrl_point2", [-2, -5, 25])
    loc_ctrl_point3 = params.get("loc_ctrl_point3", [-5, 0.75, 30])
    blade_vector    = params.get("blade_vector",    [-12, 1.5])

    # AoA coefficients
    a_AoA = params.get("a_AoA", 0.0)
    b_AoA = params.get("b_AoA", 0.0)
    c_AoA = params.get("c_AoA", 0.0)
    d_AoA = params.get("d_AoA", np.pi)
    e_AoA = params.get("e_AoA", 0.0)

    # Scaling params
    a_scX = params.get("a_scX", 1.0)
    b_scX = params.get("b_scX", 0.0)
    c_scX = params.get("c_scX", -2.0)
    d_scX = params.get("d_scX", 2.5)
    e_scX = params.get("e_scX", 3.0)

    a_scY = params.get("a_scY", 0.0)
    b_scY = params.get("b_scY", 0.0)
    c_scY = params.get("c_scY", -1.0)
    d_scY = params.get("d_scY", 0.0)
    e_scY = params.get("e_scY", 2.0)

    apply_thickness_normal = params.get("apply_thickness_normal", False)

    # ------------------------ PT 1: Define the 2D shape ------------------------

    yt = 5 * thickness * (
        0.2969 * sp.sqrt(t)
        - 0.1260 * t
        - 0.3516 * t**2
        + 0.2843 * t**3
        - 0.1036 * t**4
    )

    yc = sp.Piecewise(
        ((m / p**2) * (2 * p * t - t**2), t <= p),
        ((m / (1 - p)**2) * ((1 - 2 * p) + 2 * p * t - t**2), t > p)
    )

    if apply_thickness_normal:
        dyc_dt = sp.Piecewise(
            (2*m/p**2 * (p - t), t <= p),
            (2*m/(1 - p)**2 * (p - t), t > p)
        )
        theta_c = sp.atan(dyc_dt)
    else:
        theta_c = 0

    yu = yc + yt * sp.cos(theta_c)
    yl_1 = yc - yt * sp.cos(theta_c)
    xu = t - yt * sp.sin(theta_c)
    # lower-surface x-coordinate function
    xl_1 = t + yt * sp.sin(theta_c)

    # re-parameterize
    xl = xl_1.subs(t, (2 - t))
    yl = yl_1.subs(t, (2 - t))

    y_2D = yu * sp.Heaviside(1 - t) + yl * sp.Heaviside(t - 1)
    x_2D = xu * sp.Heaviside(1 - t) + xl * sp.Heaviside(t - 1) - 0.5

    AoA_expr = (
        a_AoA * s**4
        + b_AoA * s**3
        + c_AoA * s**2
        + d_AoA * s
        + e_AoA
    )

    rotation_matrix = sp.Matrix([
        [sp.cos(AoA_expr), -sp.sin(AoA_expr)],
        [sp.sin(AoA_expr),  sp.cos(AoA_expr)]
    ])

    XY_rotated = rotation_matrix * sp.Matrix([x_2D, y_2D])
    X_rotated = XY_rotated[0]
    Y_rotated = XY_rotated[1]

    scale_x = (
        a_scX * s**4
        + b_scX * s**3
        + c_scX * s**2
        + d_scX * s
        + e_scX
    )

    scale_y = (
        a_scY * s**4
        + b_scY * s**3
        + c_scY * s**2
        + d_scY * s
        + e_scY
    )

    # clamp scaling near ends
    clamp_scaling_ends = -1 \
        + 1/(1 + sp.exp(-200*(s - 0.01))) \
        + 1/(1 + sp.exp(-200*(0.99 - s)))

    scale_x *= clamp_scaling_ends
    scale_y *= clamp_scaling_ends

    X_rotated_scaled = X_rotated * scale_x
    Y_rotated_scaled = Y_rotated * scale_y

    # ---------------------- PT 2: Parameter domains ------------------------

    t_vals_cad = np.linspace(t_domain[0], t_domain[1], t_resolution_cad, endpoint=False)
    s_vals_cad = np.linspace(s_domain[0], s_domain[1], s_resolution_cad, endpoint=False)

    # ---------------------- PT 3: Create the 3D centerline -----------------

    scale0 = float(max(scale_x.subs(s, 0), scale_y.subs(s, 0)))
    scale1 = float(max(scale_x.subs(s, 1), scale_y.subs(s, 1)))
    scale_max = max(scale0, scale1)

    inset_ratio = 4.0 / 8.0

    blade_hub_radius = inset_ratio * hub_radius

    # first control point
    ctrl_point1 = [blade_hub_radius, 0.0, hub_length / 2.0 - 1.0]

    # last control point: offset by blade_vector (circumference + z)
    disp_theta = blade_vector[0] / (2 * np.pi * blade_hub_radius) * 2 * np.pi
    ctrl_point4 = [
        blade_hub_radius * np.cos(disp_theta),
        blade_hub_radius * np.sin(disp_theta),
        -1.0 * blade_vector[1]
    ]

    # middle control points in cylindrical coordinates
    disp_theta2 = loc_ctrl_point2[0] / (2 * np.pi * blade_hub_radius) * 2 * np.pi
    loc2_radius = blade_hub_radius + loc_ctrl_point2[2]
    ctrl_point2 = [
        loc2_radius * np.cos(disp_theta2),
        loc2_radius * np.sin(disp_theta2),
        -1.0 * loc_ctrl_point2[1]
    ]

    disp_theta3 = loc_ctrl_point3[0] / (2 * np.pi * blade_hub_radius) * 2 * np.pi
    loc3_radius = blade_hub_radius + loc_ctrl_point3[2]
    ctrl_point3 = [
        loc3_radius * np.cos(disp_theta3),
        loc3_radius * np.sin(disp_theta3),
        -1.0 * loc_ctrl_point3[1]
    ]

    control_points = np.array([ctrl_point1, ctrl_point2, ctrl_point3, ctrl_point4])

    weights = [1, 1, 1, 1]
    x_curve, y_curve, z_curve = nurbs_gen(s, control_points, weights, to_plot=False)

    if verbose:
        print("Symbolic centerline defined.")

    # ------------------- PT 4: Frenet-Serret frame ------------------------

    dx_ds = sp.diff(x_curve, s)
    dy_ds = sp.diff(y_curve, s)
    dz_ds = sp.diff(z_curve, s)

    T_vec = sp.Matrix([dx_ds, dy_ds, dz_ds])
    T_norm = sp.sqrt(T_vec.dot(T_vec))
    T_vec = T_vec / T_norm

    dT_ds = sp.diff(T_vec, s)
    N_norm = sp.sqrt(dT_ds.dot(dT_ds))
    N_vec = dT_ds / N_norm
    B_vec = T_vec.cross(N_vec)

    C_vec_sym = sp.Matrix([x_curve, y_curve, z_curve])

    X_final = C_vec_sym[0] + X_rotated_scaled * N_vec[0] + Y_rotated_scaled * B_vec[0]
    Y_final = C_vec_sym[1] + X_rotated_scaled * N_vec[1] + Y_rotated_scaled * B_vec[1]
    Z_final = C_vec_sym[2] + X_rotated_scaled * N_vec[2] + Y_rotated_scaled * B_vec[2]

    if verbose:
        print("Symbolic surface equations defined.")

    # ------------------- PT 5: Lambdify ------------------------

    C_func = sp.lambdify(s, C_vec_sym, 'numpy')
    N_func = sp.lambdify(s, N_vec, 'numpy')
    B_func = sp.lambdify(s, B_vec, 'numpy')

    X_airfoil_func = sp.lambdify((s, t), X_rotated_scaled,
                                 modules=['numpy', {'Heaviside': np.heaviside}])
    Y_airfoil_func = sp.lambdify((s, t), Y_rotated_scaled,
                                 modules=['numpy', {'Heaviside': np.heaviside}])

    if verbose:
        print("Lambdifying complete. Starting B-Rep generation with CadQuery...")

    # ------------------- PT 6: Build blade B-Rep ------------------------

    start_time = time.time()
    rib_wires = []

    if verbose:
        print("--- Starting Rib Generation (s_val, num_points, closure_dist) ---")

    for s_val in s_vals_cad:
        # Centerline + frame
        C_arr = C_func(s_val).flatten().astype(float)
        N_arr = N_func(s_val).flatten().astype(float)
        B_arr = B_func(s_val).flatten().astype(float)

        # 2D airfoil points in local plane
        x_local = X_airfoil_func(s_val, t_vals_cad)
        y_local = Y_airfoil_func(s_val, t_vals_cad)

        points_3d = []
        for (x_l, y_l) in zip(x_local, y_local):
            pt_3d = C_arr + x_l * N_arr + y_l * B_arr
            points_3d.append(cq.Vector(pt_3d[0], pt_3d[1], pt_3d[2]))

        # close loop
        points_3d.append(points_3d[0])

        # NaN check
        np_points = np.array([(v.x, v.y, v.z) for v in points_3d])
        if np.isnan(np_points).any():
            if verbose:
                print(f"!!! DEBUG: NaN detected in points at s={s_val}. Skipping rib.")
            continue

        # duplicate check (optional debug)
        duplicates = 0
        for i in range(len(points_3d) - 1):
            if (points_3d[i] - points_3d[i+1]).Length < 1e-7:
                duplicates += 1

        closure_dist = (points_3d[0] - points_3d[-1]).Length

        if verbose:
            print(
                f"  s={s_val:8.4f},  pts={len(points_3d):3d},  "
                f"closed={closure_dist:e},  dups={duplicates:3d}"
            )

        try:
            spline_edge = cq.Edge.makeSpline(points_3d)
            wire = cq.Wire.assembleEdges([spline_edge])
            rib_wires.append(wire)
        except Exception as e:
            if verbose:
                print("\n" + "!"*20 + " CRASH DETECTED " + "!"*20)
                print("Offending points (first 5):")
                for i, p in enumerate(points_3d[:5]):
                    print(f" {i}: ({p.x:8.4f}, {p.y:8.4f}, {p.z:8.4f})")
                print("Offending points (last 5):")
                for i, p in enumerate(points_3d[-5:]):
                    print(f" {len(points_3d)-5+i}: ({p.x:8.4f}, {p.y:8.4f}, {p.z:8.4f})")
            raise e

    if not rib_wires:
        raise RuntimeError("No ribs were successfully created. Loft failed.")

    blade_solid = cq.Solid.makeLoft(rib_wires)

    if verbose:
        print(f"Blade solid created in {time.time() - start_time:.2f} seconds.")

    # ---------------- PT 7: Hub + assembly -----------------------

    hub_solid = cq.Solid.makeCylinder(
        hub_radius,
        hub_length,
        pnt=cq.Vector(0, 0, -hub_length / 2.0),
        dir=cq.Vector(0, 0, 1)
    )

    propeller_solid = hub_solid

    for i in range(num_blades):
        angle_deg = i * (360.0 / num_blades)
        blade_copy = blade_solid.rotate(
            (0, 0, 0),
            (0, 0, 1),
            angle_deg
        )
        propeller_solid = propeller_solid.fuse(blade_copy)
        if verbose:
            print(f"Blade {i+1}/{num_blades} fused.")

    if verbose:
        print("Assembly complete.")

    # NOTE: we do NOT apply global_scale here; we’ll handle that in the optimizer.
    return propeller_solid


# =====================================================================================
# 2. CadQuery <-> Trimesh utilities + metrics
# =====================================================================================

def cq_solid_to_trimesh(solid, tmp_path="./temp/tmp_candidate.stl"):
    """
    Export a CadQuery solid to STL and load it as a trimesh.Trimesh.
    """
    cq.exporters.export(solid, tmp_path, "STL")
    mesh = trimesh.load_mesh(tmp_path)
    return mesh


def load_target_mesh(path):
    """
    Load the target propeller as a trimesh mesh.
    Supports STEP (via CadQuery) or STL directly.
    """
    ext = os.path.splitext(path)[1].lower()
    if ext in [".stp", ".step"]:
        # Import with CadQuery, then export to STL -> Trimesh
        wp = cq.importers.importStep(path)
        # get first solid
        target_solid = wp.val()
        tmp_path = "tmp_target.stl"
        cq.exporters.export(target_solid, tmp_path, "STL")
        mesh = trimesh.load_mesh(tmp_path)
        return mesh
    elif ext == ".stl":
        mesh = trimesh.load_mesh(path)
        return mesh
    else:
        raise ValueError(f"Unsupported target file extension: {ext}")


def iou_mesh(mesh_a, mesh_b, pitch=2.0):
    """
    Approximate 3D IoU by sampling a regular grid over the union bounding box.
    pitch: spacing between sample points in world units (bigger = faster, coarser)
    """
    bounds = np.vstack([mesh_a.bounds, mesh_b.bounds])
    mins = bounds.min(axis=0)
    maxs = bounds.max(axis=0)

    extent = maxs - mins
    steps = np.maximum((extent / pitch).astype(int), 1)

    xs = np.linspace(mins[0], maxs[0], steps[0])
    ys = np.linspace(mins[1], maxs[1], steps[1])
    zs = np.linspace(mins[2], maxs[2], steps[2])

    X, Y, Z = np.meshgrid(xs, ys, zs, indexing="ij")
    pts = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])

    inside_a = mesh_a.contains(pts)
    inside_b = mesh_b.contains(pts)

    inter = np.logical_and(inside_a, inside_b).sum()
    union = np.logical_or(inside_a, inside_b).sum()

    if union == 0:
        return 0.0
    return inter / union


def mesh_diameter(mesh):
    """
    Approximate overall propeller diameter as the max bounding-box extent.
    """
    return mesh.bounding_box.extents.max()


def mesh_volume(mesh):
    return mesh.volume


# =====================================================================================
# 3. Parameter vector mapping for optimizer
# =====================================================================================

def vector_to_params(x):
    return {
        "global_scale": float(x[0]),
        "hub_radius":   float(x[1]),
        "hub_length":   float(x[2]),
        "num_blades":   int(x[3]),
        "m":            float(x[4]),
        "p":            float(x[5]),
        "thickness":    float(x[6]),
        "loc_ctrl_point2": [float(x[7]),  float(x[8]),  float(x[9])],
        "loc_ctrl_point3": [float(x[10]), float(x[11]), float(x[12])],
        "blade_vector":    [float(x[13]), float(x[14])],
        "a_AoA": float(x[15]),
        "b_AoA": float(x[16]),
        "c_AoA": float(x[17]),
        "d_AoA": float(x[18]),
        "e_AoA": float(x[19]),
        "a_scX": float(x[20]),
        "b_scX": float(x[21]),
        "c_scX": float(x[22]),
        "d_scX": float(x[23]),
        "e_scX": float(x[24]),
        "a_scY": float(x[25]),
        "b_scY": float(x[26]),
        "c_scY": float(x[27]),
        "d_scY": float(x[28]),
        "e_scY": float(x[29]),
        "apply_thickness_normal": False,
    }


# Bounds (tune as needed; these are reasonable-ish starting points)
bounds = [
    (3.0, 12.0),    # 0: global_scale
    (2.0, 10.0),    # 1: hub_radius
    (5.0, 40.0),    # 2: hub_length
    (3.0, 3.0),     # 3: num_blades

    (0.0, 0.2),     # 4: m
    (0.1, 0.9),     # 5: p
    (0.2, 2.0),     # 6: thickness

    (-20, 20),      # 7:  loc2_x
    (-20, 20),      # 8:  loc2_y
    (0,   60),      # 9:  loc2_z

    (-20, 20),      # 10: loc3_x
    (-20, 20),      # 11: loc3_y
    (0,   60),      # 12: loc3_z

    (-40, 10),      # 13: blade_dx
    (-10, 10),      # 14: blade_dz

    (0, np.pi),  # 15: a_AoA
    (0, np.pi),  # 16: b_AoA
    (0, np.pi),  # 17: c_AoA
    (0, 2*np.pi),  # 18: d_AoA
    (0, 2*np.pi),  # 19: e_AoA

    (0.5, 5.0),    # 20: a_scX
    (0.5, 5.0),    # 21: b_scX
    (0.5, 5.0),    # 22: c_scX
    (0.5, 5.0),    # 23: d_scX
    (0.1, 5.0),     # 24: e_scX

    (0.5, 5.0),    # 25: a_scY
    (0.5, 5.0),    # 26: b_scY
    (0.5, 5.0),    # 27: c_scY
    (0.5, 5.0),    # 28: d_scY
    (0.1, 5.0),     # 29: e_scY
]


# =====================================================================================
# 4. Cost function + optimization driver
# =====================================================================================

# Set this to your real target file path:
TARGET_FILE = "C:\\Users\\varun\\Downloads\\ye et al torprop 12in rotCorrected (1).step"  # or "target_propeller.stl"

print(f"Loading target mesh from: {TARGET_FILE}")
target_mesh = load_target_mesh(TARGET_FILE)
TARGET_DIAM = mesh_diameter(target_mesh)
TARGET_VOL  = mesh_volume(target_mesh)
print(f"Target diameter ~ {TARGET_DIAM:.3f}, volume ~ {TARGET_VOL:.3f}")


def propeller_cost(x,
                   pitch=3.0,
                   w_iou=1.0,
                   w_diam=0.3,
                   w_vol=0.3):
    # Fast cooperative abort
    if _STOP.is_set():
        return 1e9

    # Increment global counter atomically
    with EVAL_COUNTER.get_lock():
        EVAL_COUNTER.value += 1
        eval_id = EVAL_COUNTER.value

    # Lightweight pre‑compute log
    if eval_id % LOG_EVERY == 0:
        elapsed = _fmt_seconds(time.time() - START_TIME)
        print(f"[{elapsed}] [PID {os.getpid()}] Eval {eval_id} START "
              f"x0(thickness)={x[0]:.3f} loc2=({x[1]:.2f},{x[2]:.2f},{x[3]:.2f}) "
              f"loc3=({x[4]:.2f},{x[5]:.2f},{x[6]:.2f}) blade=({x[7]:.2f},{x[8]:.2f}) scale={x[9]:.2f}")

    params = vector_to_params(x)
    global_scale = params["global_scale"]

    try:
        solid = build_propeller(params,
                                s_resolution_cad=15,
                                t_resolution_cad=30,
                                verbose=False)
    except Exception as e:
        print(f"[PID {os.getpid()}] Eval {eval_id} GEOMETRY FAIL: {e}")
        return 1e6

    solid_scaled = solid.scale(global_scale)
    try:
        cand_mesh = cq_solid_to_trimesh(solid_scaled, tmp_path=f"./temp/tmp_candidate_{os.getpid()}.stl")
    except Exception as e:
        print(f"[PID {os.getpid()}] Eval {eval_id} MESH EXPORT FAIL: {e}")
        return 1e6

    iou = iou_mesh(cand_mesh, target_mesh, pitch=pitch)
    diam = mesh_diameter(cand_mesh)
    vol  = mesh_volume(cand_mesh)

    diam_err = abs(diam - TARGET_DIAM) / max(TARGET_DIAM, 1e-6)
    vol_err  = abs(vol  - TARGET_VOL)  / max(TARGET_VOL,  1e-6)

    cost = (w_iou * (1.0 - iou) +
            w_diam * diam_err +
            w_vol  * vol_err)

    if eval_id % LOG_EVERY == 0:
        print(f"[PID {os.getpid()}] Eval {eval_id} DONE "
              f"IoU={iou:.3f} diam_err={diam_err:.3f} vol_err={vol_err:.3f} cost={cost:.3f}")

    if eval_id % SUMMARY_EVERY == 0:
        elapsed = time.time() - START_TIME
        rate = eval_id / max(elapsed, 1e-6)
        print(f"--- SUMMARY: {eval_id} evaluations in {_fmt_seconds(elapsed)} "
              f"({rate:.2f} eval/s) ---")

    return float(cost)


def run_optimization():
    global START_TIME
    START_TIME = time.time()
    print("Starting differential evolution...")
    print(f"Population size = popsize * dim = {10} * {len(bounds)}")
    ctx = mp.get_context("spawn")
    processes_to_spawn = os.cpu_count()
    processes_to_spawn = 5
    pool = ctx.Pool(processes=processes_to_spawn)
    try:
        result = differential_evolution(
            func=propeller_cost,
            bounds=bounds,
            maxiter=40,
            popsize=10,
            tol=1e-3,
            seed=42,
            polish=True,
            workers=pool.map,
            updating="deferred"
        )
        print("Optimization complete.")
    except KeyboardInterrupt:
        print("KeyboardInterrupt received. Stopping workers...")
        _STOP.set()
        pool.terminate()
        pool.join()
        print("All workers terminated.")
        return
    finally:
        if not pool._state == mp.pool.TERMINATE:
            pool.close()
            pool.join()

    print(f"Total evaluations: {EVAL_COUNTER.value}")
    # ...existing code...

    print("\nOptimization result:")
    print(result)

    best_x = result.x
    best_params = vector_to_params(best_x)
    print("\nBest parameter set:")
    for k, v in best_params.items():
        print(f"  {k}: {v}")

    # Build a high-res final geometry with best params
    print("\nRebuilding best geometry at high resolution...")
    solid = build_propeller(best_params,
                            s_resolution_cad=25,
                            t_resolution_cad=60,
                            verbose=True)

    solid_scaled = solid.scale(best_params["global_scale"])

    out_name = "optimized_propeller.step"
    cq.exporters.export(solid_scaled, out_name, "STEP")
    print(f"Exported optimized propeller to: {out_name}")


if __name__ == "__main__":
    run_optimization()
