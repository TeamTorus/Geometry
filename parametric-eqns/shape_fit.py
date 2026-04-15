import os
import numpy as np
import cadquery as cq
import trimesh

from scipy.optimize import differential_evolution, minimize
from scipy.spatial import cKDTree

import importlib


# ============================================================
# 0) Load initial guess from config_toroidal.py
# ============================================================

CFG_MODULE = "config_toroidal"   # config_toroidal.py in same folder
cfg = importlib.import_module(CFG_MODULE)

def params_from_config(cfg):
    return {
        "global_scale": float(cfg.global_scale),
        "hub_radius": float(cfg.hub_radius),
        "hub_length": float(cfg.hub_length),
        "num_blades": int(cfg.num_blades),

        "m": float(cfg.m),
        "p": float(cfg.p),
        "thickness": float(cfg.thickness),

        "loc_ctrl_point2": [float(cfg.loc_ctrl_point2[0]), float(cfg.loc_ctrl_point2[1]), float(cfg.loc_ctrl_point2[2])],
        "loc_ctrl_point3": [float(cfg.loc_ctrl_point3[0]), float(cfg.loc_ctrl_point3[1]), float(cfg.loc_ctrl_point3[2])],
        "blade_vector": [float(cfg.blade_vector[0]), float(cfg.blade_vector[1])],

        "a_AoA": float(cfg.a_AoA), "b_AoA": float(cfg.b_AoA), "c_AoA": float(cfg.c_AoA), "d_AoA": float(cfg.d_AoA), "e_AoA": float(cfg.e_AoA),
        "a_scX": float(cfg.a_scX), "b_scX": float(cfg.b_scX), "c_scX": float(cfg.c_scX), "d_scX": float(cfg.d_scX), "e_scX": float(cfg.e_scX),
        "a_scY": float(cfg.a_scY), "b_scY": float(cfg.b_scY), "c_scY": float(cfg.c_scY), "d_scY": float(cfg.d_scY), "e_scY": float(cfg.e_scY),

        "apply_thickness_normal": int(bool(cfg.apply_thickness_normal)),
    }

x0_params = params_from_config(cfg)


# ============================================================
# 1) Target loading (once)
# ============================================================

TARGET_FILE = "tmp_target.step"  # or .stl
TARGET_SAMPLES = 60000

CAND_SAMPLES_S = 60
CAND_SAMPLES_T = 160

# Chamfer + regularization weights
W_CHAMFER = 1.0
W_DIAM    = 0.00

# Exploit penalties (tune these)
W_SCALE_COLLAPSE = 0.40     # penalize scx/scy going tiny
W_SCALE_WILD      = 0.20     # penalize scx/scy huge oscillations
W_FRAME_BAD       = 2.00     # penalize frame singularities / NaNs
W_DELTA_REG       = 0.05     # keep params near initial guess

# Robust chamfer trimming: ignore worst X% distances (reduces "spray points" exploits)
TRIM_FRACTION = 0.10   # 0.0 = off, 0.10 means drop worst 10%

CHUNK = 200000


def load_target_mesh(path: str) -> trimesh.Trimesh:
    ext = os.path.splitext(path)[1].lower()

    if ext == ".stl":
        mesh = trimesh.load_mesh(path)
        if isinstance(mesh, trimesh.Scene):
            mesh = trimesh.util.concatenate(tuple(mesh.geometry.values()))
        return mesh

    if ext in [".step", ".stp"]:
        shape = cq.importers.importStep(path).val()
        verts, tris = shape.tessellate(0.5)
        verts = np.array([[v.x, v.y, v.z] for v in verts], dtype=np.float64)
        tris  = np.asarray(tris, dtype=np.int64)
        return trimesh.Trimesh(vertices=verts, faces=tris, process=True)

    raise ValueError(f"Unsupported file extension: {ext}")


print(f"Loading target mesh: {TARGET_FILE}")
target_mesh = load_target_mesh(TARGET_FILE)
target_mesh.remove_unreferenced_vertices()
target_mesh.process(validate=True)

target_points, _ = trimesh.sample.sample_surface(target_mesh, TARGET_SAMPLES)
target_points = target_points.astype(np.float64)

TARGET_BOUNDS = target_mesh.bounds
TARGET_DIAM = (TARGET_BOUNDS[1] - TARGET_BOUNDS[0]).max()
print(f"Target diameter (bbox max extent): {TARGET_DIAM:.4f}")


# ============================================================
# 2) Geometry helpers (same as yours)
# ============================================================

def rational_cubic_bezier(P, w, s_vals):
    s = s_vals[:, None]
    one = 1.0 - s
    B0 = one**3
    B1 = 3*s*one**2
    B2 = 3*s**2*one
    B3 = s**3
    W0 = B0 * w[0]
    W1 = B1 * w[1]
    W2 = B2 * w[2]
    W3 = B3 * w[3]
    denom = (W0 + W1 + W2 + W3)
    num = W0 * P[0] + W1 * P[1] + W2 * P[2] + W3 * P[3]
    return num / (denom + 1e-12)

def frenet_frames(curve_pts):
    dC = np.gradient(curve_pts, axis=0)
    T = dC / (np.linalg.norm(dC, axis=1, keepdims=True) + 1e-12)

    dT = np.gradient(T, axis=0)
    nrm = np.linalg.norm(dT, axis=1, keepdims=True)
    bad = (nrm[:,0] < 1e-8)
    N = dT / (nrm + 1e-12)

    B = np.cross(T, N)
    B = B / (np.linalg.norm(B, axis=1, keepdims=True) + 1e-12)

    N = np.cross(B, T)
    N = N / (np.linalg.norm(N, axis=1, keepdims=True) + 1e-12)

    return T, N, B, bad

def naca_4digit_xy(t_vals, m, p, thickness, apply_thickness_normal=False):
    t_vals = np.asarray(t_vals, dtype=np.float64)
    upper = t_vals < 1.0
    x = np.where(upper, t_vals, 2.0 - t_vals)

    yt = 5.0 * thickness * (
        0.2969*np.sqrt(np.clip(x, 0, 1))
        - 0.1260*x
        - 0.3516*x**2
        + 0.2843*x**3
        - 0.1036*x**4
    )

    yc = np.where(
        x <= p,
        (m / (p**2 + 1e-12)) * (2*p*x - x**2),
        (m / ((1-p)**2 + 1e-12)) * ((1 - 2*p) + 2*p*x - x**2)
    )

    if apply_thickness_normal:
        dyc_dx = np.where(
            x <= p,
            2*m/(p**2 + 1e-12) * (p - x),
            2*m/(((1-p)**2 + 1e-12)) * (p - x)
        )
        theta = np.arctan(dyc_dx)
    else:
        theta = 0.0

    yu = yc + yt*np.cos(theta)
    yl = yc - yt*np.cos(theta)
    xu = x - yt*np.sin(theta)
    xl = x + yt*np.sin(theta)

    y = np.where(upper, yu, yl)
    x_out = np.where(upper, xu, xl) - 0.5
    return x_out, y

def clamp_scaling(s_vals):
    s_vals = np.asarray(s_vals, dtype=np.float64)
    return -1.0 + 1.0/(1.0 + np.exp(-200.0*(s_vals - 0.01))) + 1.0/(1.0 + np.exp(-200.0*(0.99 - s_vals)))

def poly4(s_vals, a, b, c, d, e):
    return a*s_vals**4 + b*s_vals**3 + c*s_vals**2 + d*s_vals + e


# ============================================================
# 3) Candidate point cloud generation (yours + diagnostics)
# ============================================================

def candidate_pointcloud_and_diagnostics(params, n_s=CAND_SAMPLES_S, n_t=CAND_SAMPLES_T):
    """
    Returns (points, diagnostics dict)
    diagnostics includes:
      - bad_frame_frac
      - min_sc_abs
      - max_sc_abs
      - sc_smoothness (2nd diff energy)
    """
    global_scale = float(params["global_scale"])
    hub_radius   = float(params["hub_radius"])
    hub_length   = float(params["hub_length"])
    num_blades   = max(1, int(round(params["num_blades"])))

    m = float(params["m"])
    p = float(params["p"])
    thickness = float(params["thickness"])

    loc2 = params["loc_ctrl_point2"]
    loc3 = params["loc_ctrl_point3"]
    bv   = params["blade_vector"]

    a_AoA, b_AoA, c_AoA, d_AoA, e_AoA = (params["a_AoA"], params["b_AoA"], params["c_AoA"], params["d_AoA"], params["e_AoA"])
    a_scX, b_scX, c_scX, d_scX, e_scX = (params["a_scX"], params["b_scX"], params["c_scX"], params["d_scX"], params["e_scX"])
    a_scY, b_scY, c_scY, d_scY, e_scY = (params["a_scY"], params["b_scY"], params["c_scY"], params["d_scY"], params["e_scY"])
    apply_thickness_normal = bool(round(params["apply_thickness_normal"]))

    inset_ratio = 4.0/8.0
    blade_hub_radius = inset_ratio * hub_radius

    s_vals = np.linspace(0.0, 1.0, n_s, endpoint=True)
    t_vals = np.linspace(0.0, 2.0, n_t, endpoint=False)

    x2d, y2d = naca_4digit_xy(t_vals, m=m, p=p, thickness=thickness,
                             apply_thickness_normal=apply_thickness_normal)

    AoA = poly4(s_vals, a_AoA, b_AoA, c_AoA, d_AoA, e_AoA)
    scx = poly4(s_vals, a_scX, b_scX, c_scX, d_scX, e_scX)
    scy = poly4(s_vals, a_scY, b_scY, c_scY, d_scY, e_scY)

    clamp = clamp_scaling(s_vals)
    scx = scx * clamp
    scy = scy * clamp

    # diagnostics
    sc_abs = np.sqrt(scx**2 + scy**2)
    min_sc = float(np.min(sc_abs))
    max_sc = float(np.max(sc_abs))
    # smoothness via second difference energy
    dd = np.diff(sc_abs, n=2)
    smooth = float(np.mean(dd*dd)) if len(dd) > 0 else 0.0

    cosA = np.cos(AoA)[:, None]
    sinA = np.sin(AoA)[:, None]
    Xr = (cosA * x2d[None, :] - sinA * y2d[None, :]) * scx[:, None]
    Yr = (sinA * x2d[None, :] + cosA * y2d[None, :]) * scy[:, None]

    # centerline control points (relative z0 fix)
    z0 = hub_length / 2.0 - 1.0
    ctrl_point1 = np.array([blade_hub_radius, 0.0, z0], dtype=np.float64)

    theta4 = bv[0] / (blade_hub_radius + 1e-12)
    ctrl_point4 = np.array([blade_hub_radius*np.cos(theta4), blade_hub_radius*np.sin(theta4), z0 - bv[1]], dtype=np.float64)

    theta2 = loc2[0] / (blade_hub_radius + 1e-12)
    r2 = blade_hub_radius + loc2[2]
    ctrl_point2 = np.array([r2*np.cos(theta2), r2*np.sin(theta2), z0 - loc2[1]], dtype=np.float64)

    theta3 = loc3[0] / (blade_hub_radius + 1e-12)
    r3 = blade_hub_radius + loc3[2]
    ctrl_point3 = np.array([r3*np.cos(theta3), r3*np.sin(theta3), z0 - loc3[1]], dtype=np.float64)

    P = np.stack([ctrl_point1, ctrl_point2, ctrl_point3, ctrl_point4], axis=0)
    w = np.array([1.0, 1.0, 1.0, 1.0], dtype=np.float64)

    C = rational_cubic_bezier(P, w, s_vals)
    _, N, B, bad = frenet_frames(C)
    bad_frac = float(np.mean(bad))

    blade_pts = (C[:, None, :] + Xr[:, :, None] * N[:, None, :] + Yr[:, :, None] * B[:, None, :]).reshape(-1, 3)

    all_blades = []
    for i in range(num_blades):
        ang = (2.0*np.pi*i)/num_blades
        ca, sa = np.cos(ang), np.sin(ang)
        R = np.array([[ca, -sa, 0.0],
                      [sa,  ca, 0.0],
                      [0.0, 0.0, 1.0]], dtype=np.float64)
        all_blades.append(blade_pts @ R.T)
    blades_pts = np.vstack(all_blades)

    # hub sample
    n_theta = 200
    n_z = 80
    th = np.linspace(0, 2*np.pi, n_theta, endpoint=False)
    zz = np.linspace(-hub_length/2.0, hub_length/2.0, n_z)
    TH, ZZ = np.meshgrid(th, zz, indexing="xy")
    hub_pts = np.column_stack([hub_radius*np.cos(TH).ravel(), hub_radius*np.sin(TH).ravel(), ZZ.ravel()])

    pts = np.vstack([blades_pts, hub_pts]) * global_scale

    diag = dict(
        bad_frame_frac=bad_frac,
        min_sc_abs=min_sc,
        max_sc_abs=max_sc,
        sc_smoothness=smooth,
    )
    return pts, diag

# 3.5 ICP
def align_icp(candidate_pts, target_pts):
    """
    Align candidate point cloud to target using rigid ICP.
    Returns transformed candidate points.
    """
    cand = trimesh.points.PointCloud(candidate_pts)
    targ = trimesh.points.PointCloud(target_pts)

    matrix, _, _ = trimesh.registration.icp(
        cand.vertices,
        targ.vertices,
        scale=False,
        max_iterations=30
    )

    # apply transform
    pts = trimesh.transform_points(candidate_pts, matrix)
    return pts

# ============================================================
# 4) Robust Chamfer distance (trimmed)
# ============================================================

def chamfer_distance_trimmed(A, B, trim_fraction=0.0):
    treeB = cKDTree(B)
    dA_parts = []
    for i in range(0, len(A), CHUNK):
        d, _ = treeB.query(A[i:i+CHUNK], k=1, workers=-1)
        dA_parts.append(d*d)
    dA = np.concatenate(dA_parts)

    treeA = cKDTree(A)
    dB_parts = []
    for i in range(0, len(B), CHUNK):
        d, _ = treeA.query(B[i:i+CHUNK], k=1, workers=-1)
        dB_parts.append(d*d)
    dB = np.concatenate(dB_parts)

    if trim_fraction > 0.0:
        keepA = int(max(1, (1.0 - trim_fraction) * len(dA)))
        keepB = int(max(1, (1.0 - trim_fraction) * len(dB)))
        dA = np.partition(dA, keepA-1)[:keepA]
        dB = np.partition(dB, keepB-1)[:keepB]

    return float(dA.mean() + dB.mean())


def bbox_diameter(pts):
    bmin = pts.min(axis=0)
    bmax = pts.max(axis=0)
    return float((bmax - bmin).max())


# ============================================================
# 5) Delta-parameterization around x0
# ============================================================

# Vector layout (same as yours)
def params_to_vector(p):
    return np.array([
        p["global_scale"],
        p["hub_radius"],
        p["hub_length"],
        float(p["num_blades"]),   # stored as float in optimizer

        p["m"], p["p"], p["thickness"],

        p["loc_ctrl_point2"][0], p["loc_ctrl_point2"][1], p["loc_ctrl_point2"][2],
        p["loc_ctrl_point3"][0], p["loc_ctrl_point3"][1], p["loc_ctrl_point3"][2],
        p["blade_vector"][0], p["blade_vector"][1],

        p["a_AoA"], p["b_AoA"], p["c_AoA"], p["d_AoA"], p["e_AoA"],
        p["a_scX"], p["b_scX"], p["c_scX"], p["d_scX"], p["e_scX"],
        p["a_scY"], p["b_scY"], p["c_scY"], p["d_scY"], p["e_scY"],

        float(p["apply_thickness_normal"]),
    ], dtype=np.float64)

def vector_to_params(x):
    return {
        "global_scale": float(x[0]),
        "hub_radius": float(x[1]),
        "hub_length": float(x[2]),
        "num_blades": int(round(x[3])),

        "m": float(x[4]),
        "p": float(x[5]),
        "thickness": float(x[6]),

        "loc_ctrl_point2": [float(x[7]), float(x[8]), float(x[9])],
        "loc_ctrl_point3": [float(x[10]), float(x[11]), float(x[12])],
        "blade_vector": [float(x[13]), float(x[14])],

        "a_AoA": float(x[15]), "b_AoA": float(x[16]), "c_AoA": float(x[17]), "d_AoA": float(x[18]), "e_AoA": float(x[19]),
        "a_scX": float(x[20]), "b_scX": float(x[21]), "c_scX": float(x[22]), "d_scX": float(x[23]), "e_scX": float(x[24]),
        "a_scY": float(x[25]), "b_scY": float(x[26]), "c_scY": float(x[27]), "d_scY": float(x[28]), "e_scY": float(x[29]),

        "apply_thickness_normal": int(round(x[30])),
    }

x0 = params_to_vector(x0_params)

# Per-parameter "fine tune" step sizes (THIS is the main lever)
# Units match your vector entries.
STEP = np.array([
    0.05,   # global_scale
    0.05,   # hub_radius
    0.10,   # hub_length
    0.0,    # num_blades locked (keep 3)  <-- IMPORTANT

    0.005,   # m
    0.005,   # p
    0.005,   # thickness

    0.01,    # loc2_circ
    0.01,    # loc2_z
    0.01,    # loc2_rad

    0.01,    # loc3_circ
    0.01,    # loc3_z
    0.01,    # loc3_rad

    0.01,    # blade_circ
    0.01,    # blade_z

    0.025,   # a_AoA
    0.025,   # b_AoA
    0.025,   # c_AoA
    0.025,   # d_AoA
    0.025,   # e_AoA

    0.025,   # a_scX
    0.025,   # b_scX
    0.025,   # c_scX
    0.025,   # d_scX
    0.025,   # e_scX

    0.025,   # a_scY
    0.025,   # b_scY
    0.025,   # c_scY
    0.025,   # d_scY
    0.025,   # e_scY

    0.0,    # apply_thickness_normal locked
], dtype=np.float64)

# Delta bounds: each d_i in [-1, +1] means +/- STEP around x0
DELTA_BOUNDS = [(-1.0, 1.0) for _ in range(len(x0))]


def make_x_from_delta(d):
    d = np.asarray(d, dtype=np.float64)
    x = x0 + STEP * d
    return x


# ============================================================
# 6) Cost function: chamfer + diameter + anti-exploit penalties
# ============================================================

def cost_from_delta(d):
    x = make_x_from_delta(d)
    params = vector_to_params(x)

    # keep parameters sane near initial guess
    if not (0.0 <= params["m"] <= 0.20): return 1e9
    if not (0.05 <= params["p"] <= 0.95): return 1e9
    if params["hub_radius"] <= 0.2 or params["hub_length"] <= 0.5: return 1e9
    if params["thickness"] <= 0.05: return 1e9

    try:
        cand_pts, diag = candidate_pointcloud_and_diagnostics(params, n_s=CAND_SAMPLES_S, n_t=CAND_SAMPLES_T)
    except Exception:
        return 1e9

    if not np.isfinite(cand_pts).all():
        return 1e9
    
    # optional normalization
    def center_points(P):
        return P - P.mean(axis=0)

    cand_pts = center_points(cand_pts)
    target_pts = center_points(target_points)
    
    # downsample for ICP speed
    idx = np.random.choice(len(cand_pts), min(8000, len(cand_pts)), replace=False)
    cand_sample = cand_pts[idx]

    cand_sample = align_icp(cand_sample, target_points)

    # compute transform on full cloud
    matrix, _, _ = trimesh.registration.icp(
        cand_sample,
        target_points,
        scale=False,
        max_iterations=20
    )

    cand_pts = trimesh.transform_points(cand_pts, matrix)

    # main term
    cham = chamfer_distance_trimmed(cand_pts, target_points, trim_fraction=TRIM_FRACTION)

    # diameter term
    cand_d = bbox_diameter(cand_pts)
    diam_err = abs(cand_d - TARGET_DIAM) / max(TARGET_DIAM, 1e-9)

    # anti-exploit penalties
    # 1) frame instability penalty
    frame_pen = diag["bad_frame_frac"]

    # 2) collapse penalty: discourage sc going near zero
    # (blade vanishing is a common exploit)
    collapse_pen = max(0.0, (0.25 - diag["min_sc_abs"]))  # hinge at 0.25

    # 3) wild scale penalty: discourage huge scales / oscillations
    wild_pen = max(0.0, (diag["max_sc_abs"] - 20.0)) + 5.0 * diag["sc_smoothness"]

    # 4) regularize delta size (stay near x0)
    d = np.asarray(d, dtype=np.float64)
    delta_reg = float(np.mean(d*d))

    cost = (
        W_CHAMFER * cham +
        W_DIAM * diam_err +
        W_FRAME_BAD * frame_pen +
        W_SCALE_COLLAPSE * collapse_pen +
        W_SCALE_WILD * wild_pen +
        W_DELTA_REG * delta_reg
    )

    print(f"cost={cost:.4e} cham={cham:.4e} diam={diam_err:.3f} bad={frame_pen:.3f} minsc={diag['min_sc_abs']:.3f} dreg={delta_reg:.3f}")
    return float(cost)


# ============================================================
# 7) Optimize (local-ish) and report best
# ============================================================

def run():
    # Stage 1: Differential evolution in a *small local box* around x0
    # This helps avoid poor local minima but can't wander into crazy regions.
    print("\n=== Stage 1: local differential evolution around config ===")
    de = differential_evolution(
        cost_from_delta,
        bounds=DELTA_BOUNDS,
        maxiter=35,
        popsize=10,
        tol=5e-3,
        seed=42,
        workers=-1,
        updating="deferred",
        polish=False,  # we'll do our own polish step below
    )

    d1 = de.x
    print("\nDE done:", de.fun)

    # Stage 2: local polish (Powell) starting from DE result
    print("\n=== Stage 2: Powell local polish ===")
    res = minimize(
        cost_from_delta,
        d1,
        method="Powell",
        options=dict(maxiter=120, xtol=1e-3, ftol=1e-3, disp=True),
    )

    d_best = res.x
    x_best = make_x_from_delta(d_best)
    best_params = vector_to_params(x_best)

    # enforce locked discrete params
    best_params["num_blades"] = int(x0_params["num_blades"])
    best_params["apply_thickness_normal"] = int(x0_params["apply_thickness_normal"])

    print("\n=== Best params (fine-tuned around config) ===")
    for k, v in best_params.items():
        print(f"{k}: {v}")

    print("\nBest delta (scaled):")
    for i, dv in enumerate(d_best):
        if abs(dv) > 1e-3:
            print(f"  d[{i}] = {dv:.4f} (step {STEP[i]:.4f})")

    # Save best into a python file you can reuse
    save_best_config(best_params, out_path="config_toroidal_best.py")

    return best_params


def save_best_config(p, out_path="config_toroidal_best.py"):
    with open(out_path, "w") as f:
        f.write("import numpy as np\n\n")
        f.write("# Auto-generated fine-tuned parameters\n\n")
        f.write("s_domain = [0, 1]\n")
        f.write("t_domain = [0, 2]\n\n")
        f.write(f"s_resolution_cad = {getattr(cfg,'s_resolution_cad',25)}\n")
        f.write(f"t_resolution_cad = {getattr(cfg,'t_resolution_cad',60)}\n\n")
        f.write(f"global_scale = {p['global_scale']}\n")
        f.write(f"hub_radius = {p['hub_radius']}\n")
        f.write(f"hub_length = {p['hub_length']}\n")
        f.write(f"num_blades = {p['num_blades']}\n\n")
        f.write(f"m = {p['m']}\n")
        f.write(f"p = {p['p']}\n")
        f.write(f"thickness = {p['thickness']}\n\n")
        f.write(f"loc_ctrl_point2 = {p['loc_ctrl_point2']}\n")
        f.write(f"loc_ctrl_point3 = {p['loc_ctrl_point3']}\n")
        f.write(f"blade_vector = {p['blade_vector']}\n\n")
        f.write(f"a_AoA = {p['a_AoA']}\n")
        f.write(f"b_AoA = {p['b_AoA']}\n")
        f.write(f"c_AoA = {p['c_AoA']}\n")
        f.write(f"d_AoA = {p['d_AoA']}\n")
        f.write(f"e_AoA = {p['e_AoA']}\n\n")
        f.write(f"a_scX = {p['a_scX']}\n")
        f.write(f"b_scX = {p['b_scX']}\n")
        f.write(f"c_scX = {p['c_scX']}\n")
        f.write(f"d_scX = {p['d_scX']}\n")
        f.write(f"e_scX = {p['e_scX']}\n\n")
        f.write(f"a_scY = {p['a_scY']}\n")
        f.write(f"b_scY = {p['b_scY']}\n")
        f.write(f"c_scY = {p['c_scY']}\n")
        f.write(f"d_scY = {p['d_scY']}\n")
        f.write(f"e_scY = {p['e_scY']}\n\n")
        f.write(f"apply_thickness_normal = {bool(p['apply_thickness_normal'])}\n")

    print("Wrote best config to:", out_path)


if __name__ == "__main__":
    run()