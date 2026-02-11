import os
import numpy as np
import cadquery as cq
import trimesh

from scipy.optimize import differential_evolution
from scipy.spatial import cKDTree


# ============================================================
# 0) Target loading (once)
# ============================================================

TARGET_FILE = "tmp_target.step"  # or .stl
TARGET_SAMPLES = 60000                 # surface samples for target point cloud
CAND_SAMPLES_S = 60                    # curve stations (s)
CAND_SAMPLES_T = 160                   # airfoil samples (t) per station

# Chamfer weights / penalties
W_CHAMFER = 1.0
W_DIAM    = 0.15
W_VOL     = 0.00  # volume from pointcloud is rough; leave 0 unless you add better approx

# KD-tree batch size (memory/perf knob)
CHUNK = 200000


def load_target_mesh(path: str) -> trimesh.Trimesh:
    ext = os.path.splitext(path)[1].lower()

    if ext == ".stl":
        mesh = trimesh.load_mesh(path)
        if isinstance(mesh, trimesh.Scene):
            mesh = trimesh.util.concatenate(tuple(mesh.geometry.values()))
        return mesh

    if ext in [".step", ".stp"]:
        # Load STEP once with CadQuery, tessellate in memory (no STL export)
        shape = cq.importers.importStep(path).val()

        # CadQuery tessellate returns (vertices, triangles)
        # tolerance is a trade-off: smaller = more triangles
        verts, tris = shape.tessellate(0.5)

        # Convert CadQuery Vectors → numeric array
        verts = np.array([[v.x, v.y, v.z] for v in verts], dtype=np.float64)
        tris  = np.asarray(tris, dtype=np.int64)


        mesh = trimesh.Trimesh(vertices=verts, faces=tris, process=True)
        return mesh

    raise ValueError(f"Unsupported file extension: {ext}")


print(f"Loading target mesh: {TARGET_FILE}")
target_mesh = load_target_mesh(TARGET_FILE)
target_mesh.remove_unreferenced_vertices()
target_mesh.process(validate=True)

# target surface point cloud (once)
target_points, _ = trimesh.sample.sample_surface(target_mesh, TARGET_SAMPLES)
target_points = target_points.astype(np.float64)

TARGET_BOUNDS = target_mesh.bounds
TARGET_DIAM = (TARGET_BOUNDS[1] - TARGET_BOUNDS[0]).max()
print(f"Target diameter (bbox max extent): {TARGET_DIAM:.4f}")


# ============================================================
# 1) Geometry: numeric NURBS (rational cubic Bezier) + frames
# ============================================================

def rational_cubic_bezier(P, w, s_vals):
    """
    P: (4,3) control points
    w: (4,) weights
    s_vals: (N,) in [0,1]
    Returns: (N,3)
    """
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
    return num / denom


def frenet_frames(curve_pts):
    """
    curve_pts: (N,3)
    Returns T,N,B each (N,3) using numeric derivatives.
    """
    # tangent via gradient
    dC = np.gradient(curve_pts, axis=0)
    T = dC / (np.linalg.norm(dC, axis=1, keepdims=True) + 1e-12)

    dT = np.gradient(T, axis=0)
    N = dT / (np.linalg.norm(dT, axis=1, keepdims=True) + 1e-12)

    B = np.cross(T, N)
    B = B / (np.linalg.norm(B, axis=1, keepdims=True) + 1e-12)

    # Re-orthogonalize N (sometimes numeric drift)
    N = np.cross(B, T)
    N = N / (np.linalg.norm(N, axis=1, keepdims=True) + 1e-12)

    return T, N, B


# ============================================================
# 2) Airfoil: numeric NACA camber + thickness (vectorized)
# ============================================================

def naca_4digit_xy(t_vals, m, p, thickness, apply_thickness_normal=False):
    """
    Replicates your t in [0,2) parameterization:
      t in [0,1): upper surface with x=t
      t in [1,2): lower surface with x=2-t (reversed)
    Then shifts x by -0.5 like your code.
    """

    t_vals = np.asarray(t_vals, dtype=np.float64)
    upper = t_vals < 1.0
    x = np.where(upper, t_vals, 2.0 - t_vals)  # x in [0,1]

    # thickness distribution
    yt = 5.0 * thickness * (
        0.2969*np.sqrt(np.clip(x, 0, 1))
        - 0.1260*x
        - 0.3516*x**2
        + 0.2843*x**3
        - 0.1036*x**4
    )

    # camber line
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

    # upper/lower
    yu = yc + yt*np.cos(theta)
    yl = yc - yt*np.cos(theta)

    xu = x - yt*np.sin(theta)
    xl = x + yt*np.sin(theta)

    y = np.where(upper, yu, yl)
    x_out = np.where(upper, xu, xl)

    # shift x by -0.5 like your sympy code
    x_out = x_out - 0.5
    return x_out, y


def clamp_scaling(s_vals):
    # same logistic clamp you used, numeric
    s_vals = np.asarray(s_vals, dtype=np.float64)
    return -1.0 + 1.0/(1.0 + np.exp(-200.0*(s_vals - 0.01))) + 1.0/(1.0 + np.exp(-200.0*(0.99 - s_vals)))


def poly4(s_vals, a, b, c, d, e):
    s = s_vals
    return a*s**4 + b*s**3 + c*s**2 + d*s + e


# ============================================================
# 3) Candidate point cloud generation (FAST)
# ============================================================

def candidate_pointcloud(params,
                         n_s=CAND_SAMPLES_S,
                         n_t=CAND_SAMPLES_T):
    """
    Generates a surface point cloud (blade surfaces + hub surface)
    without building any B-Rep solids.

    IMPORTANT FIX:
      Control point Z offsets are applied RELATIVE to ctrl_point1.z.
    """
    # --- unpack ---
    global_scale = params["global_scale"]
    hub_radius   = params["hub_radius"]
    hub_length   = params["hub_length"]
    num_blades   = max(1, int(round(params["num_blades"])))

    m = params["m"]
    p = params["p"]
    thickness = params["thickness"]

    loc2 = params["loc_ctrl_point2"]  # [circ, z_offset, radial]
    loc3 = params["loc_ctrl_point3"]
    bv   = params["blade_vector"]     # [circ, z_offset]

    a_AoA, b_AoA, c_AoA, d_AoA, e_AoA = params["a_AoA"], params["b_AoA"], params["c_AoA"], params["d_AoA"], params["e_AoA"]
    a_scX, b_scX, c_scX, d_scX, e_scX = params["a_scX"], params["b_scX"], params["c_scX"], params["d_scX"], params["e_scX"]
    a_scY, b_scY, c_scY, d_scY, e_scY = params["a_scY"], params["b_scY"], params["c_scY"], params["d_scY"], params["e_scY"]

    apply_thickness_normal = bool(round(params["apply_thickness_normal"]))

    inset_ratio = 4.0/8.0
    blade_hub_radius = inset_ratio * hub_radius

    # --- sample parameters ---
    s_vals = np.linspace(0.0, 1.0, n_s, endpoint=True)
    t_vals = np.linspace(0.0, 2.0, n_t, endpoint=False)

    # --- airfoil in local 2D ---
    x2d, y2d = naca_4digit_xy(t_vals, m=m, p=p, thickness=thickness,
                             apply_thickness_normal=apply_thickness_normal)

    # --- along-s scaling + AoA ---
    AoA = poly4(s_vals, a_AoA, b_AoA, c_AoA, d_AoA, e_AoA)
    scx = poly4(s_vals, a_scX, b_scX, c_scX, d_scX, e_scX)
    scy = poly4(s_vals, a_scY, b_scY, c_scY, d_scY, e_scY)

    clamp = clamp_scaling(s_vals)
    scx = scx * clamp
    scy = scy * clamp

    # rotate airfoil per s: we’ll do this by broadcasting
    cosA = np.cos(AoA)[:, None]
    sinA = np.sin(AoA)[:, None]

    # shape: (n_s, n_t)
    Xr = (cosA * x2d[None, :] - sinA * y2d[None, :]) * scx[:, None]
    Yr = (sinA * x2d[None, :] + cosA * y2d[None, :]) * scy[:, None]

    # --- centerline control points (CONFIRMED FIX FOR Z OFFSETS) ---
    z0 = hub_length / 2.0 - 1.0  # ctrl_point1 z baseline

    ctrl_point1 = np.array([blade_hub_radius, 0.0, z0], dtype=np.float64)

    # arc-length -> angle: theta = arc / r
    theta4 = bv[0] / (blade_hub_radius + 1e-12)
    ctrl_point4 = np.array([
        blade_hub_radius * np.cos(theta4),
        blade_hub_radius * np.sin(theta4),
        z0 - bv[1]  # <-- FIX: relative to z0
    ], dtype=np.float64)

    theta2 = loc2[0] / (blade_hub_radius + 1e-12)
    r2 = blade_hub_radius + loc2[2]
    ctrl_point2 = np.array([
        r2 * np.cos(theta2),
        r2 * np.sin(theta2),
        z0 - loc2[1]  # <-- FIX: relative to z0
    ], dtype=np.float64)

    theta3 = loc3[0] / (blade_hub_radius + 1e-12)
    r3 = blade_hub_radius + loc3[2]
    ctrl_point3 = np.array([
        r3 * np.cos(theta3),
        r3 * np.sin(theta3),
        z0 - loc3[1]  # <-- FIX: relative to z0
    ], dtype=np.float64)

    P = np.stack([ctrl_point1, ctrl_point2, ctrl_point3, ctrl_point4], axis=0)
    w = np.array([1.0, 1.0, 1.0, 1.0], dtype=np.float64)

    # --- curve points + frames ---
    C = rational_cubic_bezier(P, w, s_vals)  # (n_s,3)
    _, N, B = frenet_frames(C)               # each (n_s,3)

    # --- surface points for one blade ---
    # C[:,None,:] + Xr[:,:,None]*N[:,None,:] + Yr[:,:,None]*B[:,None,:]
    blade_pts = (C[:, None, :]
                 + Xr[:, :, None] * N[:, None, :]
                 + Yr[:, :, None] * B[:, None, :])

    blade_pts = blade_pts.reshape(-1, 3)

    # --- replicate blades by rotation around Z ---
    all_blades = []
    for i in range(num_blades):
        ang = (2.0*np.pi*i)/num_blades
        ca, sa = np.cos(ang), np.sin(ang)
        R = np.array([[ca, -sa, 0.0],
                      [sa,  ca, 0.0],
                      [0.0, 0.0, 1.0]], dtype=np.float64)
        all_blades.append(blade_pts @ R.T)

    blades_pts = np.vstack(all_blades)

    # --- hub surface sample (cylinder) ---
    # Sample cylinder lateral surface (enough to guide optimizer)
    n_theta = 200
    n_z = 80
    th = np.linspace(0, 2*np.pi, n_theta, endpoint=False)
    zz = np.linspace(-hub_length/2.0, hub_length/2.0, n_z)
    TH, ZZ = np.meshgrid(th, zz, indexing="xy")
    hub_pts = np.column_stack([
        hub_radius * np.cos(TH).ravel(),
        hub_radius * np.sin(TH).ravel(),
        ZZ.ravel()
    ])

    # --- combine + global scale ---
    pts = np.vstack([blades_pts, hub_pts])
    pts = pts * global_scale
    return pts


# ============================================================
# 4) Distance metrics (fast KD-tree Chamfer)
# ============================================================

def chamfer_distance(A, B):
    """
    Symmetric Chamfer distance:
      mean_{a in A} min_{b in B} ||a-b||^2  +  mean_{b in B} min_{a in A} ||b-a||^2
    Uses KD trees.
    """
    treeB = cKDTree(B)
    # chunk to avoid huge memory spikes
    distsA = []
    for i in range(0, len(A), CHUNK):
        chunk = A[i:i+CHUNK]
        d, _ = treeB.query(chunk, k=1, workers=-1)
        distsA.append(d*d)
    dA = np.concatenate(distsA).mean()

    treeA = cKDTree(A)
    distsB = []
    for i in range(0, len(B), CHUNK):
        chunk = B[i:i+CHUNK]
        d, _ = treeA.query(chunk, k=1, workers=-1)
        distsB.append(d*d)
    dB = np.concatenate(distsB).mean()

    return dA + dB


def bbox_diameter(pts):
    bmin = pts.min(axis=0)
    bmax = pts.max(axis=0)
    return (bmax - bmin).max()


# ============================================================
# 5) Optimization variables: ALL tunable (like your “all params”)
# ============================================================

# vector layout
# 0 global_scale
# 1 hub_radius
# 2 hub_length
# 3 num_blades
# 4 m
# 5 p
# 6 thickness
# 7 loc2_circ, 8 loc2_z, 9 loc2_radial
# 10 loc3_circ, 11 loc3_z, 12 loc3_radial
# 13 blade_circ, 14 blade_z
# 15..19 AoA coeffs a..e
# 20..24 scX coeffs a..e
# 25..29 scY coeffs a..e
# 30 apply_thickness_normal

def vector_to_params(x):
    return {
        "global_scale": x[0],
        "hub_radius": x[1],
        "hub_length": x[2],
        "num_blades": x[3],

        "m": x[4],
        "p": x[5],
        "thickness": x[6],

        "loc_ctrl_point2": [x[7], x[8], x[9]],
        "loc_ctrl_point3": [x[10], x[11], x[12]],
        "blade_vector": [x[13], x[14]],

        "a_AoA": x[15], "b_AoA": x[16], "c_AoA": x[17], "d_AoA": x[18], "e_AoA": x[19],
        "a_scX": x[20], "b_scX": x[21], "c_scX": x[22], "d_scX": x[23], "e_scX": x[24],
        "a_scY": x[25], "b_scY": x[26], "c_scY": x[27], "d_scY": x[28], "e_scY": x[29],

        "apply_thickness_normal": x[30],
    }


# Bounds: these matter a LOT for convergence & “reasonable” shapes
bounds = [
    (3.0, 24.0),      # global_scale
    (2.0, 10.0),      # hub_radius
    (5.0, 40.0),      # hub_length
    (3.0, 3.0),       # num_blades (rounded)

    (0.0, 0.2),       # m
    (0.05, 0.95),     # p
    (0.2, 2.0),       # thickness

    (-30, 30),        # loc2_circ (arc length along hub, not angle)
    (-30, 30),        # loc2_z offset (relative to z0)
    (0.0, 60.0),      # loc2_radial

    (-30, 30),        # loc3_circ
    (-30, 30),        # loc3_z
    (0.0, 60.0),      # loc3_radial

    (-60, 60),        # blade_circ
    (-60, 60),        # blade_z

    (-4*np.pi, 4*np.pi),  # a_AoA
    (-4*np.pi, 4*np.pi),  # b_AoA
    (-4*np.pi, 4*np.pi),  # c_AoA
    (-4*np.pi, 4*np.pi),  # d_AoA
    (-4*np.pi, 4*np.pi),  # e_AoA

    (-5.0, 5.0),      # a_scX
    (-5.0, 5.0),      # b_scX
    (-5.0, 5.0),      # c_scX
    (-5.0, 5.0),      # d_scX
    (0.1, 6.0),       # e_scX

    (-5.0, 5.0),      # a_scY
    (-5.0, 5.0),      # b_scY
    (-5.0, 5.0),      # c_scY
    (-5.0, 5.0),      # d_scY
    (0.1, 6.0),       # e_scY

    (0.0, 0.0),       # apply_thickness_normal
]


# ============================================================
# 6) Cost function (FAST)
# ============================================================

def cost_fn(x):
    params = vector_to_params(x)

    # Quick sanity constraints to avoid nonsense
    # (keeps optimizer from wasting time)
    if not (0.0 <= params["m"] <= 0.25):
        return 1e9
    if not (0.01 <= params["p"] <= 0.99):
        return 1e9
    if params["hub_radius"] <= 0.1 or params["hub_length"] <= 0.1:
        return 1e9

    try:
        cand_pts = candidate_pointcloud(params, n_s=CAND_SAMPLES_S, n_t=CAND_SAMPLES_T)
    except Exception as e:
        # numerical crash -> big penalty
        return 1e9

    # Chamfer distance (main)
    cham = chamfer_distance(cand_pts, target_points)

    # Diameter penalty
    cand_d = bbox_diameter(cand_pts)
    diam_err = abs(cand_d - TARGET_DIAM) / max(TARGET_DIAM, 1e-9)

    cost = W_CHAMFER * cham + W_DIAM * diam_err

    print(f"cost={cost:.6e} cham={cham:.6e} diam_err={diam_err:.4f} global_scale={int(round(params['global_scale']))}")
    return float(cost)


# ============================================================
# 7) Optimize + export final STEP (one time)
# ============================================================

def export_best_as_step(best_params, filename="optimized_propeller.step"):
    """
    Only at the end do we build the full CadQuery solid and export STEP.
    (Still uses your style, but only once.)
    """
    # If you want, you can paste your original CadQuery build here.
    # For now, we’ll generate the final using a simplified approach:
    # - Create your original CadQuery solid once, at high resolution, then scale and export.
    #
    # Easiest: reuse your original script structure but parameterized.
    #
    # To keep this file focused, we’ll just show a minimal export:
    #
    # If you want me to integrate your original CadQuery loft+hub assembly here,
    # tell me and I’ll paste it in fully.

    print("NOTE: export_best_as_step is a stub in this fast script.")
    print("Best params:", best_params)


def run():
    result = differential_evolution(
        cost_fn,
        bounds=bounds,
        maxiter=60,
        popsize=12,
        tol=1e-3,
        seed=42,
        workers=-1,
        updating="deferred",
        polish=True,
    )

    print("\n=== Optimization finished ===")
    print(result)

    best_x = result.x
    best_params = vector_to_params(best_x)

    # Force discrete params nicely
    best_params["num_blades"] = int(round(best_params["num_blades"]))
    best_params["apply_thickness_normal"] = int(round(best_params["apply_thickness_normal"]))

    print("\nBest params:")
    for k, v in best_params.items():
        print(f"  {k}: {v}")

    # OPTIONAL: build/export CAD once at the end
    export_best_as_step(best_params, filename="optimized_propeller.step")


if __name__ == "__main__":
    run()
