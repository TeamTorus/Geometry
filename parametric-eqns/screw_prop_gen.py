# screw_propeller_gen.py
# -------------------------------------------------------------
import numpy as np
import sympy as sp
import pyvista as pv
from scipy.spatial import ConvexHull

# ---------- tiny helper: area of projected toroidal blade ----------
def _projected_area(x, y):
    """Return area of convex hull of (x,y) points."""
    pts = np.column_stack((x.ravel(), y.ravel()))
    hull = ConvexHull(pts)
    verts = pts[hull.vertices]
    # polygon area by shoelace
    xh, yh = verts[:, 0], verts[:, 1]
    return 0.5*np.abs(np.dot(xh, np.roll(yh, -1)) -
                      np.dot(yh, np.roll(xh, -1)))

# ---------- NACA 4-digit section (same as your toroidal code) ----------
def _naca_xy(t_sym, m, p, tau):
    t = sp.symbols("t", real=True)
    yt = 5*tau*(0.2969*sp.sqrt(t) - 0.1260*t - 0.3516*t**2
                + 0.2843*t**3 - 0.1036*t**4)
    yc = sp.Piecewise(((m/p**2)*(2*p*t - t**2), t <= p),
                      ((m/(1-p)**2)*((1-2*p)+2*p*t - t**2), t > p))
    xu = t; xl = 1 - t          # unit chord, 0→1 upper, 1→0 lower
    yu = yc + yt
    yl = yc - yt
    x_final = sp.Piecewise((xu, t <= 1), (xl, t > 1))
    y_final = sp.Piecewise((yu, t <= 1), (yl, t > 1))
    return sp.lambdify(t, x_final.subs(t, t_sym), "numpy"), \
           sp.lambdify(t, y_final.subs(t, t_sym), "numpy")

# ---------- screw-propeller generator ----------
def generate_equivalent_screw_propeller(
        toroidal_builder,
        toroidal_kwargs,
        naca_m=0.02, naca_p=0.4, naca_tau=0.12,
        screw_pitch_deg=25,
        s_samples=120, t_samples=120,
        constant_chord=True):
    """
    Parameters
    ----------
    toroidal_builder : callable
        Your function that returns (X_prop, Y_prop, Z_prop) for ONE blade.
        Must accept **toroidal_kwargs.
    toroidal_kwargs : dict
        Exact parameters you already pass to build the toroidal blade.
    Returns
    -------
    screw_mesh : pyvista.PolyData
        A ready-to-use mesh of the screw propeller (all blades + hub).
    info : dict
        Useful numbers: BAR_target, chord_len, R_root, R_tip …
    """
    # 1. Build one toroidal blade and work out its BAR -------------------
    X, Y, Z, hub_R, hub_L, n_blades = toroidal_builder(**toroidal_kwargs,
                                                       return_meta=True)
    # Use the first blade (columns 0-s_resolution)
    s_res = X.shape[1] // n_blades
    Xb, Yb = X[:, :s_res], Y[:, :s_res]
    A_proj = _projected_area(Xb, Yb)          # projected blade area
    R_tip = np.max(np.sqrt(X**2 + Y**2))      # outer radius
    BAR_target = A_proj / (np.pi*R_tip**2)    # per-blade BAR

    # 2. Decide screw-blade chord distribution --------------------------
    R_root = hub_R
    span = R_tip - R_root
    if constant_chord:
        chord_len = BAR_target*np.pi*R_tip**2 / span
        chord = lambda r: chord_len*np.ones_like(r)
    else:
        # elliptical planform (optional)
        C0 = BAR_target*np.pi*R_tip**2 / span
        chord = lambda r: C0*np.sqrt(1 - ((r-R_root)/span)**2)

    # 3. Build screw-blade mesh -----------------------------------------
    t_arr = np.linspace(0, 2, t_samples)         # t for NACA
    s_arr = np.linspace(0, 1, s_samples)         # span param

    # pre-compute NACA cross-section
    naca_x, naca_y = _naca_xy(t_arr, naca_m, naca_p, naca_tau)
    x2d = naca_x(t_arr)
    y2d = naca_y(t_arr)

    # local pitch / twist (constant here)
    pitch = np.deg2rad(screw_pitch_deg)

    # arrays to accumulate blade
    X_bld = np.zeros((t_samples, s_samples))
    Y_bld = np.zeros_like(X_bld)
    Z_bld = np.zeros_like(X_bld)

    for j, s in enumerate(s_arr):
        r = R_root + s*span
        c = chord(r)
        # scale 2D section to local chord
        xs = (x2d - 0.25)*c        # place 25 %-chord at r
        ys = y2d*c
        # rotate by pitch (about +z)
        xp =  xs*np.cos(pitch) - ys*np.sin(pitch)
        zp =  xs*np.sin(pitch) + ys*np.cos(pitch)
        # position in cylindrical coords then to Cartesian
        X_bld[:, j] = (r + xp)*1.0
        Y_bld[:, j] = 0.0          # first blade sits in +x half-plane
        Z_bld[:, j] = zp

    # 4. Replicate around hub -------------------------------------------
    meshes = []
    for k in range(n_blades):
        theta = 2*np.pi*k/n_blades
        Rmat = np.array([[np.cos(theta), -np.sin(theta), 0],
                         [np.sin(theta),  np.cos(theta), 0],
                         [0,              0,             1]])
        pts = np.stack((X_bld, Y_bld, Z_bld), axis=-1).reshape(-1, 3)
        pts_rot = (Rmat @ pts.T).T
        mesh = pv.PolyData(pts_rot).delaunay_2d()
        meshes.append(mesh)

    screw_blades = pv.merge(meshes)

    # 5. Quick hub (open cylinder) --------------------------------------
    theta = np.linspace(0, 2*np.pi, 80)
    z = np.linspace(-hub_L/2, hub_L/2, 2)
    th, zz = np.meshgrid(theta, z)
    Xh = hub_R*np.cos(th); Yh = hub_R*np.sin(th); Zh = zz
    hub = pv.StructuredGrid(Xh, Yh, Zh).triangulate()

    screw_mesh = pv.merge([screw_blades, hub]).clean()

    info = dict(BAR_target=BAR_target, chord=np.mean(chord(np.array([0.5]))),
                R_root=R_root, R_tip=R_tip)

    return screw_mesh, info

# -------------------------------------------------------------
if __name__ == "__main__":
    from toroidal_propeller_gen import build_one_blade  # your existing builder

    params = dict(hub_radius=5, hub_length=20, num_blades=3,
                  m=0.02, p=0.4, thickness=0.12,
                  # … the rest of your toroidal kwargs …
                  )

    screw, meta = generate_equivalent_screw_propeller(build_one_blade,
                                                      params)
    print("Target BAR =", meta["BAR_target"])
    screw.plot(show_edges=True)
    screw.save("equiv_screw_propeller.stl")
