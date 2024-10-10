import pyvista as pv
sphere_with_hole = pv.Sphere(end_theta=330)

print(sphere_with_hole.is_manifold)  # False

sphere = sphere_with_hole.fill_holes(1000)  
edges = sphere.extract_feature_edges(
    feature_edges=False, manifold_edges=False
)  

print(sphere.is_manifold)  # True

assert edges.n_cells == 0  