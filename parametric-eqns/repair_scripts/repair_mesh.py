import open3d as o3d

# Load the STL file
mesh = o3d.io.read_triangle_mesh("toroidal_propeller.stl")

# Perform basic repairs
mesh.remove_duplicated_vertices()
mesh.remove_duplicated_triangles()
mesh.remove_unreferenced_vertices()

meshee = o3d.t.geometry.TriangleMesh.from_legacy(mesh)

# Optional: Perform hole filling (only for small holes)
meshee = meshee.fill_holes()

mesh = o3d.geometry.TriangleMesh.compute_triangle_normals(meshee.to_legacy())

# Save the repaired STL
o3d.io.write_triangle_mesh("repaired_model.stl", mesh)
