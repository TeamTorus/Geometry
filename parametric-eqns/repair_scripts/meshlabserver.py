import pymeshlab

# Initialize a new MeshSet object
ms = pymeshlab.MeshSet()

# Load the STL file into the MeshSet
ms.load_new_mesh("toroidal_propeller.stl")

# Check if the mesh is manifold
print("The mesh is not manifold. Repairing...")

ms.apply_filter("meshing_repair_non_manifold_edges")
ms.apply_filter("meshing_repair_non_manifold_vertices")

ms.apply_filter("meshing_close_holes")
ms.apply_filter("meshing_merge_close_vertices")
ms.apply_filter("meshing_snap_mismatched_borders")
# ms.apply_filter("meshing_re_orient_faces_coherently")

print("Repair complete.")

# Save the repaired mesh
ms.save_current_mesh("repaired_model3.stl")
