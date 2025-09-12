# union_with_bpy.py
import sys, math
import bpy, bmesh, addon_utils
from mathutils import Vector

def union_stl_files(stl_path_a, stl_path_b, output_path):
    """
    Perform boolean union on two STL files using Blender's EXACT solver.
    
    Args:
        stl_path_a (str): Path to first STL file
        stl_path_b (str): Path to second STL file  
        output_path (str): Path for output STL file
    
    Returns:
        bool: True if successful, False otherwise
    """
    
    # --- clean scene (fresh) ---
    bpy.ops.wm.read_factory_settings(use_empty=True)

    # --- enable STL add-on (needed for import/export) ---
    addon_utils.enable("io_mesh_stl", default_set=True, persistent=True)

    def set_active(obj):
        bpy.context.view_layer.objects.active = obj
        for o in bpy.context.selected_objects:
            o.select_set(False)
        obj.select_set(True)

    def ensure_object_mode():
        if bpy.context.object and bpy.context.object.mode != 'OBJECT':
            bpy.ops.object.mode_set(mode='OBJECT')

    def mesh_cleanup(obj, merge=1e-6, triangulate=True):
        """Triangulate, weld near-duplicates, and fix normals."""
        set_active(obj); ensure_object_mode()
        if obj.type != 'MESH':
            bpy.ops.object.convert(target='MESH')
        me = obj.data
        # IMPORTANT: verbose must be keyword
        me.validate(verbose=False)

        bm = bmesh.new()
        bm.from_mesh(me)

        if triangulate:
            bmesh.ops.triangulate(bm, faces=bm.faces[:], quad_method='BEAUTY', ngon_method='BEAUTY')

        # Merge by distance (aka remove doubles)
        bmesh.ops.remove_doubles(bm, verts=bm.verts[:], dist=merge)

        # Make normals consistent
        bmesh.ops.recalc_face_normals(bm, faces=bm.faces[:])

        bm.to_mesh(me)
        bm.free()
        me.update()
        bpy.context.view_layer.update()

    def count_defects(obj):
        bm = bmesh.new()
        bm.from_mesh(obj.data)
        bm.normal_update()
        boundary = sum(1 for e in bm.edges if e.is_boundary)      # edges used by 1 face (holes)
        nonman   = sum(1 for e in bm.edges if not e.is_manifold)  # edges used by 0 or >2 faces
        bm.free()
        return boundary, nonman

    def bbox_diag(obj):
        deps = bpy.context.evaluated_depsgraph_get()
        eo = obj.evaluated_get(deps)
        pts = [eo.matrix_world @ Vector(c) for c in eo.bound_box]
        xs, ys, zs = zip(*[(p.x, p.y, p.z) for p in pts])
        dx, dy, dz = (max(xs)-min(xs)), (max(ys)-min(ys)), (max(zs)-min(zs))
        return math.sqrt(dx*dx + dy*dy + dz*dz)

    try:
        # --- import STLs ---
        bpy.ops.import_mesh.stl(filepath=stl_path_a)
        objA = bpy.context.selected_objects[0]; objA.name = "A"

        bpy.ops.import_mesh.stl(filepath=stl_path_b)
        objB = bpy.context.selected_objects[0]; objB.name = "B"

        # Apply transforms so geometry is 'real'
        for obj in (objA, objB):
            set_active(obj); ensure_object_mode()
            bpy.ops.object.transform_apply(location=False, rotation=True, scale=True)
            mesh_cleanup(obj, merge=1e-6, triangulate=True)

        print("[Pre] A defects:", *count_defects(objA))
        print("[Pre] B defects:", *count_defects(objB))

        # --- boolean union (Exact solver) applied onto A ---
        set_active(objA); ensure_object_mode()
        mod = objA.modifiers.new(name="UNION_EXACT", type='BOOLEAN')
        mod.operation = 'UNION'
        mod.solver    = 'EXACT'
        mod.object    = objB
        bpy.ops.object.modifier_apply(modifier=mod.name)

        # Remove tool object
        bpy.data.objects.remove(objB, do_unlink=True)

        # Polish & check manifoldness
        mesh_cleanup(objA, merge=1e-6, triangulate=False)
        b, n = count_defects(objA)
        print(f"[Post-boolean] boundary={b} nonmanifold={n}")

        # Optional watertight fallback: voxel remesh if not manifold
        if b > 0 or n > 0:
            set_active(objA); ensure_object_mode()
            vox = max(bbox_diag(objA)/400.0, 1e-6)  # ~400 voxels across bbox -> increase for sharper
            bpy.ops.object.voxel_remesh(voxel_size=vox, adaptivity=0.0)
            mesh_cleanup(objA, merge=1e-6, triangulate=False)
            b2, n2 = count_defects(objA)
            print(f"[After remesh] boundary={b2} nonmanifold={n2}")

        # Export STL
        set_active(objA); ensure_object_mode()
        bpy.ops.export_mesh.stl(filepath=output_path, use_selection=True, ascii=False)
        print("[Done] wrote", output_path)
        
        return True
        
    except Exception as e:
        print(f"[Error] Boolean union failed: {e}")
        return False

# Command line usage (backwards compatible)
if __name__ == "__main__" and len(sys.argv) >= 4:
    A, B, OUT = sys.argv[-3], sys.argv[-2], sys.argv[-1]
    union_stl_files(A, B, OUT)

# Example usage as function:
# union_stl_files("hub.stl", "blades.stl", "propeller.stl")
