# python .\repair_mesh_bpy.py toroidal_propeller.stl repaired_model2.stl

import bpy
import sys

def delete_all_objects():
    bpy.ops.object.select_all(action='SELECT')
    bpy.ops.object.delete(use_global=False)

def import_stl(filepath):
    bpy.ops.import_mesh.stl(filepath=filepath)

def export_stl(filepath):
    bpy.ops.export_mesh.stl(filepath=filepath)

def repair_mesh():
    obj = bpy.context.active_object
    bpy.ops.object.mode_set(mode='EDIT')
    
    # Select all geometry
    bpy.ops.mesh.select_all(action='SELECT')
    
    # Remove doubles (merge by distance)
    bpy.ops.mesh.remove_doubles(threshold=0.00001)
    
    # Recalculate normals
    bpy.ops.mesh.normals_make_consistent(inside=False)
    
    
    # Remove non-manifold edges
    bpy.ops.mesh.select_non_manifold()
    bpy.ops.mesh.delete(type='EDGE')
    

def main():
    # Parse command line arguments
    argv = sys.argv
    argv = argv[1:]
    
    input_file = argv[0]
    output_file = argv[1]

    print("Input file: ", input_file)
    print("Output file: ", output_file)
    
    # Clear existing objects
    delete_all_objects()
    
    # Import STL
    import_stl(input_file)
    
    # Select the imported object
    obj = bpy.context.selected_objects[0]
    bpy.context.view_layer.objects.active = obj
    
    # Repair the mesh
    repair_mesh()
    
    # Export the repaired STL
    export_stl(output_file)

if __name__ == "__main__":
    main()
