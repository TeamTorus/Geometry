import trimesh
import networkx as nx

def load_mesh(filename):
    mesh = trimesh.load(filename, force='mesh')
    return mesh

def identify_non_manifold(mesh):
    # Find edges with more than two adjacent faces
    edges = mesh.edges_unique
    edge_faces = mesh.edges_unique_faces
    non_manifold_edges = edges[mesh.edges_unique_faces.shape[1] > 2]
    return non_manifold_edges

def remove_non_manifold_faces(mesh):
    # Identify and remove faces that are part of non-manifold edges
    non_manifold_edges = identify_non_manifold(mesh)
    if len(non_manifold_edges) == 0:
        print("No non-manifold edges found.")
        return mesh
    
    # Find all faces adjacent to non-manifold edges
    faces_to_remove = set()
    for edge in non_manifold_edges:
        adjacent_faces = mesh.face_adjacency_edges == edge
        faces_to_remove.update(mesh.face_adjacency_edges[adjacent_faces].flatten())
    
    # Remove the identified faces
    mesh.update_faces(~mesh.faces_unique[faces_to_remove])
    mesh.remove_unreferenced_vertices()
    return mesh

def main(input_file, output_file):
    mesh = load_mesh(input_file)
    print("Original mesh is watertight:", mesh.is_watertight)
    
    mesh = remove_non_manifold_faces(mesh)
    print("After removing non-manifold faces, mesh is watertight:", mesh.is_watertight)
    
    # Attempt to fill any remaining holes
    mesh.fill_holes()
    
    # Final check
    print("Final mesh is watertight:", mesh.is_watertight)
    
    # Export the repaired mesh
    mesh.export(output_file)
    print(f"Repaired STL saved to {output_file}")

if __name__ == "__main__":
    print("Starting repair process...")
    main("toroidal_propeller.stl", "repaired_model.stl")
