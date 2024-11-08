#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/repair.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

int main() {
    // Load the STL mesh
    std::ifstream input("toroidal_propeller.stl");
    Polyhedron mesh;
    input >> mesh;

    // Remove non-manifold edges
    CGAL::Polygon_mesh_processing::remove_degenerate_faces(mesh);
    CGAL::Polygon_mesh_processing::stitch_borders(mesh);
    CGAL::Polygon_mesh_processing::remove_isolated_vertices(mesh);

    // Save the repaired mesh
    std::ofstream output("repaired_model.stl");
    output << mesh;

    return 0;
}
