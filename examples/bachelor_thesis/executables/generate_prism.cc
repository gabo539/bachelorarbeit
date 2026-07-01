#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <cmath>

#include "ovm_nd_types.hh"
#include "base_shape_generator.hh"
#include "cartesian_product.hh"
#include "nd_subdivision.hh"
#include "orientation.hh"
#include "validation.hh"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"


/**
 * @brief Builds a filled triangle (a 2-disk, D^2) using 3 vertices and 1 face.
 * Because it includes a filled face, taking the Cartesian product of this triangle with a 1D interval
 * generates a solid 3D prism
 *
 * @param r  Circumradius of the equilateral triangle.
 * @return   A Mesh4D containing the triangular 2-disk (Z = W = 0).
 */
static Geometry::Mesh4D triangular_disk(double r) {
    using VH = OpenVolumeMesh::VertexHandle;
    Geometry::Mesh4D m;
    for (int i = 0; i < 3; ++i) {
        const double theta = 2.0 * M_PI * static_cast<double>(i) / 3.0;
        m.add_vertex(Geometry::Vec4d(r * std::cos(theta),
                                     r * std::sin(theta),
                                     0.0, 0.0));
    }
    m.add_face({VH(0), VH(1), VH(2)});
    return m;
}


int main() {
    std::cout << "=== Closed triangular prism (triangle x I) ===\n";

    // Circumradius of the equilateral triangle base.
    double r = 1.0;
    std::cout << "  Triangle circumradius r [> 0]: ";
    while (!(std::cin >> r) || r <= 0.0) {
        std::cout << "  Invalid input. Please enter a number > 0: ";
        std::cin.clear(); std::cin.ignore(10000, '\n');
    }

    // Length of the interval (height of the prism).
    double length = 1.5;
    std::cout << "  Interval length (prism height) [> 0]: ";
    while (!(std::cin >> length) || length <= 0.0) {
        std::cout << "  Invalid input. Please enter a number > 0: ";
        std::cin.clear(); std::cin.ignore(10000, '\n');
    }

    // Number of subdivision iterations.
    int iterations = 0;
    std::cout << "  Subdivision iterations (0 to skip): ";
    while (!(std::cin >> iterations) || iterations < 0) {
        std::cout << "  Invalid input. Please enter an integer >= 0: ";
        std::cin.clear(); std::cin.ignore(10000, '\n');
    }

    // Subdivision scheme (only asked when iterations > 0).
    int scheme = 0;
    if (iterations > 0) {
        std::cout << "  Subdivision scheme ([0] barycentric 1->24, [1] uniform 1->8): ";
        while (!(std::cin >> scheme) || scheme < 0 || scheme > 1) {
            std::cout << "  Invalid input. Please enter 0 or 1: ";
            std::cin.clear(); std::cin.ignore(10000, '\n');
        }
    }
    std::cout << "\n";


    // Phase 1: Base factors.
    auto tri  = triangular_disk(r);
    auto iint = ShapeGen::interval(length);

    std::cout << "Triangle 2-disk: V=" << tri.n_vertices()
              << " E=" << tri.n_edges()
              << " F=" << tri.n_faces() << "\n";
    std::cout << "Interval:        V=" << iint.n_vertices()
              << " E=" << iint.n_edges() << "\n\n";


    // Phase 2: Cartesian product into R^3.
    //   Triangle gives (x, y); interval gives z.
    auto combine = [](const Geometry::Vec4d& tri_pt,
                      const Geometry::Vec4d& iv_pt) {
        return Geometry::Vec3d(tri_pt[0], tri_pt[1], iv_pt[0]);
    };

    std::cout << "Computing triangle x I ...\n";
    Geometry::Mesh3D prism =
        CartesianOps::cartesian_product<Geometry::Mesh3D>(tri, iint, combine);

    std::cout << "  V=" << prism.n_vertices()
              << " E=" << prism.n_edges()
              << " F=" << prism.n_faces()
              << " C=" << prism.n_cells()
              << "   (expected 6/9/5/1)\n\n";


    // Phase 3a: Triangulate the prism into 3 tetrahedra.
    Subdivision::triangulate_prisms_to_tets(prism);
    std::cout << "After Stage A: V=" << prism.n_vertices()
              << " E=" << prism.n_edges()
              << " F=" << prism.n_faces()
              << " C=" << prism.n_cells()
              << "   (expected C = 3)\n\n";


    // Phase 3b: Optional subdivision iterations.
    // Barycentric multiplies cell count by 24; uniform by 8.
    long expected_C = prism.n_cells();
    for (int i = 0; i < iterations; ++i) {
        const long factor = (scheme == 1) ? 8 : 24;
        expected_C *= factor;
        std::cout << "Subdivision iteration " << (i + 1) << "/" << iterations << ":\n";
        if (scheme == 1) Subdivision::uniform_subdivide(prism);
        else             Subdivision::barycentric_subdivide(prism);
        std::cout << "  V=" << prism.n_vertices()
                  << " E=" << prism.n_edges()
                  << " F=" << prism.n_faces()
                  << " C=" << prism.n_cells()
                  << "   (expected C = " << expected_C << ")\n";
    }
    std::cout << "\n";


    // Phase 3c: Consistent orientation (BFS + rebuild)
    Orientation::orient_mesh(prism);
    std::cout << "\n";


    // Phase 4: Visualise as a halfface surface mesh.
    polyscope::init();

    std::vector<std::array<double, 3>> verts;
    verts.reserve(prism.n_vertices());
    for (auto v_it = prism.vertices_begin(); v_it != prism.vertices_end(); ++v_it) {
        const auto p = prism.vertex(*v_it);
        verts.push_back({p[0], p[1], p[2]});
    }

    std::vector<std::vector<size_t>> tris_a;
    std::vector<double> tri_cell_ids;
    tris_a.reserve(static_cast<size_t>(prism.n_cells()) * 4u);
    tri_cell_ids.reserve(static_cast<size_t>(prism.n_cells()) * 4u);

    int cid = 0;
    int skipped = 0;
    for (auto c_it = prism.cells_begin(); c_it != prism.cells_end(); ++c_it) {
        // Walk the halffaces claimed by this cell. The cell-halfface
        // iterator returns them in the order they were registered with
        // add_cell during the orientation rebuild, with their hf0/hf1
        // identity reflecting the boundary-operator + bit choice.
        for (auto chf_it = prism.chf_iter(*c_it); chf_it.valid(); ++chf_it) {
            // Read vertices in the halfface's winding order by walking
            // its halfedges and taking each one's "from" vertex.
            std::vector<size_t> tri;
            for (auto hfhe_it = prism.hfhe_iter(*chf_it); hfhe_it.valid(); ++hfhe_it) {
                tri.push_back(static_cast<size_t>(
                    prism.from_vertex_handle(*hfhe_it).idx()));
            }
            if (tri.size() != 3) { ++skipped; continue; }

            tris_a.push_back(std::move(tri));
            tri_cell_ids.push_back(static_cast<double>(cid));
        }
        ++cid;
    }
    if (skipped > 0) {
        std::cout << "  [warn] Skipped " << skipped
                  << " non-triangular halfface(s) during render.\n";
    }

    // Build orientation B by reversing every triangle's winding.
    std::vector<std::vector<size_t>> tris_b = tris_a;
    for (auto& t : tris_b) std::swap(t[1], t[2]);

    std::cout << "Rendering " << tris_a.size() << " triangles per orientation ("
              << prism.n_cells() << " cells x 4 halffaces each, minus "
              << "any non-triangular faces).\n";

    auto* sm_a = polyscope::registerSurfaceMesh(
                     "Orientation A (BFS-derived)", verts, tris_a);
    sm_a->addFaceScalarQuantity("Source cell", tri_cell_ids);
    sm_a->setBackFacePolicy(polyscope::BackFacePolicy::Cull);
    sm_a->setEnabled(true);

    auto* sm_b = polyscope::registerSurfaceMesh(
                     "Orientation B (flipped)", verts, tris_b);
    sm_b->addFaceScalarQuantity("Source cell", tri_cell_ids);
    sm_b->setBackFacePolicy(polyscope::BackFacePolicy::Cull);
    sm_b->setEnabled(false);

    polyscope::show();   // Blocks until the user closes the window.
    return 0;
}
