#include <iostream>
#include <vector>
#include <array>
#include <cmath>

#include "ovm_nd_types.hh"
#include "base_shape_generator.hh"
#include "cartesian_product.hh"
#include "nd_subdivision.hh"
#include "orientation.hh"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"


/**
 * @brief Builds the 1-to-4 uniform subdivision of an equilateral triangle as
 *        a 2-disk: 3 corner vertices on the unit circle of radius r, 3 edge
 *        midpoints, and 4 triangular faces (three corner triangles plus the
 *        central inverted triangle).
 *
 * @param r  Circumradius of the equilateral triangle.
 * @return   A Mesh4D containing the subdivided 2-disk (Z = W = 0).
 */
static Geometry::Mesh4D subdivided_triangular_disk(double r) {
    using VH = OpenVolumeMesh::VertexHandle;
    Geometry::Mesh4D m;

    //create triangle on input radius
    std::array<VH, 3> v;
    for (int i = 0; i < 3; ++i) {
        const double theta = 2.0 * M_PI * static_cast<double>(i) / 3.0;
        v[i] = m.add_vertex(Geometry::Vec4d(r * std::cos(theta),
                                            r * std::sin(theta),
                                            0.0, 0.0));
    }

    // mid[i] is the midpoint of the original edge (v[i], v[(i+1)%3]).
    std::array<VH, 3> mid;
    for (int i = 0; i < 3; ++i) {
        const int j = (i + 1) % 3;
        mid[i] = m.add_vertex((m.vertex(v[i]) + m.vertex(v[j])) * 0.5);
    }

    // Corner triangle at v[i] is bounded by the two midpoints on its two
    // adjacent original edges: mid[i] and mid[(i+2)%3].
    for (int i = 0; i < 3; ++i) {
        m.add_face({v[i], mid[i], mid[(i + 2) % 3]});
    }

    // Central inverted triangle joining the three midpoints, wound the same
    // rotational sense as the corner triangles.
    m.add_face({mid[0], mid[1], mid[2]});

    return m;
}


int main() {
    std::cout << "=== Four adjacent triangular prisms (1-to-4 subdivided triangle x I) ===\n\n";


    //   Phase 1: Base factors.
    auto tri  = subdivided_triangular_disk(1.0);
    auto iint = ShapeGen::interval(1.5);

    std::cout << "Subdivided triangle 2-disk: V=" << tri.n_vertices()
              << " E=" << tri.n_edges()
              << " F=" << tri.n_faces()
              << "   (expected 6/9/4)\n";
    std::cout << "Interval:                   V=" << iint.n_vertices()
              << " E=" << iint.n_edges() << "\n\n";


    //   Phase 2: Cartesian product into R^3
    auto combine = [](const Geometry::Vec4d& tri_pt,
                      const Geometry::Vec4d& iv_pt) {
        return Geometry::Vec3d(tri_pt[0], tri_pt[1], iv_pt[0]);
    };

    std::cout << "Computing subdivided-triangle x I ...\n";
    Geometry::Mesh3D prisms =
        CartesianOps::cartesian_product<Geometry::Mesh3D>(tri, iint, combine);

    std::cout << "  V=" << prisms.n_vertices()
              << " E=" << prisms.n_edges()
              << " F=" << prisms.n_faces()
              << " C=" << prisms.n_cells()
              << "   (expected 12/24/17/4)\n\n";


    // Phase 3: Triangulate every prism into 3 tetrahedra
    Subdivision::triangulate_prisms_to_tets(prisms);
    std::cout << "After prism partitioning: V=" << prisms.n_vertices()
              << " E=" << prisms.n_edges()
              << " F=" << prisms.n_faces()
              << " C=" << prisms.n_cells()
              << "   (expected C = 12)\n\n";


    // Phase 4: BFS orientation across all 12 tets
    Orientation::orient_mesh(prisms);
    std::cout << "\n";


    // Phase 5: Visualise as a halfface surface mesh
    polyscope::init();

    std::vector<std::array<double, 3>> verts;
    verts.reserve(prisms.n_vertices());
    for (auto v_it = prisms.vertices_begin(); v_it != prisms.vertices_end(); ++v_it) {
        const auto p = prisms.vertex(*v_it);
        verts.push_back({p[0], p[1], p[2]});
    }

    std::vector<std::vector<size_t>> tris_a;
    std::vector<double> tri_cell_ids;
    tris_a.reserve(static_cast<size_t>(prisms.n_cells()) * 4u);
    tri_cell_ids.reserve(static_cast<size_t>(prisms.n_cells()) * 4u);

    int cid = 0;
    int skipped = 0;
    for (auto c_it = prisms.cells_begin(); c_it != prisms.cells_end(); ++c_it) {
        for (auto chf_it = prisms.chf_iter(*c_it); chf_it.valid(); ++chf_it) {
            std::vector<size_t> tri;
            for (auto hfhe_it = prisms.hfhe_iter(*chf_it); hfhe_it.valid(); ++hfhe_it) {
                tri.push_back(static_cast<size_t>(
                    prisms.from_vertex_handle(*hfhe_it).idx()));
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

    std::vector<std::vector<size_t>> tris_b = tris_a;
    for (auto& t : tris_b) std::swap(t[1], t[2]);

    std::cout << "Rendering " << tris_a.size() << " triangles per orientation ("
              << prisms.n_cells() << " cells x 4 halffaces each).\n";

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

    polyscope::show();
    return 0;
}
