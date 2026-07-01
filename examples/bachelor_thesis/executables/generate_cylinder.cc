/**
 * @file generate_cylinder.cc
 * @brief Low-dimensional validation of the Cartesian-product algorithm:
 *        the open cylinder S^1 x I as a 2-complex in R^3.
 *
 *
 * Parameters are prompted for at the console; pressing Enter accepts the
 * defaults shown in square brackets.
 */

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <array>

#include "ovm_nd_types.hh"
#include "base_shape_generator.hh"
#include "cartesian_product.hh"

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"


int main() {

    // Parameter prompts.
    std::cout << "=== Open cylinder S^1 x I ===\n";

    // Number of circle segments.
    int n = 8;
    std::cout << "  Number of circle segments (n >= 3) [Type number, press Enter]: ";
    while (!(std::cin >> n) || n < 3) {
        std::cout << "  Invalid input. Please enter an integer >= 3: ";
        std::cin.clear();
        std::cin.ignore(10000, '\n');
    }

    // Interval length.
    double length = 2.0;
    std::cout << "  Interval length [Type number, press Enter]: ";
    while (!(std::cin >> length) || length <= 0.0) {
        std::cout << "  Invalid input. Please enter a number > 0: ";
        std::cin.clear();
        std::cin.ignore(10000, '\n');
    }

    std::cout << "\n";


    // Phase 1: Base factors S^1 and I.
    auto s1   = ShapeGen::circle(n);
    auto iint = ShapeGen::interval(length);

    std::cout << "S^1 (n-gon): V=" << s1.n_vertices()
              << "  E=" << s1.n_edges() << "\n";
    std::cout << "I (interval): V=" << iint.n_vertices()
              << "  E=" << iint.n_edges() << "\n\n";


    // Phase 2: Combiner.
    auto combine = [](const Geometry::Vec4d& circle_pt,
                      const Geometry::Vec4d& interval_pt) {
        return Geometry::Vec3d(circle_pt[0], circle_pt[1], interval_pt[0]);
    };


    // Phase 3: Cartesian product.

    std::cout << "Computing S^1 x I ...\n";
    Geometry::Mesh3D cylinder =
        CartesianOps::cartesian_product<Geometry::Mesh3D>(s1, iint, combine);

    const long V = cylinder.n_vertices();
    const long E = cylinder.n_edges();
    const long F = cylinder.n_faces();
    const long C = cylinder.n_cells();
    const long chi = V - E + F - C;


    std::cout << "  V = " << V << "   (expected " << 2 * n << ")\n";
    std::cout << "  E = " << E << "   (expected " << 3 * n << ")\n";
    std::cout << "  F = " << F << "   (expected " << n << ")\n";
    std::cout << "  C = " << C << "   (expected 0)\n";
    std::cout << "  chi = " << chi << "   (expected 0)\n";

    const bool ok = (V == 2 * n) && (E == 3 * n) && (F == n)
                 && (C == 0)     && (chi == 0);
    std::cout << (ok ? "  [PASS]\n" : "  [FAIL]\n");


    // Phase 4: Visualise in Polyscope.
    polyscope::init();

    // Convert the OVM mesh into Polyscope's (vertices, faces) format.
    std::vector<std::array<double, 3>> vertices;
    vertices.reserve(cylinder.n_vertices());
    for (auto v_it = cylinder.vertices_begin(); v_it != cylinder.vertices_end(); ++v_it) {
        const auto p = cylinder.vertex(*v_it);
        vertices.push_back({p[0], p[1], p[2]});
    }

    std::vector<std::vector<size_t>> faces;
    faces.reserve(cylinder.n_faces());
    for (auto fa_it = cylinder.faces_begin(); fa_it != cylinder.faces_end(); ++fa_it) {
        std::vector<size_t> face_indices;
        for (auto fv_it = cylinder.fv_iter(*fa_it); fv_it.valid(); ++fv_it) {
            face_indices.push_back(static_cast<size_t>(fv_it->idx()));
        }
        faces.push_back(std::move(face_indices));
    }

    polyscope::registerSurfaceMesh("Cylinder S^1 x I", vertices, faces);
    polyscope::show();   // Blocks until the user closes the window.

    return 0;
}
