/**
 * @file generate_ditorus.cc
 * @brief Ditorus (3-torus) S^1 x S^1 x S^1 as a closed tetrahedral 3-manifold
 *        in R^4.
 *
 * Pipeline:
 *   1.  Base factors: three circles (inner, middle, outer).
 *   2a. inner x middle  -> torus quad surface in R^4 ([0,1]=inner,[2,3]=middle).
 *   2b. Triangulate the quad faces.
 *   3a. (inner x middle) x outer -> triangular prisms in R^6 ([4,5]=outer).
 *   3b. Triangulate the prisms into tetrahedra.
 *   3c. Topology + Betti checks on the coarse mesh (cheap; before subdivision).
 *   4.  Optional subdivision (barycentric 1->24 or uniform 1->8) + "block"
 *       normalisation ({2,2,2}) after each iteration.
 *   5.  Project R^6 -> R^4 via the ditorus parametrisation.
 *   6.  Consistent orientation (BFS) on the R^4 mesh.
 *   7.  Final-mesh checks: closed, non-degenerate, on-surface, volume.
 *   8.  Export mesh-quality CSVs (quality_metrics.hh).
 *   9.  Visualise in R^3 (drop w, show it as colour).
 */

#include <cctype>
#include <iostream>
#include <string>

#include "ovm_nd_types.hh"
#include "base_shape_generator.hh"
#include "cartesian_product.hh"
#include "nd_subdivision.hh"
#include "orientation.hh"
#include "validation.hh"
#include "homology.hh"
#include "projection.hh"
#include "quality_metrics.hh"
#include "visualize_4d.hh"
#include "results_dir.hh"


static int prompt_segments(const char* which) {
    int n = 8;
    std::cout << "  " << which << " circle segments (>= 3): ";
    while (!(std::cin >> n) || n < 3) {
        std::cout << "  Invalid input. Please enter an integer >= 3: ";
        std::cin.clear(); std::cin.ignore(10000, '\n');
    }
    return n;
}


static int prompt_drop_axis() {
    char c = 'w';
    std::cout << "  Which coordinate to drop for the 3D view? (x, y, z, w): ";
    while (std::cin >> c) {
        switch (std::tolower(static_cast<unsigned char>(c))) {
            case 'x': return 0;
            case 'y': return 1;
            case 'z': return 2;
            case 'w': return 3;
        }
        std::cout << "  Invalid input. Enter x, y, z, or w: ";
    }
    return 3;
}


int main(int argc, char** argv) {
    // --debug (or -d) enables the full validation + quality-metrics
    // without it the program only builds the shape and opens the viewer.
    bool debug = false;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        if (a == "--debug" || a == "-d") debug = true;
    }

    std::cout << "=== Ditorus  S^1 x S^1 x S^1  (R^6 -> R^4) ===\n";
    std::cout << (debug
        ? "(debug mode: topology, Betti and quality checks enabled)\n"
        : "(normal mode: build + visualize; pass --debug for full checks)\n");

    const int n_in  = prompt_segments("Inner");
    const int n_mid = prompt_segments("Middle");
    const int n_out = prompt_segments("Outer");


    //   Conditions: R1 > R2 + r   and   R2 > r > 0.
    double R1 = 4.0;
    std::cout << "  Outer radius R1 [> 0]: ";
    while (!(std::cin >> R1) || R1 <= 0.0) {
        std::cout << "  Invalid input. Please enter a number > 0: ";
        std::cin.clear(); std::cin.ignore(10000, '\n');
    }

    double R2 = 2.0;
    std::cout << "  Middle radius R2 (0 < R2 < R1 = " << R1 << "): ";
    while (!(std::cin >> R2) || R2 <= 0.0 || R2 >= R1) {
        std::cout << "  Invalid input. Please enter a number with 0 < R2 < R1: ";
        std::cin.clear(); std::cin.ignore(10000, '\n');
    }

    double r = 0.5;
    std::cout << "  Tube radius r (need 0 < r < R2 and R1 > R2 + r; R2 = "
              << R2 << "): ";
    while (!(std::cin >> r) || r <= 0.0 || r >= R2 || R1 <= R2 + r) {
        std::cout << "  Invalid input. Need 0 < r < R2 and R1 > R2 + r: ";
        std::cin.clear(); std::cin.ignore(10000, '\n');
    }

    int iterations = 0;
    std::cout << "  Subdivision iterations (0 to skip): ";
    while (!(std::cin >> iterations) || iterations < 0) {
        std::cout << "  Invalid input. Please enter an integer >= 0: ";
        std::cin.clear(); std::cin.ignore(10000, '\n');
    }

    int scheme = 0;
    if (iterations > 0) {
        std::cout << "  Subdivision scheme ([0] barycentric 1->24, [1] uniform 1->8): ";
        while (!(std::cin >> scheme) || scheme < 0 || scheme > 1) {
            std::cout << "  Invalid input. Please enter 0 or 1: ";
            std::cin.clear(); std::cin.ignore(10000, '\n');
        }
    }
    std::cout << "\n";


    // 1. Base factors (each circle starts on its unit S^1).
    auto s1_in  = ShapeGen::circle(n_in);   // inner tube  (radius r)
    auto s1_mid = ShapeGen::circle(n_mid);   // middle loop (radius R2)
    auto s1_out = ShapeGen::circle(n_out);   // outer loop  (radius R1)


    // 2a. inner x middle -> torus quad surface in R^4.
    //     [0,1] = inner, [2,3] = middle.
    auto comb_in_mid = [](const Geometry::Vec4d& a, const Geometry::Vec4d& b) {
        return Geometry::Vec4d(a[0], a[1], b[0], b[1]);
    };
    std::cout << "Computing inner x middle ...\n";
    Geometry::Mesh4D torus =
        CartesianOps::cartesian_product<Geometry::Mesh4D>(s1_in, s1_mid, comb_in_mid);
    Validation::print_stats(torus, "inner x middle (quad surface)");


    // 2b. Triangulate the quad faces so the next product generates triangular
    //     prisms instead of hexahedra.
    std::cout << "\n";
    Subdivision::triangulate_quad_faces(torus);
    Validation::print_stats(torus, "after quad triangulation");


    // 3a. (inner x middle) x outer -> triangular prisms in R^6.
    //     [4,5] = outer.

    auto comb_out = [](const Geometry::Vec4d& ab, const Geometry::Vec4d& c) {
        return Geometry::Vec6d(ab[0], ab[1], ab[2], ab[3], c[0], c[1]);
    };
    std::cout << "\nComputing (inner x middle) x outer ...\n";
    Geometry::Mesh6D mesh =
        CartesianOps::cartesian_product<Geometry::Mesh6D>(torus, s1_out, comb_out);
    Validation::print_stats(mesh, "Raw product (triangular prisms)");


    // 3b. Triangulate prisms -> tetrahedra.

    std::cout << "\n";
    Subdivision::triangulate_prisms_to_tets(mesh);
    Validation::print_stats(mesh, "After prism triangulation");


    // 3c. Topological validation on the mesh before subdivision (debug only).
    if (debug) {
        std::cout << "\n--- Topological checks (coarse mesh) ---\n";
        Validation::Topological::verify_topology(mesh);
        std::cout << "\n--- Homology check (Betti numbers, coarse mesh) ---\n";
        // T^3 = S^1 x S^1 x S^1: P(t) = (1 + t)^3 = 1 + 3t + 3t^2 + t^3 -> (1,3,3,1).
        Homology::verify_betti(mesh, {1, 3, 3, 1});
    }


    // 4. Optional barycentric subdivision + per-iteration block normalisation.

    for (int i = 0; i < iterations; ++i) {
        std::cout << "\nSubdivision iteration " << (i + 1) << "/" << iterations << ":\n";
        if (scheme == 1) Subdivision::uniform_subdivide(mesh);
        else             Subdivision::barycentric_subdivide(mesh);
        Subdivision::normalize_blocks(mesh, {2, 2, 2});
        Validation::print_stats(mesh, "  after subdivision");
    }


    // 5. Project R^6 -> R^4 (ditorus embedding).

    std::cout << "\n";
    Geometry::Mesh4D shape = Projection::project_ditorus(mesh, R1, R2, r);
    Validation::print_stats(shape, "Projected ditorus (R^4)");


    // 6. Consistent orientation.
    std::cout << "\n";
    Orientation::orient_mesh(shape);


    // 7. Final-mesh checks + 8. mesh-quality export (debug only).
    if (debug) {
        std::cout << "\n--- Final-mesh checks ---\n";
        Validation::Topological::verify_closed_manifold(shape);
        Validation::Geometric::verify_no_degenerate_cells(shape);
        Validation::Geometric::verify_on_ditorus(shape, R1, R2, r);
        Validation::Geometric::report_volume(
            shape, Validation::Geometric::ditorus_volume(R1, R2, r), "ditorus");

        QualityMetrics::analyze_and_export(
            shape, std::string(RESULTS_DIR) + "/quality/ditorus_"
                       + std::to_string(n_in) + "x" + std::to_string(n_mid) + "x"
                       + std::to_string(n_out) + "_iter" + std::to_string(iterations));
    }

    // 9. Visualise in R^3
    const int drop_axis = prompt_drop_axis();
    Visualize::show_4d(shape, drop_axis, "Ditorus");

    return 0;
}
