/**
 * @file generate_spheritorus.cc
 * @brief Spheritorus S^2 x S^1 as a closed tetrahedral 3-manifold in R^4.
 *
 * Pipeline (stages live in nd_subdivision.hh / orientation.hh / projection.hh):
 *   1.  Base factors:  S^2 = octahedron, S^1 = circle(m).
 *   2.  Cartesian product into R^5 (pure coordinate stacking) -> triangular
 *       prisms.
 *   3.  Triangulate the prisms into tetrahedra.
 *   3b. Topology + Betti checks on the coarse mesh (cheap; before subdivision).
 *   4.  Optional subdivision (barycentric 1->24 or uniform 1->8); after each
 *       each iteration, normalise the factor blocks {3,2} onto the round S^2 / S^1.
 *   5.  Project R^5 -> R^4 via the spheritorus parametrisation.
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

    std::cout << "=== Spheritorus  S^2 x S^1  (R^4) ===\n";
    std::cout << (debug
        ? "(debug mode: topology, Betti and quality checks enabled)\n"
        : "(normal mode: build + visualize; pass --debug for full checks)\n");

    // S^1 resolution. The S^2 (octahedron) is refined only if subdivision is applied.
    int m = 16;
    std::cout << "  Circle S^1 segments (m >= 3) [type a number, Enter]: ";
    while (!(std::cin >> m) || m < 3) {
        std::cout << "  Invalid input. Please enter an integer >= 3: ";
        std::cin.clear(); std::cin.ignore(10000, '\n');
    }

    // Major radius R (offset of the tube centre from the z=w=0 axis).
    double R = 2.0;
    std::cout << "  Major radius R (tube-centre offset) [> 0]: ";
    while (!(std::cin >> R) || R <= 0.0) {
        std::cout << "  Invalid input. Please enter a number > 0: ";
        std::cin.clear(); std::cin.ignore(10000, '\n');
    }

    // Minor radius r. The embedding requires 0 < r < R;
    double r = 0.7;
    std::cout << "  Minor radius r ( 0 < r < R = " << R << "): ";
    while (!(std::cin >> r) || r <= 0.0 || r >= R) {
        std::cout << "  Invalid input. Please enter a number with 0 < r < R = "
                  << R << ": ";
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


    // 1. Base factors.
    auto s2 = ShapeGen::octahedron(1.0);    // chi(S^2) = 2
    auto s1 = ShapeGen::circle(m);     // chi(S^1) = 0
    Validation::print_stats(s2, "Base S^2 (octahedron)");
    Validation::print_stats(s1, "Base S^1 (circle)");


    // 2. Cartesian product into R^5: S^2 -> coords [0,1,2], S^1 -> [3,4].
    auto combine = [](const Geometry::Vec4d& a, const Geometry::Vec4d& b) {
        return Geometry::Vec5d(a[0], a[1], a[2], b[0], b[1]);
    };
    std::cout << "\nComputing S^2 x S^1 ...\n";
    Geometry::Mesh5D mesh =
        CartesianOps::cartesian_product<Geometry::Mesh5D>(s2, s1, combine);
    Validation::print_stats(mesh, "Raw product (triangular prisms)");


    // 3. Triangulate prisms -> tetrahedra.
    std::cout << "\n";
    Subdivision::triangulate_prisms_to_tets(mesh);
    Validation::print_stats(mesh, "After prism triangulation");


    // 3b. Topological validation on the mesh before subdivision (debug only).
    if (debug) {
        std::cout << "\n--- Topological checks (coarse mesh) ---\n";
        Validation::Topological::verify_topology(mesh);
        std::cout << "\n--- Homology check (Betti numbers, coarse mesh) ---\n";
        // S^2 x S^1: P(t) = (1 + t^2)(1 + t) = 1 + t + t^2 + t^3 -> (1,1,1,1).
        Homology::verify_betti(mesh, {1, 1, 1, 1});
    }


    // 4. Optional barycentric subdivision, normalising the factor blocks
    //    after each iteration.
    for (int i = 0; i < iterations; ++i) {
        std::cout << "\nSubdivision iteration " << (i + 1) << "/" << iterations << ":\n";
        if (scheme == 1) Subdivision::uniform_subdivide(mesh);
        else             Subdivision::barycentric_subdivide(mesh);
        Subdivision::normalize_blocks(mesh, {3, 2});
        Validation::print_stats(mesh, "  after iteration");
    }


    // 5. Project R^5 -> R^4
    std::cout << "\n";
    Geometry::Mesh4D shape = Projection::project_spheritorus(mesh, R, r);
    Validation::print_stats(shape, "Projected spheritorus (R^4)");


    // 6. Consistent orientation.
    std::cout << "\n";
    Orientation::orient_mesh(shape);


    // 7. Final-mesh checks + 8. mesh-quality export (debug only).
    if (debug) {
        std::cout << "\n--- Final-mesh checks ---\n";
        Validation::Topological::verify_closed_manifold(shape);
        Validation::Geometric::verify_no_degenerate_cells(shape);
        Validation::Geometric::verify_on_spheritorus(shape, R, r);
        Validation::Geometric::report_volume(
            shape, Validation::Geometric::spheritorus_volume(R, r), "spheritorus");

        QualityMetrics::analyze_and_export(
            shape, std::string(RESULTS_DIR) + "/quality/spheritorus_m"
                       + std::to_string(m) + "_iter" + std::to_string(iterations));
    }

    // 9. Visualise in R^3 .
    const int drop_axis = prompt_drop_axis();
    Visualize::show_4d(shape, drop_axis, "Spheritorus");

    return 0;
}
