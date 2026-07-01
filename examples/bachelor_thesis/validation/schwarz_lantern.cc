/**
 * @file schwarz_lantern.cc
 * @brief Schwarz-lantern volume convergence study for the spheritorus and the
 *        ditorus: how the measured discrete 3-volume approaches the ideal value
 *        under two ways of refining the mesh.
 *
 *   - segments     : more circle segments, no subdivision -> clean inscription,
 *                    converges to 100% of the ideal value from below.
 *   - subdivision  : a coarse base, refined by subdivision -> can overshoot
 *                    100% (the Schwarz-lantern tilt effect).
 *
 *   Test sets
 *   spheritorus : (R=5, r=1) and (R=5, r=3)
 *   ditorus     : (R1=10, R2=5, r=1) and (R1=10, R2=5, r=4)
 *
 * Writes results/schwarz_lantern/schwarz_lantern.csv.
 */

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "ovm_nd_types.hh"
#include "base_shape_generator.hh"
#include "cartesian_product.hh"
#include "nd_subdivision.hh"
#include "projection.hh"
#include "validation.hh"
#include "results_dir.hh"

namespace {

// Silences the pipeline's progress prints while one mesh is built.
struct CoutMuter {
    std::streambuf* old_buf;
    std::ostringstream sink;
    CoutMuter() : old_buf(std::cout.rdbuf(sink.rdbuf())) {}
    ~CoutMuter() { std::cout.rdbuf(old_buf); }
};

// Builds a projected spheritorus (R^4) for the given parameters.
Geometry::Mesh4D build_spheritorus(int m, double R, double r, int iterations, int scheme) {
    CoutMuter mute;
    auto s2 = ShapeGen::octahedron(1.0);
    auto s1 = ShapeGen::circle(m);
    Geometry::Mesh5D mesh = CartesianOps::cartesian_product<Geometry::Mesh5D>(
        s2, s1, [](const Geometry::Vec4d& a, const Geometry::Vec4d& b) {
            return Geometry::Vec5d(a[0], a[1], a[2], b[0], b[1]);
        });
    Subdivision::triangulate_prisms_to_tets(mesh);
    for (int i = 0; i < iterations; ++i) {
        if (scheme == 1) Subdivision::uniform_subdivide(mesh);
        else             Subdivision::barycentric_subdivide(mesh);
        Subdivision::normalize_blocks(mesh, {3, 2});
    }
    return Projection::project_spheritorus(mesh, R, r);
}

// Builds a projected ditorus (R^4); all three circles share the segment count n.
Geometry::Mesh4D build_ditorus(int n, double R1, double R2, double r, int iterations, int scheme) {
    CoutMuter mute;
    auto s1_in  = ShapeGen::circle(n);
    auto s1_mid = ShapeGen::circle(n);
    auto s1_out = ShapeGen::circle(n);
    Geometry::Mesh4D torus = CartesianOps::cartesian_product<Geometry::Mesh4D>(
        s1_in, s1_mid, [](const Geometry::Vec4d& a, const Geometry::Vec4d& b) {
            return Geometry::Vec4d(a[0], a[1], b[0], b[1]);
        });
    Subdivision::triangulate_quad_faces(torus);
    Geometry::Mesh6D mesh = CartesianOps::cartesian_product<Geometry::Mesh6D>(
        torus, s1_out, [](const Geometry::Vec4d& ab, const Geometry::Vec4d& c) {
            return Geometry::Vec6d(ab[0], ab[1], ab[2], ab[3], c[0], c[1]);
        });
    Subdivision::triangulate_prisms_to_tets(mesh);
    for (int i = 0; i < iterations; ++i) {
        if (scheme == 1) Subdivision::uniform_subdivide(mesh);
        else             Subdivision::barycentric_subdivide(mesh);
        Subdivision::normalize_blocks(mesh, {2, 2, 2});
    }
    return Projection::project_ditorus(mesh, R1, R2, r);
}

struct Row {
    std::string shape, axis, scheme;
    int seg, iterations;
    double R1, R2, r;
    long n_cells;
    double measured, ideal;
};

const char* scheme_name(int iterations, int scheme) {
    if (iterations == 0) return "none";
    return (scheme == 1) ? "uniform" : "barycentric";
}

struct SphSet { double R, r; };
struct DitSet { double R1, R2, r; };
const SphSet kSphSets[2] = {{5.0, 1.0}, {5.0, 3.0}};
const DitSet kDitSets[2] = {{10.0, 5.0, 1.0}, {10.0, 5.0, 4.0}};

} // namespace

int main() {
    std::vector<Row> rows;

    auto add_sph = [&](const std::string& axis, int m, double R, double r,
                       int iterations, int scheme) {
        auto mesh = build_spheritorus(m, R, r, iterations, scheme);
        const double measured = Validation::Geometric::measured_volume(mesh);
        const double ideal = Validation::Geometric::spheritorus_volume(R, r);
        rows.push_back({"spheritorus", axis, scheme_name(iterations, scheme), m, iterations,
                        R, 0.0, r, static_cast<long>(mesh.n_cells()), measured, ideal});
        std::cout << "  [sph] " << axis << " m=" << m << " iter=" << iterations
                  << " " << scheme_name(iterations, scheme) << " R=" << R << " r=" << r
                  << "  ->  " << 100.0 * measured / ideal << "%  ("
                  << mesh.n_cells() << " tets)\n";
    };
    auto add_dit = [&](const std::string& axis, int n, double R1, double R2,
                       double r, int iterations, int scheme) {
        auto mesh = build_ditorus(n, R1, R2, r, iterations, scheme);
        const double measured = Validation::Geometric::measured_volume(mesh);
        const double ideal = Validation::Geometric::ditorus_volume(R1, R2, r);
        rows.push_back({"ditorus", axis, scheme_name(iterations, scheme), n, iterations,
                        R1, R2, r, static_cast<long>(mesh.n_cells()), measured, ideal});
        std::cout << "  [dit] " << axis << " n=" << n << " iter=" << iterations
                  << " " << scheme_name(iterations, scheme) << " R1=" << R1 << " R2=" << R2
                  << " r=" << r << "  ->  " << 100.0 * measured / ideal << "%  ("
                  << mesh.n_cells() << " tets)\n";
    };

    // Spheritorus
    std::cout << "Spheritorus: segment refinement (no subdivision)\n";
    for (const auto& s : kSphSets)
        for (int m : {3, 4, 6, 8, 12, 16, 24, 32}) add_sph("segments", m, s.R, s.r, 0, 0);
    std::cout << "Spheritorus: subdivision refinement (minimal base m=3)\n";
    for (const auto& s : kSphSets) {
        const int base = 3;   // minimal valid circle resolution (ShapeGen::circle requires n >= 3)
        add_sph("subdivision", base, s.R, s.r, 0, 0);                    // shared baseline
        // Uniform 1->8: per-iteration blowup x8 stays tractable out to 5
        // iterations (up to ~2.4*10^6 tets here), far enough to show the
        // overshoot converging back towards 100%.
        for (int it : {1, 2, 3, 4, 5}) add_sph("subdivision", base, s.R, s.r, it, 1);
        // Barycentric 1->24: capped at 3 iterations (~10^6 tets). 5 iterations
        // would reach ~5.7*10^8 tets from this base -- infeasible, not just slow.
        for (int it : {1, 2, 3}) add_sph("subdivision", base, s.R, s.r, it, 0);
    }

    // Ditorus
    std::cout << "Ditorus: segment refinement (no subdivision)\n";
    for (const auto& s : kDitSets)
        for (int n : {3, 4, 6, 8, 10, 12}) add_dit("segments", n, s.R1, s.R2, s.r, 0, 0);
    std::cout << "Ditorus: subdivision refinement (minimal base n=3)\n";
    for (const auto& s : kDitSets) {
        const int base = 3;
        add_dit("subdivision", base, s.R1, s.R2, s.r, 0, 0);
        for (int it : {1, 2, 3, 4, 5}) add_dit("subdivision", base, s.R1, s.R2, s.r, it, 1);
        for (int it : {1, 2, 3}) add_dit("subdivision", base, s.R1, s.R2, s.r, it, 0);
    }

    // Write CSV
    const std::string dir = std::string(RESULTS_DIR) + "/schwarz_lantern";
    std::filesystem::create_directories(dir);
    const std::string path = dir + "/schwarz_lantern.csv";
    std::ofstream f(path);
    f << "shape,axis,segments,iterations,scheme,R1,R2,r,n_cells,measured,ideal,percent\n";
    for (const auto& row : rows) {
        f << row.shape << "," << row.axis << "," << row.seg << "," << row.iterations << ","
          << row.scheme << "," << row.R1 << "," << row.R2 << "," << row.r << ","
          << row.n_cells << "," << row.measured << "," << row.ideal << ","
          << (row.ideal != 0.0 ? 100.0 * row.measured / row.ideal : 0.0) << "\n";
    }
    std::cout << "\nWrote " << rows.size() << " rows to "
              << std::filesystem::absolute(path).string() << std::endl;
    return 0;
}
