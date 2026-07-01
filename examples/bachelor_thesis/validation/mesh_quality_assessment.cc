/**
 * @file mesh_quality_assessment.cc
 * @brief Mesh-quality assessment for the spheritorus and the ditorus across
 *        three independent experiments:
 *
 *   1. Subdivision-method comparison. A fixed coarse base is refined by an
 *      increasing number of subdivision iterations, under each of the two
 *      schemes (uniform 1->8, barycentric 1->24), and compared against the
 *      unrefined mesh.
 *   2. Base-resolution refinement. The base n-gon resolution is increased
 *      directly, with no subdivision applied at all.
 *   3. Radius-ratio invariance. At a fixed, moderate resolution, the radii
 *      defining the embedding (Section 2.7) are varied.
 *
 * Per-configuration distributions are written under results/quality/, in
 * addition to the compact summary table written here.
 *
 * Quality metrics do not depend on the orientation , therefore
 * the orientation algorithm is not applied here.
 *
 * Writes results/quality/mesh_quality_assessment.csv (the compact summary
 * table) and results/quality/<configuration>_{...}.csv per configuration.
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
#include "quality_metrics.hh"
#include "results_dir.hh"

namespace {

// Silences the pipelines progress prints while one mesh is built.
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

// Builds a projected ditorus (R^4). All three circles share the segment count n.
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

const char* scheme_name(int iterations, int scheme) {
    if (iterations == 0) return "none";
    return (scheme == 1) ? "uniform" : "barycentric";
}

struct Row {
    std::string shape, strategy, scheme;
    int segments = 0, iterations = 0;
    double R1 = 0.0, R2 = 0.0, r = 0.0;
    long n_cells = 0;
    double volume_percent = 0.0;
    QualityMetrics::QualitySummary q;
};

constexpr double kPoorMeanRatioThreshold = 0.1;

} // namespace

int main() {
    std::vector<Row> rows;
    const std::string quality_dir = std::string(RESULTS_DIR) + "/quality";

    auto record_spheritorus = [&](const std::string& strategy, int m, double R, double r,
                                  int iterations, int scheme) {
        auto mesh = build_spheritorus(m, R, r, iterations, scheme);
        const std::string sname = scheme_name(iterations, scheme);
        const std::string prefix = quality_dir + "/spheritorus_" + strategy + "_" + sname
            + "_m" + std::to_string(m) + "_it" + std::to_string(iterations)
            + "_R" + std::to_string(static_cast<int>(R)) + "_r" + std::to_string(static_cast<int>(r));
        auto qs = QualityMetrics::analyze_and_export(mesh, prefix, kPoorMeanRatioThreshold);

        const double ideal = Validation::Geometric::spheritorus_volume(R, r);
        const double measured = Validation::Geometric::measured_volume(mesh);

        Row row;
        row.shape = "spheritorus"; row.strategy = strategy; row.scheme = sname;
        row.segments = m; row.iterations = iterations;
        row.R1 = R; row.R2 = 0.0; row.r = r;
        row.n_cells = static_cast<long>(mesh.n_cells());
        row.volume_percent = (ideal != 0.0) ? 100.0 * measured / ideal : 0.0;
        row.q = qs;
        rows.push_back(row);

        std::cout << "  [sph] " << strategy << " m=" << m << " it=" << iterations
                  << " " << sname << " R=" << R << " r=" << r
                  << "  -> eta(median)=" << qs.mean_ratio.median
                  << " eta(min)=" << qs.mean_ratio.min
                  << " poor=" << (100.0 * qs.poor_element_fraction) << "%"
                  << "  (" << row.n_cells << " tets)\n";
    };

    auto record_ditorus = [&](const std::string& strategy, int n, double R1, double R2,
                              double r, int iterations, int scheme) {
        auto mesh = build_ditorus(n, R1, R2, r, iterations, scheme);
        const std::string sname = scheme_name(iterations, scheme);
        const std::string prefix = quality_dir + "/ditorus_" + strategy + "_" + sname
            + "_n" + std::to_string(n) + "_it" + std::to_string(iterations)
            + "_R1" + std::to_string(static_cast<int>(R1))
            + "_R2" + std::to_string(static_cast<int>(R2))
            + "_r" + std::to_string(static_cast<int>(r));
        auto qs = QualityMetrics::analyze_and_export(mesh, prefix, kPoorMeanRatioThreshold);

        const double ideal = Validation::Geometric::ditorus_volume(R1, R2, r);
        const double measured = Validation::Geometric::measured_volume(mesh);

        Row row;
        row.shape = "ditorus"; row.strategy = strategy; row.scheme = sname;
        row.segments = n; row.iterations = iterations;
        row.R1 = R1; row.R2 = R2; row.r = r;
        row.n_cells = static_cast<long>(mesh.n_cells());
        row.volume_percent = (ideal != 0.0) ? 100.0 * measured / ideal : 0.0;
        row.q = qs;
        rows.push_back(row);

        std::cout << "  [dit] " << strategy << " n=" << n << " it=" << iterations
                  << " " << sname << " R1=" << R1 << " R2=" << R2 << " r=" << r
                  << "  -> eta(median)=" << qs.mean_ratio.median
                  << " eta(min)=" << qs.mean_ratio.min
                  << " poor=" << (100.0 * qs.poor_element_fraction) << "%"
                  << "  (" << row.n_cells << " tets)\n";
    };


    // Experiment 1: refinement-strategy comparison (uniform vs. barycentric),
    // fixed coarse base, increasing iteration count.
    std::cout << "Experiment 1: refinement-strategy comparison\n";
    {
        const double R = 5.0, r = 1.0;
        const int base_m = 4;
        record_spheritorus("subdivision", base_m, R, r, 0, 0);
        for (int it : {1, 2, 3}) record_spheritorus("subdivision", base_m, R, r, it, 1);
        for (int it : {1, 2})    record_spheritorus("subdivision", base_m, R, r, it, 0);
    }
    {
        const double R1 = 10.0, R2 = 5.0, r = 1.0;
        const int base_n = 4;
        record_ditorus("subdivision", base_n, R1, R2, r, 0, 0);
        for (int it : {1, 2, 3}) record_ditorus("subdivision", base_n, R1, R2, r, it, 1);
        for (int it : {1, 2})    record_ditorus("subdivision", base_n, R1, R2, r, it, 0);
    }


    // Experiment 2: base-resolution refinement, no subdivision. For the ditorus
    // all three circles share the segment count, so this refines the mesh
    // uniformly in every direction.
    // For the spheritorus only the S^1 factor is parametrised this way (the
    // octahedral S^2 factor is refined only by subdivision)
    std::cout << "Experiment 2: base-resolution refinement (no subdivision)\n";
    {
        const double R = 5.0, r = 1.0;
        for (int m : {4, 8, 16, 32, 64}) record_spheritorus("base_resolution", m, R, r, 0, 0);
    }
    {
        const double R1 = 10.0, R2 = 5.0, r = 1.0;
        for (int n : {3, 4, 6, 8, 12}) record_ditorus("base_resolution", n, R1, R2, r, 0, 0);
    }

    // -------------------------------------------------------------------
    // Experiment 3: radius-ratio invariance. Fixed, moderate resolution; the
    // radii are varied while the mesh combinatorics stay constant.
    // -------------------------------------------------------------------
    std::cout << "Experiment 3: radius-ratio invariance\n";
    {
        const int m = 8;
        for (double r : {1.0, 2.0, 3.0, 4.0})   // R=5 fixed; r < R required
            record_spheritorus("radius_ratio", m, 5.0, r, 0, 0);
    }
    {
        const int n = 6;
        // Tube-to-middle-radius ratio: R1, R2 fixed, r varied (r < R2 required).
        for (double r : {1.0, 2.0, 3.0, 4.0})
            record_ditorus("radius_ratio", n, 10.0, 5.0, r, 0, 0);
        // Middle-to-outer-radius ratio: R1, r fixed, R2 varied
        // (R2 > r and R1 > R2 + r required).
        for (double R2 : {2.0, 3.0, 5.0, 7.0})
            record_ditorus("radius_ratio", n, 10.0, R2, 1.0, 0, 0);
    }


    // Write the summary table.
    std::filesystem::create_directories(quality_dir);
    const std::string path = quality_dir + "/mesh_quality_assessment.csv";
    std::ofstream f(path);
    f << "shape,strategy,scheme,segments,iterations,R1,R2,r,n_cells,"
         "volume_percent,"
         "mean_ratio_min,mean_ratio_median,mean_ratio_mean,"
         "dihedral_min_deg,dihedral_max_deg,"
         "vertex_degree_mean,poor_element_fraction\n";
    for (const auto& row : rows) {
        f << row.shape << "," << row.strategy << "," << row.scheme << ","
          << row.segments << "," << row.iterations << ","
          << row.R1 << "," << row.R2 << "," << row.r << "," << row.n_cells << ","
          << row.volume_percent << ","
          << row.q.mean_ratio.min << "," << row.q.mean_ratio.median << ","
          << row.q.mean_ratio.mean << ","
          << row.q.dihedral_deg.min << "," << row.q.dihedral_deg.max << ","
          << row.q.vertex_degree.mean << "," << row.q.poor_element_fraction << "\n";
    }
    std::cout << "\nWrote " << rows.size() << " configurations to "
              << std::filesystem::absolute(path).string() << std::endl;

    return 0;
}
