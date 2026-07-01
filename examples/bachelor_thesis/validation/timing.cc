/**
 * @file timing.cc
 * @brief Runtime benchmark for the spheritorus and ditorus pipelines.
 *
 * For each configuration the construction and validation are split into stages
 * and timed separately, on the meshes the generators actually run them on:
 *   construction : product, prism triangulation, subdivision, projection, orientation
 *   validation   : verify_topology + verify_betti (coarse mesh),
 *                  closed/degenerate/on-surface checks (final mesh)
 * Each configuration is built once as warm-up and then `reps` times; per stage
 * the median over the timed runs is reported, plus the minimum total build time.
 *
 * Build Release with -O2 -DNDEBUG (asserts off). Writes results/timing/timing.csv.
 */

#include <algorithm>
#include <chrono>
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
#include "orientation.hh"
#include "validation.hh"
#include "homology.hh"
#include "projection.hh"
#include "results_dir.hh"

namespace {

using clk = std::chrono::steady_clock;
inline double ms(clk::time_point a, clk::time_point b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
}

// Silences the pipelines progress prints while one mesh is built and checked.
struct CoutMuter {
    std::streambuf* old_buf;
    std::ostringstream sink;
    CoutMuter() : old_buf(std::cout.rdbuf(sink.rdbuf())) {}
    ~CoutMuter() { std::cout.rdbuf(old_buf); }
};

struct Stages {
    double product = 0, triangulate = 0, subdivide = 0, project = 0, orient = 0;
    double verify_topo = 0, verify_betti = 0, verify_geom = 0;
    long   n_coarse = 0, n_final = 0;
};

Stages run_spheritorus_once(int m, double R, double r, int iterations, int scheme) {
    CoutMuter mute;
    Stages s;
    auto s2 = ShapeGen::octahedron(1.0);
    auto s1 = ShapeGen::circle(m);

    auto t0 = clk::now();
    Geometry::Mesh5D mesh = CartesianOps::cartesian_product<Geometry::Mesh5D>(
        s2, s1, [](const Geometry::Vec4d& a, const Geometry::Vec4d& b) {
            return Geometry::Vec5d(a[0], a[1], a[2], b[0], b[1]); });
    auto t1 = clk::now();
    Subdivision::triangulate_prisms_to_tets(mesh);
    auto t2 = clk::now();
    s.n_coarse = static_cast<long>(mesh.n_cells());

    Validation::Topological::verify_topology(mesh);
    auto t3 = clk::now();
    Homology::verify_betti(mesh, {1, 1, 1, 1});
    auto t4 = clk::now();

    for (int i = 0; i < iterations; ++i) {
        if (scheme == 1) Subdivision::uniform_subdivide(mesh);
        else             Subdivision::barycentric_subdivide(mesh);
        Subdivision::normalize_blocks(mesh, {3, 2});
    }
    auto t5 = clk::now();
    Geometry::Mesh4D shape = Projection::project_spheritorus(mesh, R, r);
    auto t6 = clk::now();
    Orientation::orient_mesh(shape);
    auto t7 = clk::now();
    s.n_final = static_cast<long>(shape.n_cells());

    Validation::Topological::verify_closed_manifold(shape);
    Validation::Geometric::verify_no_degenerate_cells(shape);
    Validation::Geometric::verify_on_spheritorus(shape, R, r);
    auto t8 = clk::now();

    s.product = ms(t0, t1); s.triangulate = ms(t1, t2);
    s.verify_topo = ms(t2, t3); s.verify_betti = ms(t3, t4);
    s.subdivide = ms(t4, t5); s.project = ms(t5, t6); s.orient = ms(t6, t7);
    s.verify_geom = ms(t7, t8);
    return s;
}

Stages run_ditorus_once(int n, double R1, double R2, double r, int iterations, int scheme) {
    CoutMuter mute;
    Stages s;
    auto s1_in  = ShapeGen::circle(n);
    auto s1_mid = ShapeGen::circle(n);
    auto s1_out = ShapeGen::circle(n);

    auto t0 = clk::now();
    Geometry::Mesh4D torus = CartesianOps::cartesian_product<Geometry::Mesh4D>(
        s1_in, s1_mid, [](const Geometry::Vec4d& a, const Geometry::Vec4d& b) {
            return Geometry::Vec4d(a[0], a[1], b[0], b[1]); });
    Subdivision::triangulate_quad_faces(torus);
    Geometry::Mesh6D mesh = CartesianOps::cartesian_product<Geometry::Mesh6D>(
        torus, s1_out, [](const Geometry::Vec4d& ab, const Geometry::Vec4d& c) {
            return Geometry::Vec6d(ab[0], ab[1], ab[2], ab[3], c[0], c[1]); });
    auto t1 = clk::now();   // product = inner-torus product + quad triangulation + outer product
    Subdivision::triangulate_prisms_to_tets(mesh);
    auto t2 = clk::now();
    s.n_coarse = static_cast<long>(mesh.n_cells());

    Validation::Topological::verify_topology(mesh);
    auto t3 = clk::now();
    Homology::verify_betti(mesh, {1, 3, 3, 1});
    auto t4 = clk::now();

    for (int i = 0; i < iterations; ++i) {
        if (scheme == 1) Subdivision::uniform_subdivide(mesh);
        else             Subdivision::barycentric_subdivide(mesh);
        Subdivision::normalize_blocks(mesh, {2, 2, 2});
    }
    auto t5 = clk::now();
    Geometry::Mesh4D shape = Projection::project_ditorus(mesh, R1, R2, r);
    auto t6 = clk::now();
    Orientation::orient_mesh(shape);
    auto t7 = clk::now();
    s.n_final = static_cast<long>(shape.n_cells());

    Validation::Topological::verify_closed_manifold(shape);
    Validation::Geometric::verify_no_degenerate_cells(shape);
    Validation::Geometric::verify_on_ditorus(shape, R1, R2, r);
    auto t8 = clk::now();

    s.product = ms(t0, t1); s.triangulate = ms(t1, t2);
    s.verify_topo = ms(t2, t3); s.verify_betti = ms(t3, t4);
    s.subdivide = ms(t4, t5); s.project = ms(t5, t6); s.orient = ms(t6, t7);
    s.verify_geom = ms(t7, t8);
    return s;
}

double median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

struct Row {
    std::string shape, scheme;
    int seg, iterations;
    double R1, R2, r;
    long n_coarse, n_final;
    double product, triangulate, subdivide, project, orient,
           total_build, total_build_min,
           verify_topo, verify_betti, verify_geom, total_verify;
};

template <typename RunOnce>
Row aggregate(const std::string& shape, const std::string& scheme,
              int seg, int iterations, double R1, double R2, double r,
              int reps, RunOnce run_once) {
    run_once();                                   // warm-up (untimed)
    std::vector<Stages> runs;
    for (int i = 0; i < reps; ++i) runs.push_back(run_once());

    auto col = [&](double Stages::*field) {
        std::vector<double> v;
        for (const auto& s : runs) v.push_back(s.*field);
        return median(v);
    };
    std::vector<double> builds, verifies;
    for (const auto& s : runs) {
        builds.push_back(s.product + s.triangulate + s.subdivide + s.project + s.orient);
        verifies.push_back(s.verify_topo + s.verify_betti + s.verify_geom);
    }

    Row row;
    row.shape = shape; row.scheme = scheme; row.seg = seg; row.iterations = iterations;
    row.R1 = R1; row.R2 = R2; row.r = r;
    row.n_coarse = runs.front().n_coarse; row.n_final = runs.front().n_final;
    row.product = col(&Stages::product);
    row.triangulate = col(&Stages::triangulate);
    row.subdivide = col(&Stages::subdivide);
    row.project = col(&Stages::project);
    row.orient = col(&Stages::orient);
    row.verify_topo = col(&Stages::verify_topo);
    row.verify_betti = col(&Stages::verify_betti);
    row.verify_geom = col(&Stages::verify_geom);
    row.total_build = median(builds);
    row.total_build_min = *std::min_element(builds.begin(), builds.end());
    row.total_verify = median(verifies);

    std::cout << "  [" << shape << "] seg=" << seg << " iter=" << iterations << " "
              << scheme << "  n=" << row.n_final << "  build " << row.total_build
              << " ms (min " << row.total_build_min << "), verify "
              << row.total_verify << " ms\n";
    return row;
}

} // namespace

int main() {
    const int reps = 7;
    std::vector<Row> rows;

    // Radii are fixed: timing is combinatorial and does not depend on them.
    const double SR = 5.0, Sr = 1.0;                    // spheritorus
    const double DR1 = 10.0, DR2 = 5.0, Dr = 1.0;      // ditorus

    std::cout << "Spheritorus: segment count variation (no subdivision)\n";
    for (int m : {4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128})
        rows.push_back(aggregate("spheritorus", "none", m, 0, SR, 0, Sr, reps,
            [&]{ return run_spheritorus_once(m, SR, Sr, 0, 0); }));

    std::cout << "Spheritorus: subdivision iteration count variation (base m=8)\n";
    {
        const int base = 8;
        rows.push_back(aggregate("spheritorus", "none", base, 0, SR, 0, Sr, reps,
            [&]{ return run_spheritorus_once(base, SR, Sr, 0, 0); }));
        for (int p : {1, 2, 3})
            rows.push_back(aggregate("spheritorus", "uniform", base, p, SR, 0, Sr, reps,
                [&]{ return run_spheritorus_once(base, SR, Sr, p, 1); }));
        for (int p : {1, 2})
            rows.push_back(aggregate("spheritorus", "barycentric", base, p, SR, 0, Sr, reps,
                [&]{ return run_spheritorus_once(base, SR, Sr, p, 0); }));
    }

    std::cout << "Ditorus: segment count variation (no subdivision)\n";
    for (int n : {3, 4, 6, 8, 10, 12, 16, 20})
        rows.push_back(aggregate("ditorus", "none", n, 0, DR1, DR2, Dr, reps,
            [&]{ return run_ditorus_once(n, DR1, DR2, Dr, 0, 0); }));

    std::cout << "Ditorus: subdivision iteration count variation (base n=4)\n";
    {
        const int base = 4;
        rows.push_back(aggregate("ditorus", "none", base, 0, DR1, DR2, Dr, reps,
            [&]{ return run_ditorus_once(base, DR1, DR2, Dr, 0, 0); }));
        for (int p : {1, 2, 3})
            rows.push_back(aggregate("ditorus", "uniform", base, p, DR1, DR2, Dr, reps,
                [&]{ return run_ditorus_once(base, DR1, DR2, Dr, p, 1); }));
        for (int p : {1, 2})
            rows.push_back(aggregate("ditorus", "barycentric", base, p, DR1, DR2, Dr, reps,
                [&]{ return run_ditorus_once(base, DR1, DR2, Dr, p, 0); }));
    }

    const std::string dir = std::string(RESULTS_DIR) + "/timing";
    std::filesystem::create_directories(dir);
    const std::string path = dir + "/timing.csv";
    std::ofstream f(path);
    f << "shape,R1,R2,r,segments,iterations,scheme,n_coarse,n_final,"
         "product_ms,triangulate_ms,subdivide_ms,project_ms,orient_ms,"
         "total_build_ms,total_build_min_ms,"
         "verify_topo_ms,verify_betti_ms,verify_geom_ms,total_verify_ms\n";
    for (const auto& r : rows) {
        f << r.shape << "," << r.R1 << "," << r.R2 << "," << r.r << ","
          << r.seg << "," << r.iterations << "," << r.scheme << ","
          << r.n_coarse << "," << r.n_final << ","
          << r.product << "," << r.triangulate << "," << r.subdivide << ","
          << r.project << "," << r.orient << ","
          << r.total_build << "," << r.total_build_min << ","
          << r.verify_topo << "," << r.verify_betti << "," << r.verify_geom << ","
          << r.total_verify << "\n";
    }
    std::cout << "\nWrote " << rows.size() << " rows to "
              << std::filesystem::absolute(path).string() << std::endl;
    return 0;
}
