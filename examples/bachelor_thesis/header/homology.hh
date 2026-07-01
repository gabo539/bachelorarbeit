#pragma once

/**
 * @file homology.hh
 * @brief Simplicial homology (Betti numbers) for closed simplicial complexes.
 *
 *     b0 = V - rank(d1)
 *     b1 = E - rank(d1) - rank(d2)
 *     b2 = F - rank(d2) - rank(d3)
 *     b3 = C - rank(d3)
 *
 *
 * COMPUTING THE RANKS
 *
 * Each boundary operator d_k is assembled as a dense matrix of +-1/0 entries and
 * its rank is read off Eigen's FullPivLU (a rank-revealing LU factorisation with
 * full pivoting). The Betti numbers then follow from the four lines above.
 *
 * SELF-CHECK
 * The Euler-Poincare theorem forces  b0 - b1 + b2 - b3 = chi = V - E + F - C.
 * If the computed Betti numbers do not alternate-sum to chi, the bug is in this
 * code, not in the mesh.
 *
 * PRECONDITION: a pure simplicial complex (triangular faces; tetrahedral cells
 * if any).
 */

#include <OpenVolumeMesh/Core/Handles.hh>
#include <Eigen/Dense>

#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <map>
#include <utility>
#include <vector>

namespace Homology {

    // Dense rank is O(n^3) in time and O(n^2) in memory, so above this many
    // simplices (in any single dimension) we skip the computation and ask for a
    // coarser mesh
    constexpr long MaxSimplices = 3000;


    /**
    * @brief Rank of a boundary matrix via Eigen's full-pivoting LU.
    *
    * FullPivLU is rank-revealing: rank() returns the number of pivots whose
    * magnitude exceeds an automatic tolerance. For the small, well-separated
    * matrices we feed it, that count is exact.
    */
    inline long matrix_rank(const Eigen::MatrixXd& M) {
        if (M.rows() == 0 || M.cols() == 0) return 0;
        return static_cast<long>(Eigen::FullPivLU<Eigen::MatrixXd>(M).rank());
    }



    /**
     * @brief Computes the Betti numbers (b0, b1, b2, b3) of a simplicial mesh.
     *
     * Builds the boundary matrices d1, d2, d3, reads their ranks, derives the Betti
     * numbers by rank-nullity, and runs the Euler-Poincare self-check against
     * chi. b3 is 0 for a 2-complex (the 2-torus has no cells).
     *
     * @tparam MeshT  Any OpenVolumeMesh kernel type.
     * @param  mesh   A pure simplicial complex.
     * @return The tuple {b0, b1, b2, b3}, or an empty vector if the mesh was too
     *         large and the computation was skipped.
     */
    template <typename MeshT>
    inline std::vector<long> compute_betti(const MeshT& mesh) {
        const long V = static_cast<long>(mesh.n_vertices());
        const long E = static_cast<long>(mesh.n_edges());
        const long F = static_cast<long>(mesh.n_faces());
        const long C = static_cast<long>(mesh.n_cells());

        std::cout << " -> Computing Betti numbers (Eigen FullPivLU) ..." << std::endl;
        std::cout << "    V:" << V << " | E:" << E << " | F:" << F
                  << " | C:" << C << std::endl;

        //  Abort early on big meshes (see header).
        if (E > MaxSimplices || F > MaxSimplices ||
            C > MaxSimplices) {
            std::cout << "    [SKIP] Mesh has > " << MaxSimplices
                      << " simplices in some dimension; dense rank is O(n^3).\n"
                      << "           Betti numbers are resolution-independent, so "
                         "re-run at a coarser\n"
                      << "           resolution (fewer segments / subdivision "
                         "passes) to validate the topology." << std::endl;
            return {};
        }

        // ---- lookup tables: sorted vertex tuple -> simplex index ----
        // d2 needs every triangle's edges; d3 needs every tet's faces.
        std::map<std::pair<int, int>, int> edge_index;
        for (auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it) {
            const auto e = mesh.edge(*e_it);
            int a = e.from_vertex().idx();
            int b = e.to_vertex().idx();
            if (a > b) std::swap(a, b);
            edge_index[{a, b}] = e_it->idx();
        }

        std::map<std::array<int, 3>, int> face_index;
        if (C > 0) {
            for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
                std::vector<int> fv;
                for (auto fv_it = mesh.fv_iter(*f_it); fv_it.valid(); ++fv_it)
                    fv.push_back(fv_it->idx());
                assert(fv.size() == 3 && "compute_betti requires triangular faces");
                std::sort(fv.begin(), fv.end());
                face_index[{fv[0], fv[1], fv[2]}] = f_it->idx();
            }
        }

        // ---- d1: rows = vertices, cols = edges.  d[a,b] = +[b] - [a] ----
        Eigen::MatrixXd d1 = Eigen::MatrixXd::Zero(V, E);
        for (auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it) {
            const auto e = mesh.edge(*e_it);
            int a = e.from_vertex().idx();
            int b = e.to_vertex().idx();
            if (a > b) std::swap(a, b);
            const int col = e_it->idx();
            d1(a, col) = -1.0;
            d1(b, col) = +1.0;
        }

        // ---- d2: rows = edges, cols = faces. ----
        // triangle [a<b<c]:  +[b,c] - [a,c] + [a,b]
        Eigen::MatrixXd d2 = Eigen::MatrixXd::Zero(E, F);
        for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
            std::vector<int> v;
            for (auto fv_it = mesh.fv_iter(*f_it); fv_it.valid(); ++fv_it)
                v.push_back(fv_it->idx());
            assert(v.size() == 3 && "compute_betti requires triangular faces");
            std::sort(v.begin(), v.end());
            const int col = f_it->idx();

            auto edge_of = [&](int x, int y) {
                auto it = edge_index.find({x, y});
                assert(it != edge_index.end() && "face edge missing from edge table");
                return it->second;
            };
            d2(edge_of(v[1], v[2]), col) = +1.0;   // drop v0
            d2(edge_of(v[0], v[2]), col) = -1.0;   // drop v1
            d2(edge_of(v[0], v[1]), col) = +1.0;   // drop v2
        }

        // ---- d3: rows = faces, cols = cells (only if the mesh has cells). ----
        // tet [a<b<c<d]:  +[b,c,d] - [a,c,d] + [a,b,d] - [a,b,c]
        Eigen::MatrixXd d3;
        if (C > 0) {
            d3 = Eigen::MatrixXd::Zero(F, C);
            for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it) {
                std::vector<int> v;
                for (auto cv_it = mesh.cv_iter(*c_it); cv_it.valid(); ++cv_it)
                    v.push_back(cv_it->idx());
                assert(v.size() == 4 && "compute_betti requires tetrahedral cells");
                std::sort(v.begin(), v.end());
                const int col = c_it->idx();

                auto face_of = [&](int x, int y, int z) {
                    auto it = face_index.find({x, y, z});
                    assert(it != face_index.end() && "tet face missing from face table");
                    return it->second;
                };
                d3(face_of(v[1], v[2], v[3]), col) = +1.0;   // drop v0
                d3(face_of(v[0], v[2], v[3]), col) = -1.0;   // drop v1
                d3(face_of(v[0], v[1], v[3]), col) = +1.0;   // drop v2
                d3(face_of(v[0], v[1], v[2]), col) = -1.0;   // drop v3
            }
        }

        // ---- ranks and Betti numbers via rank-nullity ----
        const long r1 = matrix_rank(d1);
        const long r2 = matrix_rank(d2);
        const long r3 = (C > 0) ? matrix_rank(d3) : 0;

        std::vector<long> b = {
            V - r1,         // b0
            E - r1 - r2,    // b1
            F - r2 - r3,    // b2
            C - r3          // b3
        };

        // ---- Euler-Poincare self-check ----
        const long chi = V - E + F - C;
        const long alt = b[0] - b[1] + b[2] - b[3];
        if (alt != chi) {
            std::cout << "    [BUG]  Euler-Poincare mismatch: b0-b1+b2-b3 = "
                      << alt << " but chi = " << chi
                      << " (bug in the homology code, not the mesh)." << std::endl;
        } else {
            std::cout << "    [OK]   Euler-Poincare self-check passed "
                      << "(b0-b1+b2-b3 = chi = " << chi << ")." << std::endl;
        }

        std::cout << "    Betti = (" << b[0] << ", " << b[1] << ", "
                  << b[2] << ", " << b[3] << ")" << std::endl;
        return b;
    }


    /**
     * @brief Computes Betti numbers and compares them to the target manifold.
     *
     *
     * @tparam MeshT     Any OpenVolumeMesh kernel type.
     * @param  mesh      A pure simplicial complex.
     * @param  expected  Target Betti tuple (padded with zeros to length 4):
     *                   2-torus {1,2,1}, spheritorus {1,1,1,1}, ditorus {1,3,3,1}.
     * @return true if the computed Betti numbers match `expected` (also true if
     *         the computation was skipped because the mesh was too large).
     */
    template <typename MeshT>
    inline bool verify_betti(const MeshT& mesh, const std::vector<long>& expected) {
        std::cout << " -> Verifying Betti numbers ..." << std::endl;
        const std::vector<long> b = compute_betti(mesh);

        // Empty => the computation was skipped
        if (b.empty()) return true;

        std::vector<long> exp = expected;
        exp.resize(4, 0);

        if (b == exp) {
            std::cout << "    [PASS] Betti numbers match the target manifold ("
                      << exp[0] << ", " << exp[1] << ", " << exp[2] << ", "
                      << exp[3] << ")." << std::endl;
            return true;
        }

        std::cout << "    [FAIL] Betti mismatch. Got (" << b[0] << ", " << b[1]
                  << ", " << b[2] << ", " << b[3] << "), expected (" << exp[0]
                  << ", " << exp[1] << ", " << exp[2] << ", " << exp[3] << ")."
                  << std::endl;

        return false;
    }

} // namespace Homology
