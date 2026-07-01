#pragma once

/**
 * @file orientation.hh
 * @brief Consistent-orientation pass for closed tetrahedral 3-manifolds.
 */

#include "cartesian_product.hh"
#include <OpenVolumeMesh/Core/Handles.hh>
#include <algorithm>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace Orientation {

    // 1. PRECONDITION CHECK
    /**
     * @brief Verifies that the mesh contains only tetrahedra and triangles.
     *
     * @tparam MeshT  Any OpenVolumeMesh kernel type.
     * @param  mesh   The mesh to validate.
     * @throws std::runtime_error if any cell has != 4 vertices or any face
     *         has != 3 vertices.
     */
    template <typename MeshT>
    inline void verify_pure_tet_mesh(const MeshT& mesh) {
        int non_tet_cells = 0;
        int non_tri_faces = 0;

        for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it) {
            int n = 0;
            for (auto cv_it = mesh.cv_iter(*c_it); cv_it.valid(); ++cv_it) ++n;
            if (n != 4) ++non_tet_cells;
        }

        for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
            int n = 0;
            for (auto fv_it = mesh.fv_iter(*f_it); fv_it.valid(); ++fv_it) ++n;
            if (n != 3) ++non_tri_faces;
        }

        if (non_tet_cells > 0 || non_tri_faces > 0) {
            throw std::runtime_error(
                "Orientation pass requires a pure-tetrahedral mesh. Found "
                + std::to_string(non_tet_cells) + " non-tet cell(s) and "
                + std::to_string(non_tri_faces) + " non-triangle face(s). "
                "Run Subdivision::triangulate_prisms_to_tets first.");
        }
    }



    // 2. BFS ORIENTATION

    /**
     * @brief Computes a consistent orientation bit per cell via BFS.
     *
     * @tparam MeshT  Any OpenVolumeMesh kernel type.
     * @param  mesh   A pure-tetrahedral mesh built by the pipeline.
     * @return A vector of {+1, -1} bits indexed by cell handle.
     * @throws std::runtime_error if the mesh is non-orientable, contains
     *         non-tetrahedral cells, or has unreachable cells.
     */
    template <typename MeshT>
    inline std::vector<int> compute_orientation_bits(const MeshT& mesh) {
        using VH  = OpenVolumeMesh::VertexHandle;
        using FH  = OpenVolumeMesh::FaceHandle;
        using HFH = OpenVolumeMesh::HalfFaceHandle;
        using CH  = OpenVolumeMesh::CellHandle;

        std::cout << "Computing orientation bits via BFS..." << std::endl;

        verify_pure_tet_mesh(mesh);

        // Initialize a state-tracking buffer for the BFS to prevent loops,
        // and abort to prevent crashes if the mesh is topologically empty.
        const int n_cells = static_cast<int>(mesh.n_cells());
        std::vector<int> bits(n_cells, 0);   // 0 = unvisited

        if (n_cells == 0) return bits;

        // Helper: read the 4 vertices of a cell, sorted by vertex index.
        auto sorted_vertices = [&](CH c) {
            std::vector<VH> v;
            for (auto cv_it = mesh.cv_iter(c); cv_it.valid(); ++cv_it) v.push_back(*cv_it);
            std::sort(v.begin(), v.end(),
                      [](VH a, VH b){ return a.idx() < b.idx(); });
            return v;
        };

        // Helper: among the 4 sorted vertices of a cell, return the index
        // (0..3) of the one absent from the given face. This is the i in
        // (-1)^i
        auto missing_vertex_index = [&](const std::vector<VH>& sorted_v, FH f) {
            std::set<int> face_v;
            for (auto fv_it = mesh.fv_iter(f); fv_it.valid(); ++fv_it)
                face_v.insert(fv_it->idx());
            for (int k = 0; k < 4; ++k) {
                if (face_v.find(sorted_v[k].idx()) == face_v.end()) return k;
            }
            throw std::runtime_error(
                "BFS: no vertex of this cell is absent from the given face "
                "(corrupt cell/face incidence)");
        };

        // lock the first available cell as the positive
        // global anchor (1) and queue it for propagation
        CH seed = *mesh.cells_begin();
        bits[seed.idx()] = 1;
        std::queue<CH> q;
        q.push(seed);

        long propagations = 0;

        while (!q.empty()) {
            CH c = q.front(); q.pop();

            const auto sv_c = sorted_vertices(c);

            //translating from 0-1-2 scheme to 0/1 scheme
            const int b_c_enc = (bits[c.idx()] == 1) ? 0 : 1;

            // Walk all 4 faces of the cell.
            for (auto cf_it = mesh.cf_iter(c); cf_it.valid(); ++cf_it) {
                const FH f = *cf_it;
                const int i_c = missing_vertex_index(sv_c, f);
                if (i_c < 0) {
                    throw std::runtime_error(
                        "BFS: failed to locate missing vertex (corrupt cell?)");
                }

                // Find the neighbour cell across this face by looking up
                // both halffaces' incident cells. Exactly one of them is c;
                // the other is the neighbour (or invalid for a boundary).
                const HFH hf0 = mesh.halfface_handle(f, 0);
                const HFH hf1 = mesh.halfface_handle(f, 1);
                const CH c0 = mesh.incident_cell(hf0);
                const CH c1 = mesh.incident_cell(hf1);

                CH neighbour;
                if (c0 == c)      neighbour = c1;
                else if (c1 == c) neighbour = c0;
                else {
                    throw std::runtime_error(
                        "BFS: cell is not on either side of one of its own faces");
                }

                if (!neighbour.is_valid()) {
                    // Boundary face: one halfface has no incident cell, so
                    // there is no neighbour to propagate to. Skipped here;
                    // closedness is checked separately in validation.hh.
                    continue;
                }

                const auto sv_n = sorted_vertices(neighbour);
                const int i_n = missing_vertex_index(sv_n, f);
                if (i_n < 0) {
                    throw std::runtime_error(
                        "BFS: failed to locate missing vertex on neighbour");
                }

                // Parity relation: b_n = (1 + i_c + i_n + b_c) mod 2.
                const int b_n_enc = (1 + i_c + i_n + b_c_enc) % 2;
                const int new_bit = (b_n_enc == 0) ? 1 : -1;

                if (bits[neighbour.idx()] == 0) {
                    bits[neighbour.idx()] = new_bit;
                    q.push(neighbour);
                    ++propagations;
                } else if (bits[neighbour.idx()] != new_bit) {
                    throw std::runtime_error(
                        "BFS orientation: inconsistency at cell "
                        + std::to_string(neighbour.idx())
                        + " (non-orientable mesh detected).");
                } else {
                    // Already assigned the correct bit; consistent, no action needed.
                }
            }
        }

        // Connectivity check: every cell must have received a bit
        int unvisited = 0;
        for (int i = 0; i < n_cells; ++i) {
            if (bits[i] == 0) ++unvisited;
        }
        if (unvisited > 0) {
            throw std::runtime_error(
                "BFS orientation: " + std::to_string(unvisited)
                + " cell(s) not reached from seed (mesh has multiple components).");
        }

        int positive = 0, negative = 0;
        for (int b : bits) {
            if (b == 1) ++positive;
            else if (b == -1) ++negative;
        }
        std::cout << "  -> " << positive << " positively-oriented, "
                  << negative << " negatively-oriented "
                  << "(over " << n_cells << " cells, "
                  << propagations << " propagations)." << std::endl;

        return bits;
    }


    // 3. REBUILD WITH CORRECT HALFFACE CHOICES
    /**
     * @brief Builds a new mesh whose halfface choices match the bit vector, because
     * OpenVolumeMesh cannot mutate a cell's halfface list after
     * construction.
     *
     * Each tet is re-emitted with halfface hf((i + b) mod 2) for the face
     * opposite sorted-vertex v_i.
     *
     * @tparam MeshT  Any OpenVolumeMesh kernel type.
     * @param  mesh   Source mesh (pure-tet, must match the bit vector).
     * @param  bits   Per-cell orientation bits from compute_orientation_bits.
     * @return A new mesh, topologically identical to the source but with
     *         consistently oriented halfface choices for every cell.
     */
    template <typename MeshT>
    inline MeshT rebuild_with_orientation(const MeshT& mesh,
                                          const std::vector<int>& bits) {
        using VH  = OpenVolumeMesh::VertexHandle;
        using FH  = OpenVolumeMesh::FaceHandle;
        using HFH = OpenVolumeMesh::HalfFaceHandle;

        std::cout << "Rebuilding mesh with oriented halffaces..." << std::endl;

        verify_pure_tet_mesh(mesh);

        if (static_cast<int>(bits.size()) != static_cast<int>(mesh.n_cells())) {
            throw std::runtime_error(
                "rebuild_with_orientation: bit vector size does not match "
                "cell count.");
        }

        MeshT new_mesh;

        // Carry over all vertices in their original handle order to
        // preserve the index correspondence
        std::vector<VH> v_old_to_new(mesh.n_vertices());
        for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
            v_old_to_new[v_it->idx()] = new_mesh.add_vertex(mesh.vertex(*v_it));
        }

        // Shared face dictionary across the whole rebuild
        std::map<std::vector<int>, FH> face_map;

        for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it) {
            const int bit = bits[c_it->idx()];
            if (bit != 1 && bit != -1) {
                throw std::runtime_error(
                    "rebuild_with_orientation: cell "
                    + std::to_string(c_it->idx())
                    + " has invalid orientation bit");
            }
            const int b = (bit == 1) ? 0 : 1;

            // Sorted-vertex tuple of the cell, translated into new_mesh's
            // index space.
            std::vector<VH> verts;
            for (auto cv_it = mesh.cv_iter(*c_it); cv_it.valid(); ++cv_it) {
                verts.push_back(v_old_to_new[cv_it->idx()]);
            }

            std::sort(verts.begin(), verts.end(),
                      [](VH a, VH b){ return a.idx() < b.idx(); });

            const VH v0 = verts[0], v1 = verts[1], v2 = verts[2], v3 = verts[3];

            // The four boundary triangles, indexed by the sorted-vertex
            // they omit. get_or_create_face will sort each triple
            // internally so hf0 ends up as the lexicographic cycle.
            const FH f0 = CartesianOps::detail::get_or_create_face<MeshT, VH, FH>(
                              new_mesh, {v1, v2, v3}, face_map);   // omits v_0, sign +
            const FH f1 = CartesianOps::detail::get_or_create_face<MeshT, VH, FH>(
                              new_mesh, {v0, v2, v3}, face_map);   // omits v_1, sign -
            const FH f2 = CartesianOps::detail::get_or_create_face<MeshT, VH, FH>(
                              new_mesh, {v0, v1, v3}, face_map);   // omits v_2, sign +
            const FH f3 = CartesianOps::detail::get_or_create_face<MeshT, VH, FH>(
                              new_mesh, {v0, v1, v2}, face_map);   // omits v_3, sign -

            // hf((i + b) mod 2) for face opposite v_i.
            std::vector<HFH> hfs = {
                new_mesh.halfface_handle(f0, (0 + b) % 2),
                new_mesh.halfface_handle(f1, (1 + b) % 2),
                new_mesh.halfface_handle(f2, (0 + b) % 2),   // (2+b)%2 = b
                new_mesh.halfface_handle(f3, (1 + b) % 2),   // (3+b)%2 = (1+b)%2
            };

            new_mesh.add_cell(hfs);
        }

        std::cout << "  -> Rebuilt: V=" << new_mesh.n_vertices()
                  << " E=" << new_mesh.n_edges()
                  << " F=" << new_mesh.n_faces()
                  << " C=" << new_mesh.n_cells() << std::endl;

        return new_mesh;
    }


    // 4. CONVENIENCE: ORIENT IN PLACE
    /**
     * @brief Computes the orientation bits and replaces mesh with the
     *        rebuilt, consistently oriented version.
     *
     * @tparam MeshT  Any OpenVolumeMesh kernel type.
     * @param  mesh   Input mesh; replaced in place by its oriented rebuild.
     */
    template <typename MeshT>
    inline void  orient_mesh(MeshT& mesh) {
        const auto bits = compute_orientation_bits(mesh);
        MeshT oriented = rebuild_with_orientation(mesh, bits);
        mesh = oriented;  // copy-assign from lvalue; move-assign is deleted in GeometryKernel
    }


    // 4b. FLIP ORIENTATION

    /**
     * @brief Globally flips the orientation of every cell in the mesh.
     *
     * @tparam MeshT  Any OpenVolumeMesh kernel type.
     * @param  mesh   Input mesh; replaced in place by its orientation-flipped rebuild.
     */
    template <typename MeshT>
    inline void flip_orientation(MeshT& mesh) {
        auto bits = compute_orientation_bits(mesh);
        for (auto& b : bits) b = -b;
        MeshT flipped = rebuild_with_orientation(mesh, bits);
        mesh = flipped;
    }
} // namespace Orientation
