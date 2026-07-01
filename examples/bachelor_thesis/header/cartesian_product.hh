#pragma once

/**
 * @file cartesian_product.hh
 * @brief Discrete Cartesian product of CW-complexes embedded in OpenVolumeMesh.
 *
 * Given two CW-complexes A and B with topological dimensions p and q, the
 * Cartesian product A x B is a CW-complex of dimension p+q whose k-cells
 * are formed by all pairs (sigma, tau) of cells with dim(sigma)+dim(tau) = k.
 * This module enumerates every such pair up to dimension three, which is
 * the maximum supported by OpenVolumeMesh's 3-manifold topology kernel:
 *
 *   k = 0 : V_A x V_B                          -> Vertices
 *   k = 1 : V_A x E_B  +  E_A x V_B            -> Edges
 *   k = 2 : V_A x F_B  +  F_A x V_B  +  E_A x E_B  -> Faces (Quad if E x E)
 *   k = 3 : V_A x C_B  +  C_A x V_B
 *         + E_A x F_B  +  F_A x E_B            -> Cells
 *
 */

#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <vector>
#include <map>
#include <algorithm>
#include <stdexcept>

namespace CartesianOps {

    /**
     * @brief Determines the topological dimension of an OpenVolumeMesh complex.
     *
     * @tparam MeshType  Any OpenVolumeMesh kernel type (Mesh4D, Mesh5D, Mesh6D).
     * @param  mesh      The complex whose dimension is to be evaluated.
     * @return The dimension of the highest-dimensional cell present
     *         (3 for cells, 2 for faces, 1 for edges, 0 for isolated vertices).
     */
    template <typename MeshType>
    inline int get_dim(const MeshType& mesh) {
        if (mesh.n_cells() > 0) return 3;
        if (mesh.n_faces() > 0) return 2;
        if (mesh.n_edges() > 0) return 1;
        return 0;
    }


    // INTERNAL CELL ASSEMBLY HELPER
    namespace detail {

        /**
         * @brief Retrieves or instantiates a FaceHandle for a given vertex set.
         *
         * Cell assembly produces overlapping face requests when adjacent
         * cells share a common boundary face. To make sure that each
         * geometric face is represented by exactly one FaceHandle in the
         * mesh kernel, a persistent dictionary keyed by an
         * order-independent signature, namely the sorted vector of vertex
         * indices is mantained. The first request creates the face; subsequent requests
         * for the same vertex set retrieve the cached handle.
         *
         * Triangle canonicalisation. For triangular faces (3 vertices), we
         * additionally enforce a globally deterministic winding by passing
         * the vertices to mesh.add_face in increasing index order. This
         * makes hf0 of every triangle equal the lexicographic cycle v_min ->
         * v_mid -> v_max -> v_min.
         *
         * Non-triangular faces (quads from E x E products, prism lateral
         * walls, hexahedral side faces) are left in the caller's order.
         *
         * @tparam MeshT    The OpenVolumeMesh kernel type.
         * @tparam VH       The VertexHandle type.
         * @tparam FH       The FaceHandle type.
         * @param  mesh     The mesh under construction.
         * @param  verts    Vertices defining the face in their connected order.
         * @param  face_map Dictionary of (sorted-index-key -> FH).
         * @return The cached or newly created handle for the requested face.
         */
        template <typename MeshT, typename VH, typename FH>
        inline FH get_or_create_face(MeshT& mesh, const std::vector<VH>& verts,
                                     std::map<std::vector<int>, FH>& face_map)
        {
            // 1. Build the order-independent signature.
            std::vector<int> key;
            for (auto v : verts) key.push_back(v.idx());
            std::sort(key.begin(), key.end());

            // 2. Return the existing handle if any face with
            // this exact vertex set has already been instantiated.
            auto it = face_map.find(key);
            if (it != face_map.end()) return it->second;

            // 3. Instantiate the new face.
            FH fh;
            if (verts.size() == 3) {
                // Triangle: Lexicographic order
                std::vector<VH> canonical = {VH(key[0]), VH(key[1]), VH(key[2])};
                fh = mesh.add_face(canonical);
            } else {
                // Non-triangle: preserve caller's winding.
                fh = mesh.add_face(verts);
            }
            face_map[key] = fh;
            return fh;
        }

        /**
         * @brief Selects one unclaimed half-face per boundary face so the cell
         *        can be assembled (First-Come-First-Served).
         *
         *
         * @tparam MeshT  The OpenVolumeMesh kernel type.
         * @tparam FH     The FaceHandle type.
         * @tparam HFH    The HalfFaceHandle type.
         * @param  mesh   The mesh under construction.
         * @param  faces  The faces bounding the intended cell.
         * @return Half-face handles ready to be passed to mesh.add_cell.
         * @throws std::runtime_error  If a face is already claimed on both sides.
         */
        template <typename MeshT, typename FH, typename HFH>
        inline std::vector<HFH> resolve_cell_halffaces(const MeshT& mesh, const std::vector<FH>& faces) {
            std::vector<HFH> result_halffaces;

            for(const auto& fh : faces) {
                // Both possible orientations of the face.
                auto hf0 = mesh.halfface_handle(fh, 0);
                auto hf1 = mesh.halfface_handle(fh, 1);

                // FCFS resolution.
                // If hf0 is unclaimed, take it; otherwise take hf1.
                if (!mesh.incident_cell(hf0).is_valid())
                    result_halffaces.push_back(hf0);
                else if (!mesh.incident_cell(hf1).is_valid())
                    result_halffaces.push_back(hf1);
                else
                    // Both half-faces already belong to other cells, which indicates a non-manifold configuration.
                    throw std::runtime_error("Topology Error: Face shared by >2 cells.");
            }

            return result_halffaces;
        }

    } // namespace detail



    // BINARY CARTESIAN PRODUCT (A x B)
    /**
     * @brief Computes the discrete Cartesian product of two CW-complexes.
     *
     * Enumerates every (p-cell of A) x (q-cell of B) pair with p+q in
     * {0, 1, 2, 3} and constructs the corresponding cell in the product
     * complex. Vertex pairs collapse via the bijection
     *
     *     pid(vA, vB) = vA.idx() * |V(B)| + vB.idx()
     *
     * The combiner functor combines the coordinates of the base shapes to the target shape.
     *
     * @tparam MeshOut    Output mesh type (e.g. Mesh4D, Mesh5D, Mesh6D).
     * @tparam MeshA      First input mesh type.
     * @tparam MeshB      Second input mesh type.
     * @tparam Combiner   Callable with signature (PointA, PointB) -> PointOut.
     * @param  A          First factor.
     * @param  B          Second factor.
     * @param  combine_coords  Functor fusing a vertex from A and one from B
     *                         into the ambient point of the product.
     * @return The CW-complex A x B as a MeshOut instance.
     * @throws std::runtime_error  If dim(A) + dim(B) > 3.
     */
    template <typename MeshOut, typename MeshA, typename MeshB, typename Combiner>
    inline MeshOut cartesian_product(const MeshA& A, const MeshB& B, Combiner combine_coords) {

        // OpenVolumeMesh's 3-manifold kernel cannot represent k-cells for k > 3
        if (get_dim(A) + get_dim(B) > 3) {
            throw std::runtime_error("Topology Error: Sum of topological dimensions exceeds 3.");
        }

        MeshOut C;
        using VH  = OpenVolumeMesh::VertexHandle;
        using FH  = OpenVolumeMesh::FaceHandle;
        using HFH = OpenVolumeMesh::HalfFaceHandle;

        // Persistent face dictionary
        std::map<std::vector<int>, FH> face_map;

        // Product-vertex bijection.
        // Maps a pair (vA, vB) to a single vertex handle in the product.
        auto pid = [&](OpenVolumeMesh::VertexHandle vA, OpenVolumeMesh::VertexHandle vB) {
            return VH(vA.idx() * B.n_vertices() + vB.idx());
        };


        // 1. V x V  ->  Vertices
        for (auto vA_it = A.vertices_begin(); vA_it != A.vertices_end(); ++vA_it) {
            for (auto vB_it = B.vertices_begin(); vB_it != B.vertices_end(); ++vB_it) {
                C.add_vertex(combine_coords(A.vertex(*vA_it), B.vertex(*vB_it)));
            }
        }


        // 2. V x E  and  E x V  ->  Edges
        // Replicate each edge of B at every vertex of A.
        for (auto vA_it = A.vertices_begin(); vA_it != A.vertices_end(); ++vA_it) {
            for (auto eB_it = B.edges_begin(); eB_it != B.edges_end(); ++eB_it) {
                C.add_edge(pid(*vA_it, B.edge(*eB_it).from_vertex()),
                           pid(*vA_it, B.edge(*eB_it).to_vertex()));
            }
        }
        // Replicate each edge of A at every vertex of B.
        for (auto eA_it = A.edges_begin(); eA_it != A.edges_end(); ++eA_it) {
            for (auto vB_it = B.vertices_begin(); vB_it != B.vertices_end(); ++vB_it) {
                C.add_edge(pid(A.edge(*eA_it).from_vertex(), *vB_it),
                           pid(A.edge(*eA_it).to_vertex(),   *vB_it));
            }
        }

        // ---------------------------------------------------------------
        // 3. E x E  ->  Faces (Quads)
        // ---------------------------------------------------------------
        // Sweeping a 1-cell along another 1-cell creates a 2-cell with four
        // vertices.
        // Pattern:
        //          (Top-Left)                          (Top-Right)
        //              b0 ------------------------------- b1
        //              |                                  |
        //              |                                  |
        //              |                                  |
        //              a0 ------------------------------- a1
        //          (Bottom-Left)                       (Bottom-Right)

        for (auto eA_it = A.edges_begin(); eA_it != A.edges_end(); ++eA_it) {
            for (auto eB_it = B.edges_begin(); eB_it != B.edges_end(); ++eB_it) {
                VH a0 = pid(A.edge(*eA_it).from_vertex(), B.edge(*eB_it).from_vertex());
                VH a1 = pid(A.edge(*eA_it).to_vertex(),   B.edge(*eB_it).from_vertex());
                VH b1 = pid(A.edge(*eA_it).to_vertex(),   B.edge(*eB_it).to_vertex());
                VH b0 = pid(A.edge(*eA_it).from_vertex(), B.edge(*eB_it).to_vertex());
                detail::get_or_create_face<MeshOut, VH, FH>(C, {a0, a1, b1, b0}, face_map);
            }
        }


        // 4. F x V  and  V x F  ->  Faces
        for (auto fA_it = A.faces_begin(); fA_it != A.faces_end(); ++fA_it) {
            for (auto vB_it = B.vertices_begin(); vB_it != B.vertices_end(); ++vB_it) {
                std::vector<VH> fv;
                for(auto v : A.face_vertices(*fA_it)) fv.push_back(pid(v, *vB_it));
                detail::get_or_create_face<MeshOut, VH, FH>(C, fv, face_map);
            }
        }
        for (auto vA_it = A.vertices_begin(); vA_it != A.vertices_end(); ++vA_it) {
            for (auto fB_it = B.faces_begin(); fB_it != B.faces_end(); ++fB_it) {
                std::vector<VH> fv;
                for(auto v : B.face_vertices(*fB_it)) fv.push_back(pid(*vA_it, v));
                detail::get_or_create_face<MeshOut, VH, FH>(C, fv, face_map);
            }
        }

        // 5a. F x E  ->  3-Cells (Triangular Prism)

        for (auto fA_it = A.faces_begin(); fA_it != A.faces_end(); ++fA_it) {
            for (auto eB_it = B.edges_begin(); eB_it != B.edges_end(); ++eB_it) {
                std::vector<VH> fA_v;
                for(auto fv : A.face_vertices(*fA_it)) fA_v.push_back(fv);

                auto b0 = B.edge(*eB_it).from_vertex();   // bottom layer
                auto b1 = B.edge(*eB_it).to_vertex();     // top layer

                std::vector<FH> bounds;

                // a non-triangulated face here signals an error
                if (fA_v.size() != 3) {
                    throw std::runtime_error(
                        "F x E product: expected a triangular face (triangulate "
                        "the 2-factor before forming 3-cells).");
                }

                // Triangular prism: a triangular cap on each B-layer (b0, b1)
                // plus three lateral quads, one per edge of A.

                // Bottom cap at b0, wound as in A:        v0 -> v1 -> v2.
                bounds.push_back(detail::get_or_create_face<MeshOut, VH, FH>(C,
                    {pid(fA_v[0], b0), pid(fA_v[1], b0), pid(fA_v[2], b0)}, face_map));

                // Top cap at b1, wound opposite to A:      v0 -> v2 -> v1
                bounds.push_back(detail::get_or_create_face<MeshOut, VH, FH>(C,
                    {pid(fA_v[0], b1), pid(fA_v[2], b1), pid(fA_v[1], b1)}, face_map));

                // Lateral quads, each wound bottom-edge -> up -> top-edge
                // (reversed) -> down:  (vi,b0) -> (vj,b0) -> (vj,b1) -> (vi,b1).
                bounds.push_back(detail::get_or_create_face<MeshOut, VH, FH>(C,
                    {pid(fA_v[0], b0), pid(fA_v[1], b0), pid(fA_v[1], b1), pid(fA_v[0], b1)}, face_map));  // edge v0->v1
                bounds.push_back(detail::get_or_create_face<MeshOut, VH, FH>(C,
                    {pid(fA_v[1], b0), pid(fA_v[2], b0), pid(fA_v[2], b1), pid(fA_v[1], b1)}, face_map));  // edge v1->v2
                bounds.push_back(detail::get_or_create_face<MeshOut, VH, FH>(C,
                    {pid(fA_v[2], b0), pid(fA_v[0], b0), pid(fA_v[0], b1), pid(fA_v[2], b1)}, face_map));  // edge v2->v0
                C.add_cell(detail::resolve_cell_halffaces<MeshOut, FH, HFH>(C, bounds));
            }
        }

        // 5b. E x F  ->  3-Cells

        for (auto eA_it = A.edges_begin(); eA_it != A.edges_end(); ++eA_it) {
            for (auto fB_it = B.faces_begin(); fB_it != B.faces_end(); ++fB_it) {
                std::vector<VH> fB_v;
                for(auto fv : B.face_vertices(*fB_it)) fB_v.push_back(fv);

                auto a0 = A.edge(*eA_it).from_vertex();   // bottom layer
                auto a1 = A.edge(*eA_it).to_vertex();     // top layer

                std::vector<FH> bounds;

                // Mirror of 5a with the roles of the two factors swapped
                if (fB_v.size() != 3) {
                    throw std::runtime_error(
                        "E x F product: expected a triangular face (triangulate "
                        "the 2-factor before forming 3-cells).");
                }

                // Bottom cap at a0, wound as in B:         v0 -> v1 -> v2.
                bounds.push_back(detail::get_or_create_face<MeshOut, VH, FH>(C,
                    {pid(a0, fB_v[0]), pid(a0, fB_v[1]), pid(a0, fB_v[2])}, face_map));
                // Top cap at a1, wound opposite to B:       v0 -> v2 -> v1.
                bounds.push_back(detail::get_or_create_face<MeshOut, VH, FH>(C,
                    {pid(a1, fB_v[0]), pid(a1, fB_v[2]), pid(a1, fB_v[1])}, face_map));

                // Lateral quads:  (a0,vi) -> (a0,vj) -> (a1,vj) -> (a1,vi).
                bounds.push_back(detail::get_or_create_face<MeshOut, VH, FH>(C,
                    {pid(a0, fB_v[0]), pid(a0, fB_v[1]), pid(a1, fB_v[1]), pid(a1, fB_v[0])}, face_map));  // edge v0->v1
                bounds.push_back(detail::get_or_create_face<MeshOut, VH, FH>(C,
                    {pid(a0, fB_v[1]), pid(a0, fB_v[2]), pid(a1, fB_v[2]), pid(a1, fB_v[1])}, face_map));  // edge v1->v2
                bounds.push_back(detail::get_or_create_face<MeshOut, VH, FH>(C,
                    {pid(a0, fB_v[2]), pid(a0, fB_v[0]), pid(a1, fB_v[0]), pid(a1, fB_v[2])}, face_map));  // edge v2->v0
                C.add_cell(detail::resolve_cell_halffaces<MeshOut, FH, HFH>(C, bounds));
            }
        }


        // 6. C x V  and  V x C  ->  3-Cells

        for (auto cA_it = A.cells_begin(); cA_it != A.cells_end(); ++cA_it) {
            for (auto vB_it = B.vertices_begin(); vB_it != B.vertices_end(); ++vB_it) {

                std::vector<FH> bounds;

                // Walk the cell's halffaces
                for (auto hf_it = A.chf_iter(*cA_it); hf_it.valid(); ++hf_it) {
                    std::vector<VH> fv;
                    // Translate every vertex of the halfface into the product index space.
                    for(auto v : A.face_vertices(A.face_handle(*hf_it)))
                        fv.push_back(pid(v, *vB_it));
                    bounds.push_back(detail::get_or_create_face<MeshOut, VH, FH>(C, fv, face_map));
                }
                C.add_cell(detail::resolve_cell_halffaces<MeshOut, FH, HFH>(C, bounds));
            }
        }
        for (auto vA_it = A.vertices_begin(); vA_it != A.vertices_end(); ++vA_it) {
            for (auto cB_it = B.cells_begin(); cB_it != B.cells_end(); ++cB_it) {

                std::vector<FH> bounds;
                for (auto hf_it = B.chf_iter(*cB_it); hf_it.valid(); ++hf_it) {
                    std::vector<VH> fv;
                    for(auto v : B.face_vertices(B.face_handle(*hf_it)))
                        fv.push_back(pid(*vA_it, v));
                    bounds.push_back(detail::get_or_create_face<MeshOut, VH, FH>(C, fv, face_map));
                }
                C.add_cell(detail::resolve_cell_halffaces<MeshOut, FH, HFH>(C, bounds));
            }
        }

        return C;
    }


} // namespace CartesianOps
