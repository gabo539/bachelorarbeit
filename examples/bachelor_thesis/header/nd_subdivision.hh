#pragma once

/**
 * @file nd_subdivision.hh
 * @brief Tetrahedralization and refinement of mixed polyhedral meshes.
 *
 *   triangulate_prisms_to_tets : split every triangular prism into 3 tets.
 *   triangulate_quad_faces     : split every quad face of a surface into 2 triangles.
 *   barycentric_subdivide      : barycentric subdivision (1 tet -> 24).
 *   uniform_subdivide          : uniform 1 -> 8 subdivision.
 *   normalize_blocks           : re-project each factor coordinate block onto
 *                                its unit sphere / circle.
 */

#include "base_shape_generator.hh"
#include "cartesian_product.hh"
#include <OpenVolumeMesh/Core/Handles.hh>
#include <map>
#include <set>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <cassert>

namespace Subdivision {

// 1. SHARED TET ASSEMBLY HELPER

/**
 * @brief Adds a single tetrahedron to the mesh from four vertex handles.
 *
 * The four boundary triangles are assembled in canonical order, deduplicated
 * via the shared face_map, and oriented through the FCFS resolver.
 *
 * @tparam MeshT  The OpenVolumeMesh kernel type.
 * @tparam VH     The VertexHandle type.
 * @tparam FH     The FaceHandle type.
 * @tparam HFH    The HalfFaceHandle type.
 * @param  mesh      The mesh under construction.
 * @param  v0..v3    The four vertex handles defining the tetrahedron.
 *                   They must be pairwise distinct (precondition).
 * @param  face_map  Dictionary used by get_or_create_face for
 *                   sub-face deduplication across all tets in this pass.
 */
template <typename MeshT, typename VH, typename FH, typename HFH>
inline void add_tet_to_mesh(MeshT& mesh, VH v0, VH v1, VH v2, VH v3,
                            std::map<std::vector<int>, FH>& face_map) {

    assert(v0 != v1 && v0 != v2 && v0 != v3
        && v1 != v2 && v1 != v3 && v2 != v3
        && "add_tet_to_mesh: duplicate vertex handles");

    std::vector<FH> faces;

    // Four-triangle boundary of the tetrahedron. Faces shared with
    // neighbouring tets resolve to a single FaceHandle via the face_map.
    faces.push_back(CartesianOps::detail::get_or_create_face<MeshT, VH, FH>(mesh, {v0, v1, v2}, face_map));
    faces.push_back(CartesianOps::detail::get_or_create_face<MeshT, VH, FH>(mesh, {v0, v2, v3}, face_map));
    faces.push_back(CartesianOps::detail::get_or_create_face<MeshT, VH, FH>(mesh, {v0, v3, v1}, face_map));
    faces.push_back(CartesianOps::detail::get_or_create_face<MeshT, VH, FH>(mesh, {v1, v2, v3}, face_map));

    auto hfs = CartesianOps::detail::resolve_cell_halffaces<MeshT, FH, HFH>(mesh, faces);
    assert(hfs.size() == 4 && "add_tet_to_mesh: expected 4 half-faces for a tetrahedron");
    mesh.add_cell(hfs);
}



// 2. STAGE A: TRIANGULATE PRISMS INTO TETRAHEDRA
namespace prism_detail {

    /**
     * @brief Canonical layout of a triangular prism cell.
     *
     * Stores the bottom triangle (b[0], b[1], b[2]) and the top triangle
     * (t[0], t[1], t[2]) such that b[i] and t[i] are connected by a
     * vertical edge of the prism.
     */
    struct Layout {
        OpenVolumeMesh::VertexHandle b[3];   // bottom triangle vertices
        OpenVolumeMesh::VertexHandle t[3];   // top triangle, paired by index
    };

    /**
     * @brief Recovers the (bottom, top, vertical pairing) layout of a prism.
     *
     * Identifies the two triangular faces of the cell and labels one
     * "bottom" (arbitrary). For each bottom vertex, the corresponding top
     * vertex is found by walking the cell's edges and selecting the
     * vertical edge (the one that is neither in the bottom nor the top
     * triangle).
     *
     * @tparam MeshT  The OpenVolumeMesh kernel type.
     * @param  mesh   The mesh containing the prism.
     * @param  ch     Handle to the prism cell.
     * @return Layout populated with the 6 vertices and their pairing.
     * @throws std::runtime_error if the cell is not a triangular prism.
     */
    template <typename MeshT>
    inline Layout extract_layout(const MeshT& mesh, OpenVolumeMesh::CellHandle ch) {
        using VH = OpenVolumeMesh::VertexHandle;
        using FH = OpenVolumeMesh::FaceHandle;

        // 1. Find the two triangular faces among the cell's bounding faces.
        std::vector<FH> tris;
        for (auto cf_it = mesh.cf_iter(ch); cf_it.valid(); ++cf_it) {
            int nv = 0;
            for (auto fv_it = mesh.fv_iter(*cf_it); fv_it.valid(); ++fv_it) ++nv;
            if (nv == 3) tris.push_back(*cf_it);
        }
        if (tris.size() != 2) {
            throw std::runtime_error(
                "extract_layout: expected 2 triangular faces on a prism");
        }

        Layout layout;

        // 2. Label the first triangle "bottom" and read off its vertices.
        int i = 0;
        for (auto fv_it = mesh.fv_iter(tris[0]); fv_it.valid(); ++fv_it) {
            layout.b[i++] = *fv_it;
        }
        assert(i == 3);

        // 3. Collect the top triangle's vertex set for membership testing.
        std::set<VH> top_set;
        for (auto fv_it = mesh.fv_iter(tris[1]); fv_it.valid(); ++fv_it) {
            top_set.insert(*fv_it);
        }

        // 4. For each bottom vertex, find its vertical counterpart by
        //    locating the cell edge that goes from b[j] to a top vertex.
        for (int j = 0; j < 3; ++j) {
            VH b_j = layout.b[j];
            bool found = false;
            for (auto ce_it = mesh.ce_iter(ch); ce_it.valid(); ++ce_it) {
                auto e = mesh.edge(*ce_it);
                VH from = e.from_vertex();
                VH to   = e.to_vertex();
                VH other;
                if (from == b_j)      other = to;
                else if (to == b_j)   other = from;
                else continue;
                if (top_set.count(other)) {
                    layout.t[j] = other;
                    found = true;
                    break;
                }
            }
            if (!found) {
                throw std::runtime_error(
                    "extract_layout: no vertical pair found for a bottom vertex");
            }
        }

        return layout;
    }

    /**
     * @brief Returns 0 or 1 depending on which diagonal a quad (p0,p1,p2,p3)
     *        splits along under the lowest-index rule.
     */
    inline int pick_diagonal(OpenVolumeMesh::VertexHandle p0,
                             OpenVolumeMesh::VertexHandle p1,
                             OpenVolumeMesh::VertexHandle p2,
                             OpenVolumeMesh::VertexHandle p3)
    {
        const int idx[4] = {p0.idx(), p1.idx(), p2.idx(), p3.idx()};
        int min_pos = 0;
        for (int i = 1; i < 4; ++i)
            if (idx[i] < idx[min_pos]) min_pos = i;
        return min_pos % 2;
    }

    /**
     * @brief Splits a triangular prism into 3 tetrahedra and emits them.
     *
     * @tparam MeshT     The OpenVolumeMesh kernel type.
     * @param  new_mesh  Destination mesh (the three tets are appended).
     * @param  L         The prism's layout.
     * @param  v_remap   Old-to-new vertex handle translation.
     * @param  face_map  Persistent face dictionary for the new mesh.
     */
    template <typename MeshT>
    inline void split_to_tets(MeshT& new_mesh,
                              const Layout& L,
                              const std::vector<OpenVolumeMesh::VertexHandle>& v_remap,
                              std::map<std::vector<int>, OpenVolumeMesh::FaceHandle>& face_map)
    {
        using VH  = OpenVolumeMesh::VertexHandle;
        using FH  = OpenVolumeMesh::FaceHandle;
        using HFH = OpenVolumeMesh::HalfFaceHandle;

        const VH b0 = L.b[0], b1 = L.b[1], b2 = L.b[2];
        const VH t0 = L.t[0], t1 = L.t[1], t2 = L.t[2];

        // Step 1: Pick one diagonal per lateral quad via lowest-index rule.
        // Encoding:
        //   d01 == 0 -> diagonal (b0, t1);  d01 == 1 -> diagonal (b1, t0)
        //   d12 == 0 -> diagonal (b1, t2);  d12 == 1 -> diagonal (b2, t1)
        //   d20 == 0 -> diagonal (b2, t0);  d20 == 1 -> diagonal (b0, t2)
        const int d01 = pick_diagonal(b0, b1, t1, t0);
        const int d12 = pick_diagonal(b1, b2, t2, t1);
        const int d20 = pick_diagonal(b2, b0, t0, t2);

        // Step 2: Count diagonal-incidences per vertex.
        // Each vertex sits on two lateral quads, so its count is the sum of
        // two indicators. Counts always sum to 6 (three diagonals, two
        // endpoints each); the distribution is {0, 0, 1, 1, 2, 2}.
        const int c_b0 = (1 - d01) + d20;
        const int c_b1 = d01       + (1 - d12);
        const int c_b2 = (1 - d20) + d12;
        const int c_t0 = d01       + (1 - d20);
        const int c_t1 = (1 - d01) + d12;
        const int c_t2 = (1 - d12) + d20;

        // Step 3: Categorise into spine / mid / apex on the bottom triangle
        // and on the top triangle. The categorisation is unique because
        // exactly one bottom vertex has count 2 (spine), one has count 1
        // (mid), one has count 0 (apex); same on the top.
        VH spine_b, mid_b, apex_b;
        {
            const VH bv[3] = {b0, b1, b2};
            const int cc[3] = {c_b0, c_b1, c_b2};
            for (int i = 0; i < 3; ++i) {
                if      (cc[i] == 2) spine_b = bv[i];
                else if (cc[i] == 1) mid_b   = bv[i];
                else                 apex_b  = bv[i];
            }
        }
        VH spine_t, mid_t, apex_t;
        {
            const VH tv[3] = {t0, t1, t2};
            const int cc[3] = {c_t0, c_t1, c_t2};
            for (int i = 0; i < 3; ++i) {
                if      (cc[i] == 2) spine_t = tv[i];
                else if (cc[i] == 1) mid_t   = tv[i];
                else                 apex_t  = tv[i];
            }
        }


        // Step 4: Translate to new-mesh handles and emit the three tets.
        auto R = [&](VH old) { return v_remap[old.idx()]; };
        const VH p   = R(spine_b);
        const VH q   = R(spine_t);
        const VH rb  = R(mid_b);
        const VH ab  = R(apex_b);
        const VH rt  = R(mid_t);
        const VH at  = R(apex_t);

        add_tet_to_mesh<MeshT, VH, FH, HFH>(new_mesh, p, q, rb, ab, face_map);
        add_tet_to_mesh<MeshT, VH, FH, HFH>(new_mesh, p, q, rb, rt, face_map);
        add_tet_to_mesh<MeshT, VH, FH, HFH>(new_mesh, p, q, rt, at, face_map);
    }

} // namespace prism_detail

/**
 * @brief Replaces every prism cell with three tetrahedra (in-place).
 *
 * Iterates the input mesh; tetrahedral cells are copied unchanged,
 * triangular-prism cells are decomposed into three tetrahedra
 * and any other cell type triggers a runtime error.
 *
 * @tparam MeshT  The OpenVolumeMesh kernel type.
 * @param  mesh   Input mesh of mixed cell types; replaced in-place by a
 *                pure-tetrahedral mesh.
 * @throws std::runtime_error  If any cell is a hexahedron or otherwise
 *                             unsupported.
 */
template <typename MeshT>
inline void triangulate_prisms_to_tets(MeshT& mesh) {
    using VH  = OpenVolumeMesh::VertexHandle;
    using FH  = OpenVolumeMesh::FaceHandle;
    using HFH = OpenVolumeMesh::HalfFaceHandle;

    std::cout << "Triangulating prism cells (" << mesh.n_cells()
              << " input cells) into tetrahedra..." << std::endl;

    MeshT new_mesh;

    // Carry over all original vertices; their old-handle -> new-handle
    std::vector<VH> v_old_to_new(mesh.n_vertices());
    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
        v_old_to_new[v_it->idx()] = new_mesh.add_vertex(mesh.vertex(*v_it));
    }

    std::map<std::vector<int>, FH> face_map;

    for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it) {
        // Count the cell's vertices to decide its type.
        int nv = 0;
        std::vector<VH> cv;
        for (auto cv_it = mesh.cv_iter(*c_it); cv_it.valid(); ++cv_it) {
            cv.push_back(*cv_it); ++nv;
        }

        if (nv == 4) {
            // Tetrahedron: copy as-is. Vertex order is preserved
            add_tet_to_mesh<MeshT, VH, FH, HFH>(new_mesh,
                v_old_to_new[cv[0].idx()],
                v_old_to_new[cv[1].idx()],
                v_old_to_new[cv[2].idx()],
                v_old_to_new[cv[3].idx()],
                face_map);
        }
        else if (nv == 6) {
            // Triangular prism: recover its layout and split into 3 tets.
            const auto layout = prism_detail::extract_layout(mesh, *c_it);
            prism_detail::split_to_tets(new_mesh, layout, v_old_to_new, face_map);
        }
        else {
            throw std::runtime_error(
                "triangulate_prisms_to_tets: unsupported cell with "
                + std::to_string(nv) + " vertices");
        }
    }

    mesh = new_mesh;

    std::cout << "  -> Result: " << mesh.n_cells() << " tetrahedra." << std::endl;
}



// 2b. SURFACE TRIANGULATION: SPLIT QUAD FACES INTO TRIANGLES
/**
 * @brief Splits every quadrilateral face of a cell-free surface into two
 *        triangles.
 *
 * @tparam MeshT  The OpenVolumeMesh kernel type.
 * @param  mesh   A cell-free surface; replaced in place by its triangulation.
 * @throws std::runtime_error if the mesh has cells, or any face has fewer
 *         than 3 or more than 4 vertices.
 */
template <typename MeshT>
inline void triangulate_quad_faces(MeshT& mesh) {
    using VH = OpenVolumeMesh::VertexHandle;

    if (mesh.n_cells() != 0) {
        throw std::runtime_error(
            "triangulate_quad_faces: expected a cell-free surface, but the "
            "mesh has " + std::to_string(mesh.n_cells()) + " cell(s).");
    }

    std::cout << "Triangulating quad faces (" << mesh.n_faces()
              << " input faces)..." << std::endl;

    MeshT new_mesh;

    // Carry over all vertices in their original handle order.
    std::vector<VH> v_old_to_new(mesh.n_vertices());
    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
        v_old_to_new[v_it->idx()] = new_mesh.add_vertex(mesh.vertex(*v_it));
    }

    long quads = 0, tris = 0;

    for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
        // The face's vertex cycle, translated into the new index space.
        std::vector<VH> v;
        for (auto fv : mesh.face_vertices(*f_it)) v.push_back(v_old_to_new[fv.idx()]);

        if (v.size() == 3) {
            new_mesh.add_face({v[0], v[1], v[2]});
            ++tris;
        } else if (v.size() == 4) {
            // Lowest-index diagonal: both triangles inherit the quad winding.
            //   diag (v0,v2): (v0,v1,v2) + (v0,v2,v3)
            //   diag (v1,v3): (v1,v2,v3) + (v1,v3,v0)
            const int d = prism_detail::pick_diagonal(v[0], v[1], v[2], v[3]);
            if (d == 0) {
                new_mesh.add_face({v[0], v[1], v[2]});
                new_mesh.add_face({v[0], v[2], v[3]});
            } else {
                new_mesh.add_face({v[1], v[2], v[3]});
                new_mesh.add_face({v[1], v[3], v[0]});
            }
            ++quads;
        } else {
            throw std::runtime_error(
                "triangulate_quad_faces: unsupported face with "
                + std::to_string(v.size()) + " vertices");
        }
    }

    mesh = new_mesh;

    std::cout << "  -> Split " << quads << " quad(s), copied " << tris
              << " triangle(s); result has " << mesh.n_faces()
              << " triangular faces." << std::endl;
}


// 3. STAGE B: BARYCENTRIC SUBDIVISION OF A PURE TET MESH
/**
 * @brief Replaces a pure tetrahedral mesh with its barycentric subdivision.
 *
 * Precondition: every cell of the input mesh is a tetrahedron. This is
 * checked at runtime; the function aborts with a runtime error if any
 * cell has a different vertex count.
 *
 * @tparam MeshT  The OpenVolumeMesh kernel type.
 * @param  mesh   Pure-tetrahedral input mesh; replaced in-place by its
 *                barycentric subdivision.
 * @throws std::runtime_error  If any cell is not a tetrahedron.
 */
template <typename MeshT>
inline void barycentric_subdivide(MeshT& mesh) {
    std::cout << "Barycentrically subdividing tet mesh ("
              << mesh.n_cells() << " tetrahedra)..." << std::endl;

    using VH     = OpenVolumeMesh::VertexHandle;
    using FH     = OpenVolumeMesh::FaceHandle;
    using HFH    = OpenVolumeMesh::HalfFaceHandle;
    using PointT = typename MeshT::PointT;

    // Precondition check: Mesh has to be tetrahedral
    for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it) {
        int nv = 0;
        for (auto cv_it = mesh.cv_iter(*c_it); cv_it.valid(); ++cv_it) ++nv;
        if (nv != 4) {
            throw std::runtime_error(
                "barycentric_subdivide: input mesh must be pure-tet "
                "(call triangulate_prisms_to_tets first)");
        }
    }

    MeshT new_mesh;


    // Stage 1a: Carry over the original 0-cells.
    std::vector<VH> old_v_to_new(mesh.n_vertices());
    for(auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
        old_v_to_new[v_it->idx()] = new_mesh.add_vertex(mesh.vertex(*v_it));
    }

    // Centroid helper. Initialising the accumulator from a real vertex and
    // subtracting it from itself yields a zero of PointT in a way that does
    // not rely on Eigen's setZero or OVM's vectorize, neither of which is
    // portable across our 4D / 5D / 6D mesh kernels.
    auto centroid = [&](auto first_it) {
        assert(first_it.valid() && "centroid: empty vertex circulator");
        PointT sum = mesh.vertex(*first_it);
        sum -= sum;
        int count = 0;
        for (auto it = first_it; it.valid(); ++it) {
            sum += mesh.vertex(*it);
            ++count;
        }
        return sum / static_cast<double>(count);
    };


    // Stage 1b: Edge midpoints, face centroids, cell centroids.
    std::vector<VH> edge_midpoints(mesh.n_edges());
    for (auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it) {
        auto he = mesh.halfedge_handle(*e_it, 0);
        PointT mid = (mesh.vertex(mesh.from_vertex_handle(he))
                    + mesh.vertex(mesh.to_vertex_handle(he))) * 0.5;
        edge_midpoints[e_it->idx()] = new_mesh.add_vertex(mid);
    }

    std::vector<VH> face_barycenters(mesh.n_faces());
    for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
        face_barycenters[f_it->idx()] = new_mesh.add_vertex(centroid(mesh.fv_iter(*f_it)));
    }

    std::vector<VH> cell_barycenters(mesh.n_cells());
    for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it) {
        cell_barycenters[c_it->idx()] = new_mesh.add_vertex(centroid(mesh.cv_iter(*c_it)));
    }


    // Stage 2: The flag loop. One tet per flag (v, e, f, c).
    std::map<std::vector<int>, FH> face_map;

    for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it) {
        VH v_Cell = cell_barycenters[c_it->idx()];

        for (auto hf_it = mesh.chf_iter(*c_it); hf_it.valid(); ++hf_it) {
            VH v_Face = face_barycenters[mesh.face_handle(*hf_it).idx()];


            for (auto he_it = mesh.hfhe_iter(*hf_it); he_it.valid(); ++he_it) {
                VH v_Edge  = edge_midpoints[mesh.edge_handle(*he_it).idx()];
                VH v_Start = old_v_to_new[mesh.from_vertex_handle(*he_it).idx()];
                VH v_End   = old_v_to_new[mesh.to_vertex_handle(*he_it).idx()];

                add_tet_to_mesh<MeshT, VH, FH, HFH>(new_mesh,
                    v_Start, v_Edge, v_Face, v_Cell, face_map);
                add_tet_to_mesh<MeshT, VH, FH, HFH>(new_mesh,
                    v_Edge,  v_End,  v_Face, v_Cell, face_map);
            }
        }
    }


    // Stage 3: Replace the input mesh with its subdivision.
    mesh = new_mesh;

    std::cout << "  -> Result: " << mesh.n_cells() << " tetrahedra." << std::endl;
}



// 3b. STAGE B (ALTERNATIVE): UNIFORM 1->8 REFINEMENT OF A PURE TET MESH
/**
 * @brief Replaces a pure tetrahedral mesh with its uniform 1->8 refinement.
 *
 * @tparam MeshT  The OpenVolumeMesh kernel type.
 * @param  mesh   Pure-tetrahedral input mesh; replaced in-place by its uniform
 *                1->8 refinement.
 * @throws std::runtime_error  If any cell is not a tetrahedron.
 */
template <typename MeshT>
inline void uniform_subdivide(MeshT& mesh) {
    using VH     = OpenVolumeMesh::VertexHandle;
    using FH     = OpenVolumeMesh::FaceHandle;
    using HFH    = OpenVolumeMesh::HalfFaceHandle;
    using PointT = typename MeshT::PointT;

    std::cout << "Uniformly (1->8) subdividing tet mesh ("
              << mesh.n_cells() << " tetrahedra)..." << std::endl;

    // Precondition: pure-tet input.
    for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it) {
        int nv = 0;
        for (auto cv_it = mesh.cv_iter(*c_it); cv_it.valid(); ++cv_it) ++nv;
        if (nv != 4) {
            throw std::runtime_error(
                "uniform_subdivide: input mesh must be pure-tet "
                "(call triangulate_prisms_to_tets first)");
        }
    }

    MeshT new_mesh;

    // Carry over the original vertices (handle order preserved).
    std::vector<VH> old_v_to_new(mesh.n_vertices());
    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
        old_v_to_new[v_it->idx()] = new_mesh.add_vertex(mesh.vertex(*v_it));
    }

    // One midpoint per edge, shared by every incident cell. Keyed by the
    // sorted original-vertex pair so a cell can look its midpoints up directly.
    std::map<std::pair<int, int>, VH> mid_of;
    for (auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it) {
        auto he = mesh.halfedge_handle(*e_it, 0);
        const VH from = mesh.from_vertex_handle(he);
        const VH to   = mesh.to_vertex_handle(he);
        const PointT mid = (mesh.vertex(from) + mesh.vertex(to)) * 0.5;
        int a = from.idx(), b = to.idx();
        if (a > b) std::swap(a, b);
        mid_of[{a, b}] = new_mesh.add_vertex(mid);
    }

    auto mid = [&](VH u, VH w) {
        int a = u.idx(), b = w.idx();
        if (a > b) std::swap(a, b);
        return mid_of.at({a, b});
    };

    std::map<std::vector<int>, FH> face_map;

    for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it) {
        // Original vertices, sorted by index for a deterministic split layout.
        std::vector<VH> cv;
        for (auto cv_it = mesh.cv_iter(*c_it); cv_it.valid(); ++cv_it) cv.push_back(*cv_it);
        std::sort(cv.begin(), cv.end(), [](VH a, VH b){ return a.idx() < b.idx(); });

        const VH v0 = old_v_to_new[cv[0].idx()];
        const VH v1 = old_v_to_new[cv[1].idx()];
        const VH v2 = old_v_to_new[cv[2].idx()];
        const VH v3 = old_v_to_new[cv[3].idx()];

        const VH m01 = mid(cv[0], cv[1]), m02 = mid(cv[0], cv[2]), m03 = mid(cv[0], cv[3]);
        const VH m12 = mid(cv[1], cv[2]), m13 = mid(cv[1], cv[3]), m23 = mid(cv[2], cv[3]);

        // 4 corner tets.
        add_tet_to_mesh<MeshT, VH, FH, HFH>(new_mesh, v0, m01, m02, m03, face_map);
        add_tet_to_mesh<MeshT, VH, FH, HFH>(new_mesh, v1, m01, m12, m13, face_map);
        add_tet_to_mesh<MeshT, VH, FH, HFH>(new_mesh, v2, m02, m12, m23, face_map);
        add_tet_to_mesh<MeshT, VH, FH, HFH>(new_mesh, v3, m03, m13, m23, face_map);

        // 4 inner tets: central octahedron cut along the (m01, m23) diagonal.
        add_tet_to_mesh<MeshT, VH, FH, HFH>(new_mesh, m01, m23, m02, m12, face_map);
        add_tet_to_mesh<MeshT, VH, FH, HFH>(new_mesh, m01, m23, m12, m13, face_map);
        add_tet_to_mesh<MeshT, VH, FH, HFH>(new_mesh, m01, m23, m13, m03, face_map);
        add_tet_to_mesh<MeshT, VH, FH, HFH>(new_mesh, m01, m23, m03, m02, face_map);
    }
    mesh = new_mesh;

    std::cout << "  -> Result: " << mesh.n_cells() << " tetrahedra." << std::endl;
}



// 4. NORMALISATION

/**
 * @brief Re-normalises consecutive coordinate blocks of every vertex to unit
 *        length, in place.
 * @tparam MeshT       The OpenVolumeMesh kernel type.
 * @param  mesh        Mesh whose vertices carry product coordinates.
 * @param  block_sizes Sizes of the consecutive factor blocks; must sum to the
 *                     ambient dimension of the mesh's point type.
 * @throws std::runtime_error if the block sizes do not sum to the point
 *         dimension.
 */
template <typename MeshT>
inline void normalize_blocks(MeshT& mesh, const std::vector<int>& block_sizes) {
    using PointT = typename MeshT::PointT;

    int total = 0;
    for (int s : block_sizes) total += s;
    if (total != static_cast<int>(PointT::size())) {
        throw std::runtime_error(
            "normalize_blocks: block sizes sum to " + std::to_string(total)
            + " but the point dimension is "
            + std::to_string(static_cast<int>(PointT::size())));
    }

    long touched = 0;
    for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
        PointT p = mesh.vertex(*v_it);
        int offset = 0;
        for (int sz : block_sizes) {
            double n2 = 0.0;
            for (int k = 0; k < sz; ++k) n2 += p[offset + k] * p[offset + k];
            const double n = std::sqrt(n2);
            if (n > 1e-12) {
                for (int k = 0; k < sz; ++k) p[offset + k] /= n;
            }
            offset += sz;
        }
        mesh.set_vertex(*v_it, p);
        ++touched;
    }

    std::cout << "Normalising factor blocks {";
    for (size_t i = 0; i < block_sizes.size(); ++i)
        std::cout << block_sizes[i] << (i + 1 < block_sizes.size() ? "," : "");
    std::cout << "} on " << touched << " vertices." << std::endl;
}

} // namespace Subdivision
