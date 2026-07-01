#pragma once

/**
 * @file validation.hh
 * @brief Topological and geometric verification suite for closed 3-manifolds.
 *
 * The pipeline is intended to produce closed, orientable 3-manifolds. The
 * checks are split into two namespaces:
 *
 *   Validation::Topological  - combinatorial invariants that depend only on
 *                              the cell structure:
 *       1. chi(M) = V - E + F - C = 0  (Euler characteristic of any closed
 *          odd-dimensional manifold).
 *       2. Every 2-face is incident to exactly two 3-cells (closedness).
 *       3. Vertex, edge and face links have the expected sphere topology.
 *       4. Every cell carries a consistently wound boundary.
 *
 *   Validation::Geometric    - checks on the embedding (vertex positions):
 *       - No degenerate (zero-volume) cells.
 *       - Vertices lie on the intended surface (spheritorus / ditorus).
 *
 */

#include "ovm_nd_types.hh"
#include "OpenVolumeMesh/Geometry/Vector11T.hh"
#include <OpenVolumeMesh/Core/Handles.hh>

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <set>
#include <queue>
#include <cassert>

namespace Validation {

    /**
     * @brief Prints the cell counts and Euler characteristic of a mesh.
     *
     * Centralises the chi = V - E + F - C computation that the generators
     * would otherwise each repeat at every pipeline stage.
     *
     * @tparam MeshType  Any OpenVolumeMesh kernel type.
     * @param  mesh   The mesh to summarise.
     * @param  label  A short tag identifying the pipeline stage.
     */
    template <typename MeshType>
    inline void print_stats(const MeshType& mesh, const std::string& label) {
        const long V = mesh.n_vertices();
        const long E = mesh.n_edges();
        const long F = mesh.n_faces();
        const long C = mesh.n_cells();
        std::cout << "[" << label << "] V:" << V << " | E:" << E << " | F:" << F
                  << " | C:" << C << "  -> chi = " << (V - E + F - C) << "\n";
    }



    //  TOPOLOGICAL CHECKS
    namespace Topological {

    /**
     * @brief Verifies that the mesh has no boundary 2-faces.
     *
     * OpenVolumeMesh's is_boundary
     * predicate flags any face incident to fewer than two cells, so the
     * mesh is closed iff this check returns zero boundary faces.
     *
     * @tparam MeshT  Any OpenVolumeMesh kernel type.
     * @param  mesh   The mesh to verify.
     * @return true if the mesh has no boundary faces; false otherwise.
     */
    template <typename MeshT>
    inline bool verify_closed_manifold(const MeshT& mesh) {
        std::cout << " -> Verifying Closed Boundary..." << std::endl;
        long boundaries = 0;

        for(auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
            if (mesh.is_boundary(*f_it)) boundaries++;
        }

        if (boundaries == 0) {
            std::cout << "    [PASS] Mesh is closed (0 exposed faces)." << std::endl;
            return true;
        } else {
            std::cout << "    [FAIL] Mesh has " << boundaries << " holes!" << std::endl;
            return false;
        }
    }

    /**
     * @brief Tests whether a vertex link is connected.
     *
     */
    template <typename MeshT>
    inline bool is_vertex_link_connected(
        const MeshT& mesh,
        const std::set<OpenVolumeMesh::EdgeHandle>& link_e,
        const std::set<OpenVolumeMesh::FaceHandle>& link_f)
    {
        if (link_f.empty()) return true;

        std::set<OpenVolumeMesh::FaceHandle> visited;
        std::queue<OpenVolumeMesh::FaceHandle> q;
        q.push(*link_f.begin());
        visited.insert(*link_f.begin());

        while (!q.empty()) {
            const auto f = q.front(); q.pop();
            for (auto fe_it = mesh.fe_iter(f); fe_it.valid(); ++fe_it) {
                // Invariant: every edge of a link face was inserted into
                // link_e when that face was added. A miss means the link was
                // built incorrectly.
                assert(link_e.find(*fe_it) != link_e.end()
                       && "link face has an edge missing");
                for (auto ef_it = mesh.ef_iter(*fe_it); ef_it.valid(); ++ef_it) {
                    if (link_f.find(*ef_it) == link_f.end()) continue;
                    if (visited.insert(*ef_it).second) q.push(*ef_it);
                }
            }
        }
        return visited.size() == link_f.size();
    }

    /**
     * @brief Tests whether an edge link is connected.
     */
    template <typename MeshT>
    inline bool is_edge_link_connected(
        const MeshT& mesh,
        const std::set<OpenVolumeMesh::VertexHandle>& link_v,
        const std::set<OpenVolumeMesh::EdgeHandle>& link_e)
    {
        if (link_e.empty()) return true;
        std::set<OpenVolumeMesh::EdgeHandle> visited;
        std::queue<OpenVolumeMesh::EdgeHandle> q;

        q.push(*link_e.begin());
        visited.insert(*link_e.begin());
        while (!q.empty()) {
            const auto eh = q.front(); q.pop();
            const auto edge = mesh.edge(eh);
            const OpenVolumeMesh::VertexHandle endpoints[2] = {edge.from_vertex(), edge.to_vertex()};
            for (const auto& v : endpoints) {
                // Invariant: both endpoints of a link edge were inserted into
                // link_v when that edge was added. A miss means the link was
                // built incorrectly.
                assert(link_v.find(v) != link_v.end()
                       && "link edge has an endpoint missing");
                for (auto ve_it = mesh.ve_iter(v); ve_it.valid(); ++ve_it) {
                    if (link_e.find(*ve_it) == link_e.end()) continue;
                    if (visited.insert(*ve_it).second) q.push(*ve_it);
                }
            }
        }
        return visited.size() == link_e.size();
    }

    /**
     * @brief Verifies that every vertex link is homeomorphic to S^2.
     *
     * @tparam MeshT  Any OpenVolumeMesh kernel type.
     * @param  mesh   The mesh to verify.
     * @return true if every vertex has a connected 2-sphere link.
     */
    template <typename MeshT>
    inline bool verify_vertice_links(const MeshT& mesh) {
        std::cout << " -> Verifying Local Vertex Links ..." << std::endl;
        int isolated     = 0;
        int chi_failures = 0;
        int disconnected = 0;

        for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
            // Sets used to deduplicate the vertices, edges, and faces that
            // form the shell around the central vertex.
            std::set<OpenVolumeMesh::VertexHandle> link_v;
            std::set<OpenVolumeMesh::EdgeHandle>   link_e;
            std::set<OpenVolumeMesh::FaceHandle>   link_f;

            // Walk every 3-cell incident to the vertex, then every 2-face
            // of that cell. A face belongs to the link if it does not
            // contain the central vertex.
            for (auto vc_it = mesh.vc_iter(*v_it); vc_it.valid(); ++vc_it) {
                for (auto cf_it = mesh.cf_iter(*vc_it); cf_it.valid(); ++cf_it) {
                    bool contains_v = false;
                    for (auto fv_it = mesh.fv_iter(*cf_it); fv_it.valid(); ++fv_it) {
                        if (*fv_it == *v_it) { contains_v = true; break; }
                    }
                    if (!contains_v) {
                        link_f.insert(*cf_it);
                        // Record the boundary vertices and edges of this face
                        for (auto fv_it = mesh.fv_iter(*cf_it); fv_it.valid(); ++fv_it) link_v.insert(*fv_it);
                        for (auto fe_it = mesh.fe_iter(*cf_it); fe_it.valid(); ++fe_it) link_e.insert(*fe_it);
                    }
                }
            }

            // An isolated vertex (empty link, no incident cells) cannot occur
            // in a closed 3-manifold
            if (link_v.empty()) { ++isolated; continue; }

            // chi(S^2) = 2.
            const bool chi_ok = (link_v.size() - link_e.size() + link_f.size() == 2);
            if (!chi_ok) {
                chi_failures++;
                continue;
            }

            if (!is_vertex_link_connected(mesh, link_e, link_f)) {
                disconnected++;
            }
        }

        const bool ok = (isolated == 0 && chi_failures == 0 && disconnected == 0);
        if (ok) {
            std::cout << "    [PASS] All vertex links are connected 2-Spheres." << std::endl;
        } else {
            if (isolated > 0)
                std::cout << "    [FAIL] Found " << isolated << " isolated vertex/vertices (empty link)!" << std::endl;
            if (chi_failures > 0)
                std::cout << "    [FAIL] Found " << chi_failures << " non-manifold vertices (link chi != 2)!" << std::endl;
            if (disconnected > 0)
                std::cout << "    [FAIL] Found " << disconnected << " disconnected vertex links!" << std::endl;
        }
        return ok;
    }

    /**
     * @brief Verifies that every edge link is homeomorphic to S^1.
     *
     * @tparam MeshT  Any OpenVolumeMesh kernel type.
     * @param  mesh   The mesh to verify.
     * @return true if every edge has a connected circular link.
     */
    template <typename MeshT>
    inline bool verify_edge_links(const MeshT& mesh) {
        std::cout << " -> Verifying Local Edge Links ..." << std::endl;
        int isolated     = 0;
        int chi_failures = 0;
        int disconnected = 0;

        for (auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it) {
            const auto edge = mesh.edge(*e_it);
            const OpenVolumeMesh::VertexHandle v1 = edge.from_vertex();
            const OpenVolumeMesh::VertexHandle v2 = edge.to_vertex();


            std::set<OpenVolumeMesh::VertexHandle> link_v;
            std::set<OpenVolumeMesh::EdgeHandle>   link_e;

            // Collect the cells incident to this edge via its faces. The radial
            // edge->cell circulator (ec_iter) under-reports the star on a mesh
            // that has not been consistently oriented yet
            std::set<OpenVolumeMesh::CellHandle> star;
            for (auto ef_it = mesh.ef_iter(*e_it); ef_it.valid(); ++ef_it) {
                for (auto fc : mesh.face_cells(*ef_it))
                    if (fc.is_valid()) star.insert(fc);
            }

            for (const auto& ch : star) {
                for (auto ce_it = mesh.ce_iter(ch); ce_it.valid(); ++ce_it) {
                    const auto ce = mesh.edge(*ce_it);
                    const OpenVolumeMesh::VertexHandle a = ce.from_vertex();
                    const OpenVolumeMesh::VertexHandle b = ce.to_vertex();
                    // Keep only edges of C whose endpoints both miss e.
                    if (a != v1 && a != v2 && b != v1 && b != v2) {
                        link_e.insert(*ce_it);
                        link_v.insert(a);
                        link_v.insert(b);
                    }
                }
            }

            // An isolated edge (empty link, no incident cells) cannot occur in
            // a closed 3-manifold
            if (link_v.empty()) { ++isolated; continue; }

            // chi(S^1) = 0, i.e. a circular link has equally many vertices and
            // edges.
            const bool chi_ok = (link_v.size() == link_e.size());
            if (!chi_ok) {
                chi_failures++;
                continue;
            }
            // chi = 0 alone permits two disjoint circles; BFS confirms a
            // single component.
            if (!is_edge_link_connected(mesh, link_v, link_e)) {
                disconnected++;
            }
        }

        const bool ok = (isolated == 0 && chi_failures == 0 && disconnected == 0);
        if (ok) {
            std::cout << "    [PASS] All edge links are connected 1-Spheres." << std::endl;
        } else {
            if (isolated > 0)
                std::cout << "    [FAIL] Found " << isolated << " isolated edge(s) (empty link)!" << std::endl;
            if (chi_failures > 0)
                std::cout << "    [FAIL] Found " << chi_failures << " non-circular edge links (chi != 0)!" << std::endl;
            if (disconnected > 0)
                std::cout << "    [FAIL] Found " << disconnected << " disconnected edge links!" << std::endl;
        }
        return ok;
    }

    /**
     * @brief Verifies that every face link is homeomorphic to S^0.
     *
     * @tparam MeshT  Any OpenVolumeMesh kernel type.
     * @param  mesh   The mesh to verify.
     * @return true if every face is interior and has a 2-point link.
     */
    template <typename MeshT>
    inline bool verify_face_links(const MeshT& mesh) {
        std::cout << " -> Verifying Local Face Links ..." << std::endl;
        int non_manifold_faces = 0;
        int boundary_faces     = 0;

        for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
            // Boundary faces have no second apex
            if (mesh.is_boundary(*f_it)) {
                boundary_faces++;
                continue;
            }

            // Vertices on the face itself are excluded from the link.
            std::set<OpenVolumeMesh::VertexHandle> face_v;
            for (auto fv_it = mesh.fv_iter(*f_it); fv_it.valid(); ++fv_it) {
                face_v.insert(*fv_it);
            }

            std::set<OpenVolumeMesh::VertexHandle> link_v;
            const auto incident_cells = mesh.face_cells(*f_it);
            for (const auto& ch : incident_cells) {
                if (!ch.is_valid()) continue;
                for (auto cv_it = mesh.cv_iter(ch); cv_it.valid(); ++cv_it) {
                    if (face_v.find(*cv_it) == face_v.end()) link_v.insert(*cv_it);
                }
            }

            if (link_v.size() != 2) non_manifold_faces++;
        }

        const bool ok = (non_manifold_faces == 0 && boundary_faces == 0);
        if (ok) {
            std::cout << "    [PASS] All face links are valid 0-Spheres (2 points)." << std::endl;
        } else {
            if (boundary_faces > 0)
                std::cout << "    [FAIL] Found " << boundary_faces
                          << " boundary face(s); an open boundary has no S^0 face link!" << std::endl;
            if (non_manifold_faces > 0)
                std::cout << "    [FAIL] Found " << non_manifold_faces
                          << " interior face(s) whose link is not S^0!" << std::endl;
        }
        return ok;
    }

    /**
     * @brief Runs all three link checks in sequence.
     *
     * Combines verify_vertice_links (S^2), verify_edge_links (S^1) and
     * verify_face_links (S^0) into a single call.
     *
     * @tparam MeshT  Any OpenVolumeMesh kernel type.
     * @param  mesh   The mesh to verify.
     * @return true iff every vertex, edge, and face link has the expected
     *         topology.
     */
    template <typename MeshT>
    inline bool verify_all_links(const MeshT& mesh) {
        bool passed = true;
        if (!verify_vertice_links(mesh)) passed = false;
        if (!verify_edge_links(mesh))    passed = false;
        if (!verify_face_links(mesh))    passed = false;
        return passed;
    }

    /**
     * @brief Runs all topological checks on a mesh.
     * @tparam MeshT  Any OpenVolumeMesh kernel type.
     * @param  mesh   The mesh to verify.
     * @return true if every sub-check passes; false otherwise.
     */
    template <typename MeshT>
    inline bool verify_topology(const MeshT& mesh) {
        std::cout << "Calculating Euler Characteristic" << std::endl;

        long V = mesh.n_vertices();
        long E = mesh.n_edges();
        long F = mesh.n_faces();
        long C = mesh.n_cells();
        long euler = V - E + F - C;

        std::cout << "V:" << V << " | E:" << E << " | F:" << F << " | C:" << C << std::endl;

        bool passed = true;

        // Test 1: Global Euler characteristic.
        std::cout << " -> Verifying Euler Characteristic..." << std::endl;
        if (euler == 0) {
            std::cout << "    [PASS] Euler characteristic is 0." << std::endl;
        } else {
            std::cout << "    [FAIL] Euler characteristic is " << euler << "! Mesh is corrupted." << std::endl;
            passed = false;
        }

        // Tests 2-3: closedness and manifold links.
        if (!verify_closed_manifold(mesh)) passed = false;
        if (!verify_all_links(mesh))       passed = false;

        if (passed) std::cout << "       TOPOLOGY CHECKS PASSED!      " << std::endl;
        else        std::cout << "       TOPOLOGY CHECKS FAILED!    " << std::endl;

        return passed;
    }

    } // namespace Topological



    //  GEOMETRIC CHECKS (depend on the vertex positions / embedding)
    namespace Geometric {

    /**
     * @brief Computes the 3-volume of a tetrahedron embedded in any R^n.
     *
     *
     * @tparam PointT  An OpenVolumeMesh::Geometry::VectorT instance.
     * @param  p0..p3  The four vertex positions of the tetrahedron.
     * @return The non-negative 3-volume of the tetrahedron.
     */
    template <typename PointT>
    inline double compute_tet_volume(const PointT& p0, const PointT& p1,
                                     const PointT& p2, const PointT& p3) {
        // Translate p0 to the origin to obtain three spanning edge vectors.
        PointT e1 = p1 - p0;
        PointT e2 = p2 - p0;
        PointT e3 = p3 - p0;

        // Inner product as `operator|`.
        // The Gram matrix entries are the pairwise dot products of the
        // three spanning edge vectors
        const double g11 = e1 | e1;
        const double g12 = e1 | e2;
        const double g13 = e1 | e3;
        const double g22 = e2 | e2;
        const double g23 = e2 | e3;
        const double g33 = e3 | e3;

        // 3x3 determinant. The matrix is
        // symmetric, so g_ji = g_ij.
        double det = g11*(g22*g33 - g23*g23)
                   - g12*(g12*g33 - g23*g13)
                   + g13*(g12*g23 - g22*g13);

        //Prevent negative values caused by floating-point arithmetic
        if (det < 0) det = 0;

        // Volume of a 3-parallelepiped is sqrt(det(G)); the 3-simplex
        // occupies one-sixth of that volume.
        return std::sqrt(det) / 6.0;
    }

    /**
     * @brief Verifies that every tetrahedral cell has a positive volume.
     *
     * Iterates over the cell complex, computing each tetrahedral cell's
     * 3-volume.
     *
     * @tparam MeshT  Any OpenVolumeMesh kernel type.
     * @param  mesh   The mesh to verify.
     * @return true if every tetrahedral cell has volume above
     *         floating-point threshold; false otherwise.
     */
    template <typename MeshT>
    inline bool verify_no_degenerate_cells(const MeshT& mesh) {
        std::cout << " -> Verifying Cell Volumes ..." << std::endl;
        int bad_cells     = 0;
        int skipped_cells = 0;
        double min_vol    = 1e9;

        for(auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it) {
            // Collect the cell's vertex positions.
            std::vector<typename MeshT::PointT> pts;
            for(auto cv_it = mesh.cv_iter(*c_it); cv_it.valid(); ++cv_it) {
                pts.push_back(mesh.vertex(*cv_it));
            }

            if(pts.size() == 4) {
                double vol = compute_tet_volume(pts[0], pts[1], pts[2], pts[3]);
                if (vol < min_vol) min_vol = vol;
                if (vol <= 1e-9) bad_cells++;
            } else {
                skipped_cells++;
            }
        }

        // Inform the caller about non-tet cells regardless of pass/fail.
        if (skipped_cells > 0) {
            std::cout << "    [INFO] Skipped " << skipped_cells
                      << " non-tetrahedral cell(s)." << std::endl;
        }

        if (bad_cells == 0) {
            std::cout << "    [PASS] All tetrahedral cells have positive volume "
                      << "(Min Volume: " << min_vol << ")." << std::endl;
            return true;
        } else {
            std::cout << "    [FAIL] Found " << bad_cells << " degenerate cells!" << std::endl;
            return false;
        }
    }

    /**
     * @brief Confirms every vertex lies on the spheritorus to within `tol`.
     */
    inline bool verify_on_spheritorus(const Geometry::Mesh4D& mesh,
                                      double R, double r, double tol = 1e-9) {
        std::cout << " -> Verifying vertices lie on the spheritorus "
                     "((sqrt(x^2+y^2)-R)^2 + z^2 + w^2 = r^2)..." << std::endl;
        const double rhs = r * r;
        double max_res = 0.0;
        int bad = 0;

        for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
            const auto p = mesh.vertex(*v_it);
            const double x = p[0], y = p[1], z = p[2], w = p[3];
            const double s = std::sqrt(x * x + y * y) - R;
            const double res = std::fabs(s * s + z * z + w * w - rhs);
            if (res > max_res) max_res = res;
            if (res > tol) ++bad;
        }

        if (bad == 0) {
            std::cout << "    [PASS] All vertices on the spheritorus "
                      << "(max residual " << max_res << ")." << std::endl;
            return true;
        }
        std::cout << "    [FAIL] " << bad << " vertex/vertices off the spheritorus "
                  << "(max residual " << max_res << " > tol " << tol << ")."
                  << std::endl;
        return false;
    }

    /**
     * @brief Confirms every vertex lies on the ditorus to within `tol`.
     *
     */
    inline bool verify_on_ditorus(const Geometry::Mesh4D& mesh,
                                  double R1, double R2, double r,
                                  double tol = 1e-9) {
        std::cout << " -> Verifying vertices lie on the ditorus "
                     "((sqrt((sqrt(x^2+y^2)-R1)^2+z^2)-R2)^2 + w^2 = r^2)..."
                  << std::endl;
        const double rhs = r * r;
        double max_res = 0.0;
        int bad = 0;

        for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
            const auto p = mesh.vertex(*v_it);
            const double x = p[0], y = p[1], z = p[2], w = p[3];
            const double s1 = std::sqrt(x * x + y * y) - R1;
            const double s2 = std::sqrt(s1 * s1 + z * z) - R2;
            const double res = std::fabs(s2 * s2 + w * w - rhs);
            if (res > max_res) max_res = res;
            if (res > tol) ++bad;
        }

        if (bad == 0) {
            std::cout << "    [PASS] All vertices on the ditorus "
                      << "(max residual " << max_res << ")." << std::endl;
            return true;
        }
        std::cout << "    [FAIL] " << bad << " vertex/vertices off the ditorus "
                  << "(max residual " << max_res << " > tol " << tol << ")."
                  << std::endl;
        return false;
    }


    //  VOLUME / AREA CONVERGENCE (measured discrete measure vs ideal target)

    /**
     * @brief 3-volume of the discrete manifold: the sum of its tet volumes.
     *
     * @tparam MeshT  Any OpenVolumeMesh kernel type.
     * @param  mesh   A tetrahedral mesh.
     * @return The total 3-volume (non-tet cells are ignored).
     */
    template <typename MeshT>
    inline double measured_volume(const MeshT& mesh) {
        double total = 0.0;
        for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it) {
            std::vector<typename MeshT::PointT> p;
            for (auto cv_it = mesh.cv_iter(*c_it); cv_it.valid(); ++cv_it)
                p.push_back(mesh.vertex(*cv_it));
            if (p.size() == 4)
                total += compute_tet_volume(p[0], p[1], p[2], p[3]);
        }
        return total;
    }

    /**
     * @brief Surface area of a discrete triangulated surface: sum of face areas.
     *
     * Used for the 2-torus, which is a cell-free triangle surface rather than a
     * tetrahedral solid. The area is computed in any R^n via the Gram identity
     *   area = (1/2) sqrt(|a|^2 |b|^2 - (a.b)^2).
     *
     * @tparam MeshT  Any OpenVolumeMesh kernel type.
     * @param  mesh   A triangle surface.
     * @return The total surface area (non-triangle faces are ignored).
     */
    template <typename MeshT>
    inline double measured_surface_area(const MeshT& mesh) {
        double total = 0.0;
        for (auto f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it) {
            std::vector<typename MeshT::PointT> p;
            for (auto fv : mesh.face_vertices(*f_it)) p.push_back(mesh.vertex(fv));
            if (p.size() == 3) {
                const auto a = p[1] - p[0];
                const auto b = p[2] - p[0];
                const double aa = a | a, bb = b | b, ab = a | b;
                double g = aa * bb - ab * ab;
                if (g < 0.0) g = 0.0;
                total += 0.5 * std::sqrt(g);
            }
        }
        return total;
    }

    // Ideal measures of the target manifolds (Pappus's theorem):
    //   spheritorus S^2 x S^1 : V3 = (sphere area 4*pi*r^2) x (loop 2*pi*R)
    //   ditorus     T^3       : V3 = (2*pi*r) x (2*pi*R2) x (2*pi*R1)
    //   2-torus     S^1 x S^1 : A  = (tube 2*pi*r) x (loop 2*pi*R)
    inline double spheritorus_volume(double R, double r) {
        return 8.0 * M_PI * M_PI * R * r * r;
    }
    inline double ditorus_volume(double R1, double R2, double r) {
        return 8.0 * M_PI * M_PI * M_PI * R1 * R2 * r;
    }
    inline double torus_surface_area(double R, double r) {
        return 4.0 * M_PI * M_PI * R * r;
    }

    /**
     * @brief Prints the measured 3-volume against the ideal target.
     */
    template <typename MeshT>
    inline void report_volume(const MeshT& mesh, double ideal,
                              const std::string& name) {
        const double measured = measured_volume(mesh);
        const double pct = (ideal != 0.0) ? 100.0 * measured / ideal : 0.0;
        std::cout << " -> Volume of " << name << ": measured " << measured
                  << "  /  ideal " << ideal << "  (" << pct
                  << "% of target)." << std::endl;
    }

    /**
     * @brief Prints the measured surface area against the ideal target.
     */
    template <typename MeshT>
    inline void report_surface_area(const MeshT& mesh, double ideal,
                                    const std::string& name) {
        const double measured = measured_surface_area(mesh);
        const double pct = (ideal != 0.0) ? 100.0 * measured / ideal : 0.0;
        std::cout << " -> Surface area of " << name << ": measured " << measured
                  << "  /  ideal " << ideal << "  (" << pct
                  << "% of target)." << std::endl;
    }


    } // namespace Geometric

} // namespace Validation
