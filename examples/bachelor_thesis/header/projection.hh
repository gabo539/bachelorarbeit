#pragma once

/**
 * @file projection.hh
 * @brief Final embedding of a product complex into R^4.
 *
 * The matching on-surface checks (verify_on_spheritorus / verify_on_ditorus)
 * live in validation.hh under Validation::Geometric.
 */

#include "ovm_nd_types.hh"
#include "nd_subdivision.hh"
#include <OpenVolumeMesh/Core/Handles.hh>
#include <cmath>
#include <iostream>
#include <map>
#include <stdexcept>
#include <vector>

namespace Projection {

    /**
     * @brief Rebuilds a pure-tetrahedral mesh in a new ambient space, remapping
     *        every vertex coordinate through a caller-supplied functor.
     *
     * @tparam MeshOut  Target mesh kernel (e.g. Geometry::Mesh4D).
     * @tparam MeshIn   Source mesh kernel (e.g. Geometry::Mesh5D / Mesh6D).
     * @tparam CoordFn  Callable (MeshIn::PointT) -> MeshOut::PointT.
     * @param  src      Pure-tet source mesh.
     * @param  coord    Coordinate map applied to every vertex.
     * @return The remapped mesh in MeshOut's ambient space.
     * @throws std::runtime_error if any cell is not a tetrahedron.
     */
    template <typename MeshOut, typename MeshIn, typename CoordFn>
    inline MeshOut remap_coordinates(const MeshIn& src, CoordFn coord) {
        using VH  = OpenVolumeMesh::VertexHandle;
        using FH  = OpenVolumeMesh::FaceHandle;
        using HFH = OpenVolumeMesh::HalfFaceHandle;

        MeshOut out;

        // Carry over vertices in handle order so that vmap[i].idx() == i
        std::vector<VH> vmap(src.n_vertices());
        for (auto v_it = src.vertices_begin(); v_it != src.vertices_end(); ++v_it) {
            vmap[v_it->idx()] = out.add_vertex(coord(src.vertex(*v_it)));
        }

        std::map<std::vector<int>, FH> face_map;

        for (auto c_it = src.cells_begin(); c_it != src.cells_end(); ++c_it) {
            std::vector<VH> cv;
            for (auto cv_it = src.cv_iter(*c_it); cv_it.valid(); ++cv_it) {
                cv.push_back(vmap[cv_it->idx()]);
            }
            if (cv.size() != 4) {
                throw std::runtime_error(
                    "remap_coordinates: expected a pure-tetrahedral mesh (call "
                    "Subdivision::triangulate_prisms_to_tets first).");
            }
            Subdivision::add_tet_to_mesh<MeshOut, VH, FH, HFH>(
                out, cv[0], cv[1], cv[2], cv[3], face_map);
        }

        return out;
    }

    /**
     * @brief Projects the R^5 product S^2 x S^1 onto the spheritorus in R^4.
     *
     * Reads the sphere block from coordinates [0,1,2] = (a,b,c) and the circle
     * block from [3,4] = (d,f), then applies
     *   x = (R + r a) d,  y = (R + r a) f,  z = r b,  w = r c.
     *
     * @param src  Pure-tet S^2 x S^1 mesh whose factor blocks are unit-length
     *             (run Subdivision::normalize_blocks(mesh, {3,2}) after any
     *             subdivision so this holds).
     * @param R    Major radius (tube-centre offset); require R > r > 0.
     * @param r    Minor radius (tube thickness).
     */
    inline Geometry::Mesh4D project_spheritorus(const Geometry::Mesh5D& src,
                                                double R, double r) {
        std::cout << "Projecting S^2 x S^1 (R^5) -> spheritorus (R^4)  [R=" << R
                  << ", r=" << r << "] ..." << std::endl;
        return remap_coordinates<Geometry::Mesh4D>(
            src, [R, r](const Geometry::Vec5d& p) {
                const double a = p[0], b = p[1], c = p[2], d = p[3], f = p[4];
                const double rho = R + r * a;        // distance to the z=w=0 axis
                return Geometry::Vec4d(rho * d, rho * f, r * b, r * c);
            });
    }

    /**
     * @brief Projects the R^6 product S^1 x S^1 x S^1 onto the ditorus in R^4.
     *
     * Reads the inner circle from [0,1] = (a,b), the middle circle from
     * [2,3] = (c,d), and the outer circle from [4,5] = (e,f), then applies
     *   x = (R1 + (R2 + r a) c) e,  y = (R1 + (R2 + r a) c) f,
     *   z = (R2 + r a) d,           w = r b.
     *
     * @param src  Pure-tet 3-torus mesh whose factor blocks are unit-length
     *             (run Subdivision::normalize_blocks(mesh, {2,2,2}) after any
     *             subdivision so this holds).
     * @param R1   Outer radius; require R1 > R2 + r.
     * @param R2   Middle radius; require R2 > r.
     * @param r    Inner tube radius.
     */
    inline Geometry::Mesh4D project_ditorus(const Geometry::Mesh6D& src,
                                            double R1, double R2, double r) {
        std::cout << "Projecting S^1 x S^1 x S^1 (R^6) -> ditorus (R^4)  [R1="
                  << R1 << ", R2=" << R2 << ", r=" << r << "] ..." << std::endl;
        return remap_coordinates<Geometry::Mesh4D>(
            src, [R1, R2, r](const Geometry::Vec6d& p) {
                const double a = p[0], b = p[1];   // inner circle  (radius r)
                const double c = p[2], d = p[3];   // middle circle (radius R2)
                const double e = p[4], f = p[5];   // outer circle  (radius R1)
                const double mid   = R2 + r * a;
                const double outer = R1 + mid * c;
                return Geometry::Vec4d(outer * e, outer * f, mid * d, r * b);
            });
    }

} // namespace Projection
