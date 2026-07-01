#pragma once

/**
 * @file base_shape_generator.hh
 * @brief Creation of elementary topological sub-complexes used as
 *        Cartesian-product factors.
 *
 *  This header provides three sub-complexes:
 *
 *   - circle(n)       : S^1 as an n-gon on the unit circle.
 *
 *   - interval(L)     : I = [0, L]
 *
 *   - octahedron(r)   : S^2 as a minimal triangulated 2-sphere (6 vertices,
 *                       8 faces)
 *
 * Every function returns a Mesh4D so that products of any pair of sub-complexes
 * can re-use the same underlying type. Coordinates not intrinsic to the
 * factor's dimension are zero-padded (e.g. S^1 leaves Z = W = 0).
 */

#include "ovm_nd_types.hh"
#include <cmath>
#include <vector>
#include <cassert>

namespace ShapeGen {

    using Mesh4D = Geometry::Mesh4D;
    using Vec4d  = Geometry::Vec4d;
    using VH     = OpenVolumeMesh::VertexHandle;

    /**
     * @brief Generates a discretised 1-sphere (S^1) as a regular n-gon.
     *
     * Vertices are placed on the unit circle of the X/Y plane at uniform
     * angular intervals theta_i = 2*pi*i/n, creating n vertices and n edges
     * forming a closed cycle.

     * @param n  Number of vertices (must be at least 3 to form a valid cycle).
     * @return   A Mesh4D containing the S^1 sub-complex with Z = W = 0.
     */
    inline Mesh4D circle(int n) {

        // At least three distinct vertices needed to form a
        // non-self-intersecting closed polygonal cycle.
        assert(n >= 3 && "circle: need at least 3 segments");

        Mesh4D c;

        // 1. Vertex Placement
        // Distribute n vertices uniformly along the unit circle in the X/Y plane
        for (int i = 0; i < n; ++i) {
            const double theta = 2.0 * M_PI * static_cast<double>(i) / n;
            c.add_vertex(Vec4d(std::cos(theta), std::sin(theta), 0.0, 0.0));
        }

        // 2. Edge Connectivity
        // Connect each vertex to its successor modulo n
        for (int i = 0; i < n; ++i) {
            c.add_edge(VH(i), VH((i + 1) % n));
        }

        return c;
    }

    /**
     * @brief Generates a discretised closed interval I = [0, L].
     *
     * @param length  The Euclidean length L of the interval along the X axis.
     * @return        A Mesh4D containing two vertices and one edge.
     */
    inline Mesh4D interval(double length) {
        Mesh4D c;

        // Place the two endpoints along the X axis at 0 and L.
        auto v0 = c.add_vertex(Vec4d(0.0,    0.0, 0.0, 0.0));
        auto v1 = c.add_vertex(Vec4d(length, 0.0, 0.0, 0.0));

        // 1-cell connecting them.
        c.add_edge(v0, v1);

        return c;
    }

    /**
     * @brief Generates a discretised 2-sphere (S^2) as a regular octahedron.
     *
     * The octahedron is the minimal triangulation of S^2 that respects axis
     * alignment: 6 vertices, 12 edges, 8 triangular faces.
     *
     * Face windings are oriented counter-clockwise as seen from outside the
     * sphere.
     *
     * @param r  Radius of the sphere.
     * @return   A Mesh4D containing the octahedral S^2 with W = 0.
     */
    inline Mesh4D octahedron(double r) {
        Mesh4D c;
        // 1. Vertex Placement
        VH v0 = c.add_vertex(Vec4d( r,  0,  0, 0));   // +X pole
        VH v1 = c.add_vertex(Vec4d(-r,  0,  0, 0));   // -X pole
        VH v2 = c.add_vertex(Vec4d( 0,  r,  0, 0));   // +Y pole
        VH v3 = c.add_vertex(Vec4d( 0, -r,  0, 0));   // -Y pole
        VH v4 = c.add_vertex(Vec4d( 0,  0,  r, 0));   // +Z pole (apex)
        VH v5 = c.add_vertex(Vec4d( 0,  0, -r, 0));   // -Z pole (nadir)

        // 2. Face Generation (Edges are inferred by OpenVolumeMesh)

        // Upper cap: triangles around the +Z apex
        c.add_face({v0, v2, v4});   // (+X, +Y, +Z)
        c.add_face({v2, v1, v4});   // (-X, +Y, +Z)
        c.add_face({v1, v3, v4});   // (-X, -Y, +Z)
        c.add_face({v3, v0, v4});   // (+X, -Y, +Z)

        // Lower cap: triangles around the -Z nadir, with reversed winding
        c.add_face({v2, v0, v5});   // (+X, +Y, -Z)
        c.add_face({v1, v2, v5});   // (-X, +Y, -Z)
        c.add_face({v3, v1, v5});   // (-X, -Y, -Z)
        c.add_face({v0, v3, v5});   // (+X, -Y, -Z)

        return c;
    }

} // namespace ShapeGen
