#pragma once

/**
 * @file ovm_nd_types.hh
 * @brief Centralised mesh type aliases.
 *
 * The construction pipeline operates in mixed dimensions: base sub-complexes
 * (S^1, S^2, I, ...) are stored in a 4D "cargo" mesh whose unused coordinates
 * are zero-padded, while their Cartesian products may live in R^5 or R^6
 * before being projected down to R^4. This header centralises the
 * vector and mesh kernel typedefs so every translation unit operates on
 * exactly the same instantiation of OpenVolumeMesh::GeometryKernel.
 *
 * The Mesh3D alias is provided as an output target for low-dimensional
 * algorithm-validation runs (S^1 x I, S^1 x S^1, etc.) where the natural
 * ambient space is R^3 and the result is directly visualizable.
 */

#include <OpenVolumeMesh/Geometry/VectorT.hh>
#include <OpenVolumeMesh/Geometry/Vector11T.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>

namespace Geometry {

    // 3D Types
    // ------------------------------------------------------------------------
    // Used for low-dimensional validation runs of the pipeline. Examples:
    //   - S^1 x I = open cylinder (2-complex with quad faces).
    //   - S^1 x S^1 = torus (closed 2-manifold).
    using Vec3d  = OpenVolumeMesh::Geometry::Vec3d;
    using Mesh3D = OpenVolumeMesh::GeometricPolyhedralMeshV3d;


    // 4D Types
    // ------------------------------------------------------------------------
    // Used as the "ambient" type for base sub-complexes. Components beyond the
    // intrinsic dimension of the shape are zero-padded; for instance, S^1 only
    // populates components [0] and [1], leaving [2] and [3] at zero until a
    // Cartesian product re-projects them. Zero-padding into one shared type is
    // a deliberate choice over distinct 1D-3D point types: it keeps the
    // product and combiner code uniform (it just reads the coordinates it
    // needs)
    using Vec4d  = OpenVolumeMesh::Geometry::Vec4d;
    using Mesh4D = OpenVolumeMesh::GeometricPolyhedralMeshV4d;


    // 5D Types
    // ------------------------------------------------------------------------
    // Required as an intermediate ambient space for
    // S^2 x S^1
    using Vec5d  = OpenVolumeMesh::Geometry::Vec5d;
    using Mesh5D = OpenVolumeMesh::GeometryKernel<Vec5d, OpenVolumeMesh::TopologyKernel>;


    // 6D Types (OpenVolumeMesh-native)
    // ------------------------------------------------------------------------
    // Required as an intermediate ambient space for
    // S^1 x S^1 x S^1
    using Vec6d  = OpenVolumeMesh::Geometry::Vec6d;
    using Mesh6D = OpenVolumeMesh::GeometryKernel<Vec6d, OpenVolumeMesh::TopologyKernel>;

} // namespace Geometry
