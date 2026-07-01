#pragma once

/**
 * @file visualize_4d.hh
 * @brief Polyscope view of a 4D tetrahedral mesh.
 *
 * A closed 3-manifold in R^4 cannot be drawn directly, so one coordinate axis
 * is colour-scaled
 *
 * The solid is rendered the same way as generate_prism: we emit the oriented
 * boundary triangles of every cell (its half-faces, walked in their stored
 * winding). Back-face culling is enabled. Toggle the back-face option in the polyscope UI to
 * compare. Every interior face is intentionally emitted twice.
 */

#include "ovm_nd_types.hh"
#include <OpenVolumeMesh/Core/Handles.hh>

#include <array>
#include <iostream>
#include <string>
#include <vector>

#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/curve_network.h"
#include "polyscope/slice_plane.h"

namespace Visualize {

    /**
     * @brief Shows a 4D tetrahedral mesh in polyscope.
     *
     * @tparam MeshT      An OpenVolumeMesh kernel with 4D points.
     * @param  mesh       The tetrahedral mesh to display (already oriented).
     * @param  drop_axis  Axis 0..3 used as the colour scalar; the other three
     *                    become the 3D position.
     * @param  name       Display name for the registered structures.
     */
    template <typename MeshT>
    inline void show_4d(const MeshT& mesh, int drop_axis, const std::string& name) {
        std::cout << "\n--- Launching Polyscope (" << name << ", axis "
                  << drop_axis << " as colour) ---" << std::endl;
        polyscope::init();

        // 3D positions (the three kept axes) and the dropped axis as colour.
        std::vector<std::array<double, 3>> pos;
        std::vector<double> color;
        pos.reserve(mesh.n_vertices());
        color.reserve(mesh.n_vertices());

        std::vector<int> vmap(mesh.n_vertices());
        int idx = 0;
        for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
            const auto p = mesh.vertex(*v_it);
            double xyz[3]; int k = 0;
            for (int i = 0; i < 4; ++i) if (i != drop_axis) xyz[k++] = p[i];
            pos.push_back({xyz[0], xyz[1], xyz[2]});
            color.push_back(p[drop_axis]);
            vmap[v_it->idx()] = idx++;
        }


        std::vector<std::vector<size_t>> faces;
        int skipped = 0;
        for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it) {
            for (auto chf_it = mesh.chf_iter(*c_it); chf_it.valid(); ++chf_it) {
                std::vector<size_t> tri;
                for (auto hfhe_it = mesh.hfhe_iter(*chf_it); hfhe_it.valid(); ++hfhe_it)
                    tri.push_back(static_cast<size_t>(
                        vmap[mesh.from_vertex_handle(*hfhe_it).idx()]));
                if (tri.size() != 3) { ++skipped; continue; }
                faces.push_back(std::move(tri));
            }
        }
        std::cout << "  Emitted " << faces.size() << " oriented boundary triangles"
                  << (skipped ? " (skipped " + std::to_string(skipped)
                                + " non-triangular)" : "") << "." << std::endl;

        auto* sm = polyscope::registerSurfaceMesh(name, pos, faces);
        sm->addVertexScalarQuantity("dropped coordinate", color)->setEnabled(true);
        sm->setBackFacePolicy(polyscope::BackFacePolicy::Cull);   // toggle in UI

        // Wireframe skeleton, also coloured by the dropped coordinate.
        std::vector<std::array<size_t, 2>> edges;
        edges.reserve(mesh.n_edges());
        for (auto e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it) {
            edges.push_back({static_cast<size_t>(vmap[mesh.edge(*e_it).from_vertex().idx()]),
                             static_cast<size_t>(vmap[mesh.edge(*e_it).to_vertex().idx()])});
        }
        auto* cn = polyscope::registerCurveNetwork(name + " skeleton", pos, edges);
        cn->addNodeScalarQuantity("dropped coordinate", color)->setEnabled(true);
        cn->setRadius(0.003);
        cn->setEnabled(false);   // off by default; the surface is the main view

        polyscope::addSceneSlicePlane();
        polyscope::show();
    }

} // namespace Visualize
