#pragma once

/**
 * @file quality_metrics.hh
 * @brief Mesh-quality metrics for tetrahedral meshes, exported as CSV.
 *
 * Four distributions are written, one value per line:
 *
 *   <prefix>_volumes.csv      : the 3-volume of every tetrahedral cell.
 *   <prefix>_mean_ratio.csv   : the mean-ratio shape indicator of every cell.
 *   <prefix>_degrees.csv      : the degree (incident-edge count) of every vertex.
 *   <prefix>_dihedrals.csv    : the six dihedral angles (degrees) of every tet.
 *
 * A <prefix>_summary.csv with count / min / max / mean per metric is also
 * written.
 */

#include "validation.hh"
#include <OpenVolumeMesh/Core/Handles.hh>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

namespace QualityMetrics {

    /**
     * @brief The six interior dihedral angles (in degrees) of a tetrahedron.
     */
    template <typename PointT>
    inline std::array<double, 6> tet_dihedral_angles(
        const PointT& p0, const PointT& p1, const PointT& p2, const PointT& p3) {

        const PointT p[4] = {p0, p1, p2, p3};
        std::array<double, 6> angles{};
        int idx = 0;

        for (int i = 0; i < 4; ++i) {
            for (int j = i + 1; j < 4; ++j) {
                // The two vertices not on this edge.
                int k = -1, l = -1;
                for (int m = 0; m < 4; ++m) {
                    if (m == i || m == j) continue;
                    if (k < 0) k = m; else l = m;
                }

                const PointT e = p[j] - p[i];
                const double e2 = e | e;

                PointT ak = p[k] - p[i];
                PointT al = p[l] - p[i];
                if (e2 > 1e-24) {            // strip the component along the edge
                    ak = ak - e * ((ak | e) / e2);
                    al = al - e * ((al | e) / e2);
                }

                const double na = std::sqrt(ak | ak);
                const double nb = std::sqrt(al | al);
                double c = (na > 1e-12 && nb > 1e-12) ? (ak | al) / (na * nb) : 0.0;
                c = std::max(-1.0, std::min(1.0, c));
                angles[idx++] = std::acos(c) * 180.0 / M_PI;
            }
        }
        return angles;
    }

    /**
     * @brief Mean-ratio shape quality of a tetrahedron (1 = regular, 0 = degenerate).
     */
    template <typename PointT>
    inline double tet_mean_ratio(const PointT& p0, const PointT& p1,
                                 const PointT& p2, const PointT& p3) {
        const double V = Validation::Geometric::compute_tet_volume(p0, p1, p2, p3);
        if (V <= 0.0) return 0.0;

        const PointT p[4] = {p0, p1, p2, p3};
        double sum_sq = 0.0;
        for (int i = 0; i < 4; ++i)
            for (int j = i + 1; j < 4; ++j) {
                const PointT e = p[j] - p[i];
                sum_sq += (e | e);          // squared edge length
            }
        if (sum_sq <= 0.0) return 0.0;

        // 12 (3V)^{2/3} = 12 * cbrt((3V)^2) = 12 * cbrt(9 V^2).
        return 12.0 * std::cbrt(9.0 * V * V) / sum_sq;
    }

    /**
     * @brief Count / min / max / mean / median of one distribution.
     */
    struct Summary {
        long   count  = 0;
        double min    = 0.0;
        double max    = 0.0;
        double mean   = 0.0;
        double median = 0.0;
    };

    template <typename T>
    inline Summary summarize(std::vector<T> values) {
        Summary s;
        if (values.empty()) return s;
        std::sort(values.begin(), values.end());
        s.count = static_cast<long>(values.size());
        s.min   = static_cast<double>(values.front());
        s.max   = static_cast<double>(values.back());
        s.mean  = std::accumulate(values.begin(), values.end(), 0.0)
                  / static_cast<double>(s.count);
        const size_t mid = values.size() / 2;
        s.median = (values.size() % 2 == 0)
            ? 0.5 * (static_cast<double>(values[mid - 1]) + static_cast<double>(values[mid]))
            : static_cast<double>(values[mid]);
        return s;
    }

    /**
     * @brief Aggregate summaries returned by analyze_and_export, for callers
     *        driving a parameter study that need the numbers programmatically
     *        rather than by re-parsing the per-configuration CSV files.
     */
    struct QualitySummary {
        Summary volume;
        Summary mean_ratio;
        Summary dihedral_deg;
        Summary vertex_degree;
        double  poor_element_fraction = 0.0;  ///< fraction with mean_ratio < threshold
    };

    namespace detail {

        // Writes one value per line under a single-column header.
        template <typename T>
        inline void write_column(const std::string& path, const std::string& header,
                                 const std::vector<T>& values) {
            std::ofstream f(path);
            if (!f) { std::cerr << "  [warn] could not open " << path << std::endl; return; }
            f << header << "\n";
            for (const auto& v : values) f << v << "\n";
        }

        // Appends a "metric,count,min,max,mean" summary row.
        inline void summary_row(std::ofstream& f, const std::string& name,
                                const Summary& s) {
            if (s.count == 0) { f << name << ",0,,,\n"; return; }
            f << name << "," << s.count << "," << s.min << "," << s.max << "," << s.mean << "\n";
        }

    } // namespace detail

    /**
     * @brief Computes the four quality distributions and writes them as CSV.
     *
     * @tparam MeshT   Any OpenVolumeMesh kernel type.
     * @param  mesh    A tetrahedral mesh (non-tet cells are skipped).
     * @param  prefix  Output path prefix; the parent directory is created if
     *                 needed. Files <prefix>_{volumes,mean_ratio,degrees,
     *                 dihedrals,summary}.csv are produced.
     */
    template <typename MeshT>
    inline QualitySummary analyze_and_export(const MeshT& mesh, const std::string& prefix,
                                             double poor_threshold = 0.1) {
        using PointT = typename MeshT::PointT;

        std::vector<double> volumes;       // one per tet
        std::vector<double> mean_ratios;   // one per tet
        std::vector<double> dihedrals;     // six per tet
        std::vector<int>    degrees;       // one per vertex

        long poor_count = 0;

        for (auto c_it = mesh.cells_begin(); c_it != mesh.cells_end(); ++c_it) {
            std::vector<PointT> p;
            for (auto cv_it = mesh.cv_iter(*c_it); cv_it.valid(); ++cv_it)
                p.push_back(mesh.vertex(*cv_it));
            if (p.size() != 4) continue;

            volumes.push_back(
                Validation::Geometric::compute_tet_volume(p[0], p[1], p[2], p[3]));

            const double eta = tet_mean_ratio(p[0], p[1], p[2], p[3]);
            mean_ratios.push_back(eta);
            if (eta < poor_threshold) ++poor_count;

            const std::array<double, 6> dih =
                tet_dihedral_angles(p[0], p[1], p[2], p[3]);
            for (double a : dih) dihedrals.push_back(a);
        }

        for (auto v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it) {
            int deg = 0;
            for (auto ve_it = mesh.ve_iter(*v_it); ve_it.valid(); ++ve_it) ++deg;
            degrees.push_back(deg);
        }

        const std::filesystem::path parent =
            std::filesystem::path(prefix).parent_path();
        if (!parent.empty()) std::filesystem::create_directories(parent);

        detail::write_column(prefix + "_volumes.csv",      "TetVolume",       volumes);
        detail::write_column(prefix + "_mean_ratio.csv",   "MeanRatio",       mean_ratios);
        detail::write_column(prefix + "_degrees.csv",      "VertexDegree",    degrees);
        detail::write_column(prefix + "_dihedrals.csv",    "DihedralAngleDeg",dihedrals);

        QualitySummary qs;
        qs.volume          = summarize(volumes);
        qs.mean_ratio       = summarize(mean_ratios);
        qs.dihedral_deg     = summarize(dihedrals);
        qs.vertex_degree    = summarize(std::vector<double>(degrees.begin(), degrees.end()));
        qs.poor_element_fraction =
            mean_ratios.empty() ? 0.0
                : static_cast<double>(poor_count) / static_cast<double>(mean_ratios.size());

        std::ofstream s(prefix + "_summary.csv");
        if (s) {
            s << "Metric,Count,Min,Max,Mean\n";
            detail::summary_row(s, "TetVolume",        qs.volume);
            detail::summary_row(s, "MeanRatio",        qs.mean_ratio);
            detail::summary_row(s, "VertexDegree",     qs.vertex_degree);
            detail::summary_row(s, "DihedralAngleDeg", qs.dihedral_deg);
        }

        const auto abs = std::filesystem::absolute(prefix);
        std::cout << "Quality metrics written to " << abs.parent_path().string()
                  << "/  (" << abs.filename().string()
                  << "_{volumes,mean_ratio,degrees,dihedrals,summary}.csv; "
                  << volumes.size() << " tets, " << degrees.size()
                  << " vertices)." << std::endl;

        return qs;
    }
} // namespace QualityMetrics
