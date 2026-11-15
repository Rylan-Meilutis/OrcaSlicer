#include "ThermoelasticSections.hpp"
#include "ThermoelasticExactSections.hpp"
#include "SectionIntegration.hpp"
#include "SectionGrowthStream.hpp"
#include "ThermoelasticThermalStrain.hpp"
#include "ThermoelasticSectionContact.hpp"
#include "ThermoelasticSparseFactor.hpp"
#include "BoundingBox.hpp"
#include "ClipperUtils.hpp"
#include "ExPolygon.hpp"
#include "Point.hpp"
#include "Polygon.hpp"
#include "ThermoelasticLayers.hpp"
#include "ThermoelasticMaterial.hpp"
#include "ThermoelasticMesh.hpp"
#include "ThermoelasticSectionTopology.hpp"
#include "ThermoelasticLocalSections.hpp"
#include "ThermoelasticLocalTopology.hpp"
#include "libslic3r.h"
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <Eigen/OrderingMethods>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

namespace Slic3r::Thermoelastic {
namespace {
using Index    = std::array<int, 4>; // XYZ plus physical component identity.
using Matrix24 = Eigen::Matrix<double, 24, 24>;
using Vector24 = Eigen::Matrix<double, 24, 1>;
using BMatrix  = Eigen::Matrix<double, 6, 24>;
Eigen::Matrix3d strain_tensor(const Eigen::Matrix<double,6,1>& strain)
{
    Eigen::Matrix3d tensor;
    tensor << strain[0],.5*strain[3],.5*strain[5],
              .5*strain[3],strain[1],.5*strain[4],
              .5*strain[5],.5*strain[4],strain[2];
    return tensor;
}
Eigen::Vector3d compatible_displacement(const SectionMaterial& material,
    const Eigen::Vector3d& center, const Eigen::Vector3d& point)
{
    const Eigen::Vector3d r = point-center;
    Eigen::Vector3d u = strain_tensor(material.eigenstrain+material.eigenstrain_gradient_per_m*center)*r;
    std::array<Eigen::Matrix3d,3> derivative;
    for (unsigned k = 0; k < 3; ++k)
        derivative[k] = strain_tensor(material.eigenstrain_gradient_per_m.col(k));
    // u_i,jk = epsilon_ij,k + epsilon_ik,j - epsilon_jk,i. A linear
    // symmetric strain has zero compatibility derivatives, so this quadratic
    // displacement is an exact stress-free solution on a genuinely free body.
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j)
            for (unsigned k = 0; k < 3; ++k)
                u[i] += .5*(derivative[k](i,j)+derivative[j](i,k)-derivative[i](j,k))*r[j]*r[k];
    if (!u.allFinite()) throw std::runtime_error("Nonfinite compatible section displacement.");
    return u;
}
struct Cell
{
    Matrix24 stiffness    = Matrix24::Zero();
    Vector24 force        = Vector24::Zero();
    double initial_energy = 0., volume = 0., bed_area = 0.;
    std::array<std::size_t, 8> nodes;
    ExPolygons bed_trace;
    std::vector<ExactSection> exact_bed_traces;
};
struct CellMaterialMoments {
    std::array<double, 27> integral{};
    std::array<coord_t,3> anchor{}, minimum{}, maximum{};
    bool initialized=false;
    double bed_area = 0.;
    ExPolygons bed_trace;
    std::vector<ExactSection> exact_bed_traces;
};
std::size_t point_count(const ExPolygons& polygons)
{
    std::size_t count = 0;
    for (const auto& polygon : polygons) {
        count += polygon.contour.points.size();
        for (const auto& hole : polygon.holes)
            count += hole.points.size();
    }
    return count;
}
std::array<BMatrix, 8> strain_polynomials(const Eigen::Vector3d& spacing)
{
    std::array<BMatrix, 8> result;
    for (auto& value : result)
        value.setZero();
    for (unsigned node = 0; node < 8; ++node)
        for (unsigned axis = 0; axis < 3; ++axis) {
            for (unsigned monomial = 0; monomial < 8; ++monomial) {
                if (monomial & (1u << axis))
                    continue;
                double coefficient = (node & (1u << axis) ? 1. : -1.) / spacing[axis];
                for (unsigned other = 0; other < 3; ++other)
                    if (other != axis) {
                        const bool high = node & (1u << other), power = monomial & (1u << other);
                        coefficient *= high ? (power ? 1. : 0.) : (power ? -1. : 1.);
                    }
                for (unsigned component = 0; component < 3; ++component) {
                    unsigned row = axis;
                    if (component != axis)
                        row = (axis + component == 1) ? 3 : ((axis + component == 3) ? 4 : 5);
                    result[monomial](row, 3 * node + component) += coefficient;
                }
            }
        }
    return result;
}
unsigned moment_index(unsigned a, unsigned b = 0)
{ return (a & 1) + (b & 1) + 3 * (((a >> 1) & 1) + ((b >> 1) & 1)) + 9 * (((a >> 2) & 1) + ((b >> 2) & 1)); }
using PieceProducer = std::function<std::vector<std::size_t>(const std::function<void(const SectionPiece&)>&)>;
SectionMechanicalResult solve_section_source(double bed_z_m,
                                             const PieceProducer& produce,
                                             const std::vector<SectionMaterial>& materials,
                                             const Eigen::Vector3d& origin,
                                             const Eigen::Vector3d& spacing,
                                             const SectionSolveLimits& limits,
                                             const std::function<void()>& checkpoint,
                                             const std::vector<SectionSampleRequest>& samples,
                                             const StreamIntegrationLimits* integration = nullptr,
                                             StreamIntegrationStatistics* statistics = nullptr,
                                             const std::vector<double>& z_planes = {},
                                             const std::vector<SectionBedReleaseRegion>& released_regions = {},
                                             const std::optional<SectionCellContact>& cell_contact = {},
                                             const SectionLoadField& load = {})
{
    if (statistics) *statistics = {};
    if (cell_contact && (released_regions.empty() || !cell_contact->max_iterations || cell_contact->max_iterations > 16))
        throw std::invalid_argument("Invalid section cell-contact study.");
    const std::size_t grid_multiplier = limits.sparse_factor ? 16 : 1;
    if (limits.sparse_factor && (!limits.sparse_factor->max_entries || limits.sparse_factor->max_entries > 64000000 ||
        !limits.sparse_factor->max_work || limits.sparse_factor->max_work > 256000000))
        throw std::invalid_argument("Invalid section sparse factor limits.");
    const std::size_t maximum_integration_work = limits.sparse_factor ? 3200000000u : 600000000;
    if (integration && (!integration->max_work || integration->max_work > maximum_integration_work ||
        !integration->max_total_points || integration->max_total_points > (limits.sparse_factor ? 64000000u : 32000000u) ||
        !integration->max_live_points || integration->max_live_points > 200000 ||
        !integration->max_groups || integration->max_groups > (limits.sparse_factor ? 125000u : 31250u) ||
        !integration->max_memberships || integration->max_memberships > 2000000 ||
        !integration->max_pieces || integration->max_pieces > (limits.sparse_factor ? 500000u : 250000u) ||
        !integration->max_provisional_components || integration->max_provisional_components > 250000))
        throw std::invalid_argument("Invalid streaming integration limits.");
    if (samples.size() > 200000 || !produce || !std::isfinite(bed_z_m) || std::abs(bed_z_m) > 10. || materials.empty() ||
        materials.size() > 256 || !origin.allFinite() || !spacing.allFinite() || (spacing.array() <= 0.).any() ||
        origin.cwiseAbs().maxCoeff() > 10. || spacing.maxCoeff() > 10. || spacing.minCoeff() < 1e-9 || !limits.max_cells ||
        limits.max_cells > 4000 * grid_multiplier || !limits.max_nodes || limits.max_nodes > 6000 * grid_multiplier ||
        !limits.max_free_dofs || limits.max_free_dofs > 12000 * grid_multiplier ||
        !limits.max_points || limits.max_points > 200000 || !limits.max_work || limits.max_work > 2000000 * grid_multiplier ||
        !limits.max_triplets || limits.max_triplets > 1200000 * grid_multiplier ||
        !limits.max_bed_trace_points || limits.max_bed_trace_points > 2000000)
        throw std::invalid_argument("Invalid section mechanics geometry or limits.");
    std::array<coord_t, 3> native_origin, native_spacing;
    for (unsigned axis = 0; axis < 3; ++axis) {
        // Polygon clipping and the Q1 basis must share exactly the same grid.
        for (double value : {origin[axis], spacing[axis]}) {
            const double native = value * 1e9;
            if (std::abs(native - std::round(native)) > 1e-5)
                throw std::invalid_argument("Section background grid must use native nanometre coordinates.");
        }
        native_origin[axis]  = std::llround(origin[axis] * 1e9);
        native_spacing[axis] = std::llround(spacing[axis] * 1e9);
    }
    if (std::llround(bed_z_m * 1e9) != native_origin[2])
        throw std::invalid_argument("Section grid and topology use different bed planes.");
    if (released_regions.size() > 64)
        throw std::invalid_argument("Too many section bed release regions.");
    std::vector<std::array<coord_t,4>> release_cells;
    for (const auto& region : released_regions) {
        if (!region.minimum_m.allFinite() || !region.maximum_m.allFinite() ||
            region.minimum_m.cwiseAbs().maxCoeff() > 10. || region.maximum_m.cwiseAbs().maxCoeff() > 10. ||
            (region.maximum_m.array() <= region.minimum_m.array()).any())
            throw std::invalid_argument("Invalid section bed release region.");
        std::array<coord_t,4> box;
        for (unsigned axis = 0; axis < 2; ++axis)
            for (unsigned end = 0; end < 2; ++end) {
                const double value = (end ? region.maximum_m[axis] : region.minimum_m[axis]) * 1e9;
                const coord_t native = std::llround(value);
                if (std::abs(value-double(native)) > 1e-5 || (native-native_origin[axis]) % native_spacing[axis] != 0)
                    throw std::invalid_argument("Section bed release regions must align with basis cell faces.");
                box[axis+2*end] = (native-native_origin[axis])/native_spacing[axis];
            }
        release_cells.push_back(box);
    }
    std::vector<coord_t> native_z;
    if (!z_planes.empty()) {
        if (z_planes.size() < 2 || z_planes.size() > 128)
            throw std::invalid_argument("Invalid graded section grid plane count.");
        for (const double plane : z_planes) {
            if (!std::isfinite(plane) || std::abs(plane) > 10. ||
                std::abs(plane * 1e9 - std::round(plane * 1e9)) > 1e-5)
                throw std::invalid_argument("Graded section planes must use native nanometre coordinates.");
            const coord_t value = std::llround(plane * 1e9);
            if (!native_z.empty() && value <= native_z.back())
                throw std::invalid_argument("Graded section planes must increase strictly.");
            native_z.push_back(value);
        }
        if (native_z.front() != native_origin[2])
            throw std::invalid_argument("Graded section grid and topology use different bed planes.");
    }
    const auto z_coordinate = [&](int z) {
        return native_z.empty() ? native_origin[2] + z * native_spacing[2] : native_z.at(std::size_t(z));
    };
    const auto z_size = [&](int z) { return z_coordinate(z+1) - z_coordinate(z); };
    std::size_t work = 0, integration_work = 0, assembly_work = 0;
    bool integrating = true;
    const auto charge = [&](std::size_t amount = 1) {
        if (checkpoint)
            checkpoint();
        auto& consumed = integration ? (integrating ? integration_work : assembly_work) : work;
        const auto maximum = integration && integrating ? integration->max_work : limits.max_work;
        if (amount > maximum - consumed)
            throw ResourceLimit("Section mechanics work budget exceeded: phase=" + std::string(integrating ? "integration" : "assembly") +
                ", used=" + std::to_string(consumed) + ", requested=" + std::to_string(amount) +
                ", limit=" + std::to_string(maximum) + ".");
        consumed += amount;
        if (!integration) (integrating ? integration_work : assembly_work) += amount;
        if (statistics) {
            statistics->integration_work = integration_work;
            statistics->assembly_work = assembly_work;
        }
    };
    const auto live_points = integration ? std::min(limits.max_points, integration->max_live_points) : limits.max_points;
    const auto total_points = integration ? integration->max_total_points : limits.max_points;
    const auto max_groups = integration ? integration->max_groups : limits.max_work / 64;
    const auto max_memberships = integration ? integration->max_memberships : limits.max_work;
    const auto max_pieces = integration ? integration->max_pieces : 10000;
    const auto max_provisional_components = integration ? integration->max_provisional_components : 10000;
    std::vector<SectionMaterial> constitutive = materials;
    for (auto& material : constitutive) {
        material.elasticity_pa = validated_elasticity(material.elasticity_pa);
        if (!material.eigenstrain.allFinite() || !material.eigenstrain_gradient_per_m.allFinite())
            throw std::invalid_argument("Invalid section eigenstrain.");
    }
    // Both opaque topology factories union each interval/material before
    // moving its individual polygons into pieces. Borrow that certified
    // geometry: repeating the union neither adds authority nor needs a copy.
    std::size_t points = 0, piece_count = 0, bed_trace_points = 0;
    const auto validate_piece = [&](const SectionPiece& piece) {
        if (piece_count == max_pieces)
            throw ResourceLimit("Section mechanics processed piece budget exceeded: limit=" + std::to_string(max_pieces) + ".");
        if (piece.component >= max_provisional_components)
            throw ResourceLimit("Section mechanics provisional component budget exceeded: id=" + std::to_string(piece.component) +
                ", limit=" + std::to_string(max_provisional_components) + ".");
        ++piece_count;
        if (statistics) statistics->piece_count = piece_count;
        const auto& slab = piece.slab;
        if (!std::isfinite(slab.bottom_m) || !std::isfinite(slab.top_m) || slab.top_m <= slab.bottom_m || std::abs(slab.top_m) > 10. ||
            slab.material >= materials.size())
            throw std::invalid_argument("Invalid mechanical slab.");
        if (std::llround(slab.bottom_m * 1e9) < native_origin[2])
            throw std::invalid_argument("Mechanical slab is below the bed.");
        const auto n = slab.footprint.vertex_count();
        if (n > live_points)
            throw ResourceLimit("Section mechanics live input point budget exceeded.");
        if (n > total_points - points)
            throw ResourceLimit("Section mechanics input point budget exceeded.");
        points += n;
        if (statistics) {
            statistics->processed_points = points;
            statistics->peak_input_points = std::max(statistics->peak_input_points, n);
        }
        if (!slab.footprint.is_exact()) for (const auto& polygon : slab.footprint.native()) {
            const auto valid = [](const Polygon& ring) {
                if (ring.points.size() < 3)
                    throw std::invalid_argument("Invalid section polygon.");
                for (const auto& p : ring.points)
                    if (std::abs(double(p.x())) > 1e10 || std::abs(double(p.y())) > 1e10)
                        throw std::invalid_argument("Section polygon coordinate outside supported range.");
            };
            valid(polygon.contour);
            for (const auto& hole : polygon.holes)
                valid(hole);
        }
        charge(n + 1);
    };
    const auto regular_polynomial = strain_polynomials(spacing);
    std::map<coord_t, std::array<BMatrix,8>> graded_polynomials;
    const auto polynomial_at = [&](int z) -> const std::array<BMatrix,8>& {
        if (native_z.empty()) return regular_polynomial;
        const auto height = z_size(z);
        const auto found = graded_polynomials.find(height);
        if (found != graded_polynomials.end()) return found->second;
        auto size = spacing;
        size.z() = double(height) * 1e-9;
        return graded_polynomials.emplace(height,strain_polynomials(size)).first->second;
    };
    std::map<Index, Cell> cells;
    // Integrals are additive over disjoint sampled slabs. Keep material and
    // physical-component identity: neither voids nor constitutive interfaces
    // are replaced by a density-scaled cell. Each insertion follows a >=64
    // work charge below, bounding storage by max_work/64 groups.
    std::map<std::pair<Index, std::size_t>, CellMaterialMoments> grouped_moments;
    // Membership is spatial, not a global Z-band lookup: tile-local geometry
    // may have different vertical partitions in neighboring XY columns.
    std::map<std::array<int, 3>, std::vector<std::size_t>> sample_requests;
    std::vector<std::vector<std::size_t>> sample_components(samples.size());
    std::size_t membership_records = 0;
    const auto index               = [&](coord_t coordinate, unsigned axis) {
        if (axis == 2 && !native_z.empty())
            return int(std::upper_bound(native_z.begin(),native_z.end(),coordinate)-native_z.begin())-1;
        const coord_t delta = coordinate - native_origin[axis];
        const coord_t value = delta / native_spacing[axis] - (delta % native_spacing[axis] < 0 ? 1 : 0);
        if (value < -10000000 || value > 10000000)
            throw ResourceLimit("Section background grid range exceeded.");
        return int(value);
    };
    for (std::size_t i = 0; i < samples.size(); ++i) {
        charge();
        const auto& request = samples[i];
        if (!request.position_m.allFinite() || request.position_m.cwiseAbs().maxCoeff() > 10. || request.material >= materials.size())
            throw std::invalid_argument("Invalid section field sample request.");
        sample_requests[{index(std::llround(request.position_m.x() * 1e9), 0), index(std::llround(request.position_m.y() * 1e9), 1),
                         index(std::llround(request.position_m.z() * 1e9), 2)}]
            .push_back(i);
    }
    const auto component_mapping = produce([&](const SectionPiece& piece) {
        validate_piece(piece);
        const auto& footprint = piece.slab.footprint;
        const bool exact = footprint.is_exact();
        const ExPolygons empty_geometry;
        const auto& geometry = exact ? empty_geometry : footprint.native();
        if (footprint.empty())
            return;
        const auto box = *footprint.bounds();
        const auto exact_budget = [&] {
            const auto remaining = integration ? integration->max_work-integration_work : limits.max_work-work;
            if (!remaining) throw ResourceLimit("Exact mechanical integration budget exhausted.");
            return std::min<std::size_t>(4000000,remaining);
        };
        const auto& slab = piece.slab;
        const coord_t slab_bottom = std::llround(slab.bottom_m * 1e9), slab_top = std::llround(slab.top_m * 1e9);
        if (!native_z.empty() && (slab_bottom < native_z.front() || slab_top > native_z.back()))
            throw std::invalid_argument("Graded section grid does not cover the complete solid.");
        for (int x = index(box.min.x(), 0); x <= index(box.max.x(), 0); ++x)
            for (int y = index(box.min.y(), 1); y <= index(box.max.y(), 1); ++y) {
                charge(footprint.vertex_count() + 1);
                const coord_t left = native_origin[0] + x * native_spacing[0], right = left + native_spacing[0];
                const coord_t bottom = native_origin[1] + y * native_spacing[1], top = bottom + native_spacing[1];
                if (right <= left || top <= bottom)
                    throw std::invalid_argument("Unresolved section grid spacing.");
                ExactSection exact_clipped;
                if (exact) {
                    const ExPolygons cell_geometry{ExPolygon(Polygon{{left,bottom},{right,bottom},{right,top},{left,top}})};
                    const auto cell = make_exact_section(cell_geometry,exact_budget(),[&] { charge(0); });
                    charge(cell.work);
                    const auto cut = boolean_exact_sections(footprint.exact(),cell.section,ExactSectionOperation::Intersection,
                        exact_budget(),[&] { charge(0); });
                    charge(cut.work);
                    exact_clipped = cut.section;
                }
                const auto clipped = exact ? ExPolygons{} : intersection_ex(geometry,
                    Polygons{Polygon{{left,bottom},{right,bottom},{right,top},{left,top}}});
                const auto clipped_points = exact ? exact_clipped.vertex_count() : point_count(clipped);
                if (clipped_points > live_points)
                    throw ResourceLimit("Section clip point budget exceeded.");
                if (statistics) statistics->peak_clip_points = std::max(statistics->peak_clip_points, clipped_points);
                if (exact ? exact_clipped.empty() : clipped.empty())
                    continue;
                for (int z = index(slab_bottom, 2); z <= index(slab_top-1, 2); ++z) {
                    charge(64 + 9 * clipped_points);
                    const std::array<coord_t, 3> start{left, bottom, z_coordinate(z)};
                    auto size = native_spacing;
                    size[2] = z_size(z);
                    const coord_t lo = std::max(start[2], slab_bottom), hi = std::min(start[2] + size[2], slab_top);
                    if (hi <= lo)
                        continue;
                    // Outward metadata bounds choose a conservative conditioning
                    // support only; integration still uses the exact clipped set.
                    const auto clipped_box=exact ? *exact_clipped.bounds() : get_extents(clipped);
                    const std::array<coord_t,3> piece_min{clipped_box.min.x(),clipped_box.min.y(),lo};
                    const std::array<coord_t,3> piece_max{clipped_box.max.x(),clipped_box.max.y(),hi};
                    const Index key{x, y, z, int(piece.component)};
                    const auto group_key = std::make_pair(key, slab.material);
                    if (!grouped_moments.count(group_key) && grouped_moments.size() >= max_groups)
                        throw ResourceLimit("Section mechanics moment-group budget exceeded.");
                    auto& group = grouped_moments[group_key];
                    if (!group.initialized) {
                        group.anchor=group.minimum=piece_min;
                        group.maximum=piece_max;
                        group.initialized=true;
                    }
                    for (unsigned axis=0;axis<3;++axis) {
                        group.minimum[axis]=std::min(group.minimum[axis],piece_min[axis]);
                        group.maximum[axis]=std::max(group.maximum[axis],piece_max[axis]);
                    }
                    // Integrate near occupied support, not a distant background
                    // corner. Rebasing later preserves all 27 physical moments.
                    SectionIntegration::Moments moment;
                    if (exact) {
                        // The cell is the clipping domain, not the group anchor.
                        // Translate the resulting polynomial moments afterwards.
                        const auto integrated = exact_section_volume_moments(exact_clipped,start,size,lo,hi,
                            exact_budget(),[&] { charge(0); });
                        charge(integrated.work);
                        Eigen::Vector3d offset;
                        for (unsigned axis=0;axis<3;++axis)
                            offset[axis]=double(start[axis]-group.anchor[axis])/double(size[axis]);
                        moment=SectionIntegration::rebase_moments(integrated.values,offset,Eigen::Vector3d::Ones());
                    } else {
                        moment=SectionIntegration::volume_moments(clipped, group.anchor, size, lo, hi);
                    }
                    if (!(moment[0] > 0.) || !std::isfinite(moment[0]))
                        throw std::invalid_argument("Invalid clipped section volume.");
                    const auto requests = sample_requests.find({x, y, z});
                    if (requests != sample_requests.end())
                        for (const auto i : requests->second) {
                            charge();
                            const auto& request    = samples[i];
                            const coord_t sample_z = std::llround(request.position_m.z() * 1e9);
                            if (request.material != slab.material || sample_z < slab_bottom || sample_z >= slab_top)
                                continue;
                            charge(footprint.vertex_count());
                            const Point xy(static_cast<coord_t>(std::llround(request.position_m.x() * 1e9)),
                                           static_cast<coord_t>(std::llround(request.position_m.y() * 1e9)));
                            bool contained;
                            if (exact) {
                                const auto membership=exact_closed_contains_point(footprint.exact(),xy,exact_budget(),[&] { charge(0); });
                                charge(membership.work);
                                contained=membership.contained && !membership.boundary;
                            } else {
                                contained=std::any_of(geometry.begin(), geometry.end(), [&](const ExPolygon& p) { return p.contains(xy, false); });
                            }
                            if (contained) {
                                if (membership_records == max_memberships)
                                    throw ResourceLimit("Section sample membership budget exceeded.");
                                ++membership_records;
                                if (statistics) statistics->peak_memberships = membership_records;
                                sample_components[i].push_back(piece.component);
                            }
                        }
                    if (statistics) statistics->peak_groups = std::max(statistics->peak_groups, grouped_moments.size());
                    for (std::size_t i = 0; i < moment.size(); ++i)
                        group.integral[i] += moment[i];
                    if (lo == native_origin[2]) {
                        group.bed_area += moment[0] / (double(hi - lo) * 1e-9);
                        if (!released_regions.empty()) {
                            const auto count = exact ? exact_clipped.vertex_count() : point_count(clipped);
                            if (count > limits.max_bed_trace_points-bed_trace_points)
                                throw ResourceLimit("Section bed trace point budget exceeded.");
                            bed_trace_points += count;
                            if (exact) group.exact_bed_traces.push_back(exact_clipped);
                            else group.bed_trace.insert(group.bed_trace.end(),clipped.begin(),clipped.end());
                        }
                    }
                }
            }
    });
    if (statistics) statistics->provisional_mapping_size = component_mapping.size();
    if (component_mapping.size() > max_provisional_components) {
        if (integration) throw ResourceLimit("Section mechanics provisional mapping budget exceeded.");
        throw std::invalid_argument("Missing completed section connectivity.");
    }
    if (!piece_count || component_mapping.empty())
        throw std::invalid_argument("Missing completed section connectivity.");
    // Components may merge only after a later neighboring tile arrives. Keep
    // provisional storage bounded by moment work, then apply physical cell
    // and node limits to the final certified component mapping, not its prefix.
    decltype(grouped_moments) final_moments;
    for (auto& [key, group] : grouped_moments) {
        charge(64);
        if (std::size_t(key.first[3]) >= component_mapping.size() || component_mapping[key.first[3]] >= 10000)
            throw std::invalid_argument("Invalid final section component mapping.");
        auto mapped     = key;
        mapped.first[3] = int(component_mapping[key.first[3]]);
        auto& final     = final_moments[mapped];
        if (!final.initialized) {
            final.anchor=final.minimum=group.minimum;
            final.maximum=group.maximum;
            final.initialized=true;
        }
        Eigen::Vector3d offset;
        for (unsigned axis=0;axis<3;++axis) {
            offset[axis]=double(group.anchor[axis]-final.anchor[axis])/
                double(axis==2 ? z_size(key.first[2]) : native_spacing[axis]);
            final.minimum[axis]=std::min(final.minimum[axis],group.minimum[axis]);
            final.maximum[axis]=std::max(final.maximum[axis],group.maximum[axis]);
        }
        const auto rebased=SectionIntegration::rebase_moments(group.integral,offset,Eigen::Vector3d::Ones());
        for (std::size_t i = 0; i < rebased.size(); ++i) final.integral[i] += rebased[i];
        final.bed_area += group.bed_area;
        for (auto& polygon : group.bed_trace) final.bed_trace.push_back(std::move(polygon));
        for (auto& trace : group.exact_bed_traces) final.exact_bed_traces.push_back(std::move(trace));
    }
    grouped_moments.swap(final_moments);
    final_moments.clear();
    integrating = false;
    // A detached component contained in one cell has no shared basis with
    // another cell. Put its eight Q1 coordinates on its occupied bounding box:
    // this is an invertible change of the same polynomial space, not geometry
    // simplification, stiffness regularization or an added physical support.
    struct LocalBasis { Eigen::Vector3d origin, size; };
    std::map<Index,LocalBasis> local_bases;
    std::map<int,Index> component_cell;
    std::set<int> ineligible;
    std::map<int,std::array<coord_t,6>> component_bounds;
    for (const auto& [key,group] : grouped_moments) {
        charge();
        const auto inserted=component_cell.emplace(key.first[3],key.first);
        if ((!inserted.second && inserted.first->second!=key.first) || group.bed_area>0.)
            ineligible.insert(key.first[3]);
        auto found=component_bounds.emplace(key.first[3],std::array<coord_t,6>{
            group.minimum[0],group.minimum[1],group.minimum[2],group.maximum[0],group.maximum[1],group.maximum[2]});
        for (unsigned axis=0;axis<3;++axis) {
            found.first->second[axis]=std::min(found.first->second[axis],group.minimum[axis]);
            found.first->second[axis+3]=std::max(found.first->second[axis+3],group.maximum[axis]);
        }
    }
    for (const auto& [component,key] : component_cell) if (!ineligible.count(component)) {
        charge();
        const auto& bounds=component_bounds.at(component);
        LocalBasis basis;
        for (unsigned axis=0;axis<3;++axis) {
            basis.origin[axis]=double(bounds[axis])*1e-9;
            basis.size[axis]=double(bounds[axis+3]-bounds[axis])*1e-9;
        }
        if ((basis.size.array()<=0.).any()) throw std::runtime_error("Empty isolated component support.");
        local_bases.emplace(key,basis);
    }
    std::map<int, Eigen::Matrix<double, 6, 1>> component_strain;
    std::map<int,SectionMaterial> component_material;
    using StrainTerms = std::array<Eigen::Matrix<double,6,1>,8>;
    std::map<std::pair<Index,std::size_t>,StrainTerms> spatial_loads;
    std::set<int> nonuniform_strain;
    double prescribed_stress_bound = 0.;
    for (auto& [key, group] : grouped_moments) {
        if (integration) charge(64);
        else if (checkpoint) checkpoint();
        if (!cells.count(key.first) && cells.size() >= limits.max_cells) {
            std::size_t required = 0;
            std::optional<Index> previous;
            for (const auto& [planned, unused] : grouped_moments) {
                charge();
                if (!previous || *previous != planned.first) ++required;
                previous = planned.first;
            }
            std::ostringstream message;
            message << "Section mechanics cell budget exceeded: required=" << required << " limit=" << limits.max_cells;
            throw ResourceLimit(message.str());
        }
        auto& cell           = cells[key.first];
        Eigen::Vector3d basis_origin(double(native_origin[0]+key.first[0]*native_spacing[0])*1e-9,
            double(native_origin[1]+key.first[1]*native_spacing[1])*1e-9,double(z_coordinate(key.first[2]))*1e-9);
        Eigen::Vector3d basis_size(double(native_spacing[0])*1e-9,double(native_spacing[1])*1e-9,double(z_size(key.first[2]))*1e-9);
        const Eigen::Vector3d original_origin=basis_origin, original_size=basis_size;
        const auto local_basis=local_bases.find(key.first);
        if (local_basis!=local_bases.end()) {
            basis_origin=local_basis->second.origin;
            basis_size=local_basis->second.size;
        }
        Eigen::Vector3d moment_offset;
        for (unsigned axis=0;axis<3;++axis)
            moment_offset[axis]=double(group.anchor[axis]-std::llround(basis_origin[axis]*1e9))*1e-9/basis_size[axis];
        group.integral=SectionIntegration::rebase_moments(group.integral,moment_offset,original_size.cwiseQuotient(basis_size));
        const auto polynomial = local_basis==local_bases.end() ? polynomial_at(key.first[2]) : strain_polynomials(basis_size);
        auto material = constitutive[key.second];
        std::optional<TrilinearThermalStrain> trilinear;
        if (load) {
            const Eigen::Vector3d minimum(double(native_origin[0]+key.first[0]*native_spacing[0])*1e-9,
                double(native_origin[1]+key.first[1]*native_spacing[1])*1e-9,double(z_coordinate(key.first[2]))*1e-9);
            const Eigen::Vector3d size(double(native_spacing[0])*1e-9,double(native_spacing[1])*1e-9,
                double(z_size(key.first[2]))*1e-9);
            const auto value = load({key.second,std::size_t(key.first[3]),
                {key.first[0],key.first[1],key.first[2]},minimum,minimum+size,group.integral[0]});
            if (const auto* affine = std::get_if<LinearThermalStrain>(&value)) {
                if (!affine->at_world_zero.allFinite() || !affine->gradient_per_m.allFinite())
                    throw std::invalid_argument("Invalid spatial section load.");
                material.eigenstrain = affine->at_world_zero;
                material.eigenstrain_gradient_per_m = affine->gradient_per_m;
            } else {
                trilinear = std::get<TrilinearThermalStrain>(value);
                for (const auto& corner : trilinear->at_corners)
                    if (!corner.allFinite()) throw std::invalid_argument("Invalid trilinear section load.");
                if (local_basis!=local_bases.end()) {
                    const auto original=*trilinear;
                    for (unsigned corner=0;corner<8;++corner) {
                        Eigen::Vector3d point=basis_origin;
                        for (unsigned axis=0;axis<3;++axis) if (corner&(1u<<axis)) point[axis]+=basis_size[axis];
                        const Eigen::Vector3d normalized=(point-original_origin).cwiseQuotient(original_size);
                        trilinear->at_corners[corner].setZero();
                        for (unsigned source=0;source<8;++source) {
                            double weight=1.;
                            for (unsigned axis=0;axis<3;++axis)
                                weight*=source&(1u<<axis) ? normalized[axis] : 1.-normalized[axis];
                            trilinear->at_corners[corner]+=weight*original.at_corners[source];
                        }
                    }
                }
                // General Q1 eigenstrain need not be compatible with a free
                // affine/quadratic displacement. Use the equilibrium solve.
                nonuniform_strain.insert(key.first[3]);
            }
        }
        const auto strain = component_strain.emplace(key.first[3], material.eigenstrain);
        const auto reference = component_material.emplace(key.first[3],material);
        if (!reference.second && (reference.first->second.eigenstrain_gradient_per_m.array() !=
                                   material.eigenstrain_gradient_per_m.array()).any())
            nonuniform_strain.insert(key.first[3]);
        if (!strain.second && (strain.first->second.array() != material.eigenstrain.array()).any())
            nonuniform_strain.insert(key.first[3]);
        const bool gradient = trilinear.has_value() || !material.eigenstrain_gradient_per_m.isZero(0.);
        if (gradient) charge(128); // Additional strain-moment products and stress bounds.
        StrainTerms strain_terms, stress_terms;
        for (auto& term : strain_terms) term.setZero();
        strain_terms[0] = material.eigenstrain;
        const unsigned terms = gradient ? 8 : 1;
        if (trilinear) {
            strain_terms = trilinear->at_corners;
            // Tensor-product finite differences turn corner values into
            // coefficients of 1,x,y,xy,z,xz,yz,xyz on the unit cell.
            for (unsigned axis = 0; axis < 3; ++axis)
                for (unsigned term = 0; term < 8; ++term)
                    if (term & (1u << axis)) strain_terms[term] -= strain_terms[term ^ (1u << axis)];
        } else if (gradient) {
            strain_terms[0] += material.eigenstrain_gradient_per_m*basis_origin;
            for (unsigned axis = 0; axis < 3; ++axis)
                strain_terms[1u << axis] = material.eigenstrain_gradient_per_m.col(axis)*basis_size[axis];
        }
        if (load) spatial_loads.emplace(key,strain_terms);
        for (unsigned t = 0; t < terms; ++t)
            stress_terms[t] = material.elasticity_pa*strain_terms[t];
        for (unsigned corner = 0; corner < (gradient ? 8u : 1u); ++corner) {
            Eigen::Matrix<double,6,1> stress = stress_terms[0];
            for (unsigned term = 1; term < terms; ++term)
                if ((corner & term) == term) stress += stress_terms[term];
            if (!stress.allFinite() || !std::isfinite(stress.norm()))
                throw std::invalid_argument("Nonfinite prescribed section stress.");
            prescribed_stress_bound = std::max(prescribed_stress_bound,stress.norm());
        }
        std::array<BMatrix, 8> db;
        for (unsigned i = 0; i < 8; ++i)
            db[i] = material.elasticity_pa * polynomial[i];
        cell.volume += group.integral[0];
        cell.bed_area += group.bed_area;
        for (auto& polygon : group.bed_trace) cell.bed_trace.push_back(std::move(polygon));
        for (auto& trace : group.exact_bed_traces) cell.exact_bed_traces.push_back(std::move(trace));
        // B and prescribed strain are multilinear. Their products
        // fit the same degree-two moments used for stiffness, even in cut cells.
        for (unsigned t = 0; t < terms; ++t)
            for (unsigned s = 0; s < terms; ++s)
                cell.initial_energy += .5*strain_terms[t].dot(stress_terms[s])*
                    group.integral[moment_index(t,s)];
        for (unsigned a = 0; a < 8; ++a) {
            for (unsigned t = 0; t < terms; ++t)
                cell.force.noalias() += group.integral[moment_index(a,t)] * polynomial[a].transpose() * stress_terms[t];
            for (unsigned b = 0; b < 8; ++b)
                cell.stiffness.noalias() += group.integral[moment_index(a, b)] * polynomial[a].transpose() * db[b];
        }
    }
    grouped_moments.clear();
    if (cells.empty())
        throw std::invalid_argument("Empty occupied section geometry.");
    SectionMechanicalResult result;
    result.uses_spatial_load = bool(load);
    result.prescribed_stress_bound_pa = prescribed_stress_bound;
    result.basis_origin_m = origin;
    result.basis_spacing_m = spacing;
    std::map<Index, std::size_t> nodes;
    std::set<std::size_t> fixed;
    std::set<std::size_t> footprint_nodes;
    std::set<int> geometrically_attached_components;
    std::array<std::set<int>, 3> occupied_axis_cells;
    for (auto& [key, cell] : cells) {
        charge();
        result.occupied_volume_m3 += cell.volume;
        result.prescribed_strain_energy_j += cell.initial_energy;
        result.bed_area_m2 += cell.bed_area;
        bool released = false;
        if (cell.bed_area > 0.) {
            geometrically_attached_components.insert(key[3]);
            for (const auto& box : release_cells) {
                charge();
                if (key[0] >= box[0] && key[0] < box[2] && key[1] >= box[1] && key[1] < box[3])
                    released = true;
            }
            if (released) result.released_bed_area_m2 += cell.bed_area;
        }
        for (unsigned axis = 0; axis < 3; ++axis)
            occupied_axis_cells[axis].insert(key[axis]);
        for (unsigned corner = 0; corner < 8; ++corner) {
            Index node = key;
            for (unsigned axis = 0; axis < 3; ++axis)
                node[axis] += ((corner >> axis) & 1);
            if (!nodes.count(node) && nodes.size() >= limits.max_nodes)
                throw ResourceLimit("Section mechanics node budget exceeded.");
            auto found         = nodes.emplace(node, nodes.size());
            cell.nodes[corner] = found.first->second;
            if (found.second) {
                auto position = (origin + Eigen::Vector3d(node[0] * spacing.x(), node[1] * spacing.y(), node[2] * spacing.z())).eval();
                if (!native_z.empty()) position.z() = double(z_coordinate(node[2])) * 1e-9;
                if (const auto local=local_bases.find(key);local!=local_bases.end()) {
                    position=local->second.origin;
                    for (unsigned axis=0;axis<3;++axis) if (corner&(1u<<axis)) position[axis]+=local->second.size[axis];
                }
                result.nodes_m.push_back(position);
            }
            if (cell.bed_area > 0. && !released && !(corner & 4))
                fixed.insert(found.first->second);
            if (cell.bed_area > 0. && !(corner & 4))
                footprint_nodes.insert(found.first->second);
        }
    }
    if (!(result.bed_area_m2 > 0.) && !limits.allow_free_components)
        throw std::invalid_argument("Section solid has no positive-area bed attachment.");
    std::map<int, std::vector<std::size_t>> component_nodes;
    std::set<int> anchored_components;
    for (const auto& [key, id] : nodes) {
        charge();
        component_nodes[key[3]].push_back(id);
        if (fixed.count(id)) anchored_components.insert(key[3]);
    }
    if (!released_regions.empty() && !(result.released_bed_area_m2 > 0.))
        throw std::invalid_argument("Section bed release regions contain no occupied bed area.");
    for (const auto component : geometrically_attached_components)
        if (!anchored_components.count(component))
            throw std::runtime_error("Section bed release fully detaches a component.");
    result.bed_connected_nodes.assign(nodes.size(), false);
    std::vector<Eigen::Vector3d> analytic_displacement(nodes.size(), Eigen::Vector3d::Zero());
    std::set<std::size_t> analytic_nodes;
    std::set<std::size_t> reference_dofs;
    std::map<int,Eigen::Vector3d> quadratic_centers;
    for (const auto& [component, ids] : component_nodes) {
        if (anchored_components.count(component)) {
            for (const auto id : ids) result.bed_connected_nodes[id] = true;
            continue;
        }
        if (!limits.allow_free_components)
            throw std::invalid_argument("Section solid contains an unattached component.");
        ++result.free_components;
        const Eigen::Vector3d center = result.nodes_m[ids.front()];
        const auto& material = component_material.at(component);
        if (!nonuniform_strain.count(component) && material.eigenstrain_gradient_per_m.isZero(0.)) {
            // Uniform eigenstrain admits u=epsilon*(x-origin), exactly in Q1,
            // regardless of shape or elasticity. No constrained mode or small
            // cut-cell stiffness solve is needed. Verify full equilibrium below.
            const auto& strain = component_strain.at(component);
            Eigen::Matrix3d tensor;
            tensor << strain[0], .5*strain[3], .5*strain[5],
                      .5*strain[3], strain[1], .5*strain[4],
                      .5*strain[5], .5*strain[4], strain[2];
            for (const auto id : ids) {
                charge();
                analytic_nodes.insert(id);
                analytic_displacement[id] = tensor * (result.nodes_m[id]-center);
            }
            ++result.affine_free_components;
            continue;
        }
        if (!nonuniform_strain.count(component)) {
            quadratic_centers.emplace(component,center);
            for (const auto id : ids) {
                charge(27);
                analytic_nodes.insert(id);
                analytic_displacement[id] = compatible_displacement(material,center,result.nodes_m[id]);
            }
            ++result.quadratic_free_components;
            continue;
        }
        double length = 0.;
        for (const auto id : ids) {
            charge();
            length = std::max(length, (result.nodes_m[id] - center).norm());
        }
        if (!std::isfinite(length) || length <= 0.)
            throw std::runtime_error("Free section reference frame is unresolved.");
        Eigen::MatrixXd rigid(6, 3 * ids.size());
        for (std::size_t i = 0; i < ids.size(); ++i) {
            charge(6);
            const Eigen::Vector3d r = (result.nodes_m[ids[i]] - center) / length;
            rigid.col(3*i)     << 1., 0., 0., 0., r.z(), -r.y();
            rigid.col(3*i + 1) << 0., 1., 0., -r.z(), 0., r.x();
            rigid.col(3*i + 2) << 0., 0., 1., r.y(), -r.x(), 0.;
        }
        // Select six independent reference coordinates, not physical supports.
        // C*R invertible means ANY displacement can satisfy C*u=0 by adding
        // only a rigid motion. Hence this removes no strain mode and adds no
        // elastic stiffness. Full equilibrium at these coordinates is checked
        // below; their reactions are never reported as bed forces.
        Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(rigid);
        qr.setThreshold(1e-12);
        if (qr.rank() != 6)
            throw std::runtime_error("Free section rigid-body modes are unresolved.");
        for (unsigned i = 0; i < 6; ++i) {
            const auto column = std::size_t(qr.colsPermutation().indices()[i]);
            reference_dofs.insert(3 * ids[column / 3] + column % 3);
        }
    }
    result.reference_dofs = reference_dofs.size();
    for (unsigned axis = 0; axis < 3; ++axis)
        result.occupied_axis_cells[axis] = occupied_axis_cells[axis].size();
    std::set<std::size_t> contacting;
    std::vector<double> contact_scale(nodes.size(),0.);
    const auto solve_contact_set = [&] {
        result.elastic_energy_j  = 0.;
        result.relative_residual = 0.;
        result.componentwise_relative_residual = 0.;
        std::fill(contact_scale.begin(), contact_scale.end(), 0.);
        std::vector<int> free(nodes.size() * 3, -1);
        int free_count = 0;
        for (std::size_t n = 0; n < nodes.size(); ++n)
            if (!fixed.count(n) && !analytic_nodes.count(n))
                for (unsigned a = 0; a < 3; ++a)
                    if (!reference_dofs.count(3 * n + a) && !(a == 2 && contacting.count(n)))
                        free[3 * n + a] = free_count++;
        if (std::size_t(free_count) > limits.max_free_dofs)
            throw ResourceLimit("Section mechanics free-DOF budget exceeded.");
        result.solved_dofs  = std::size_t(free_count);
        Eigen::VectorXd rhs = Eigen::VectorXd::Zero(free_count);
        std::vector<Eigen::Triplet<double>> triplets;
        for (const auto& [key, cell] : cells) {
            charge();
            if (!cell.stiffness.allFinite() || !cell.force.allFinite())
                throw std::invalid_argument("Nonfinite section stiffness.");
            for (unsigned a = 0; a < 24; ++a) {
                const int row = free[3 * cell.nodes[a / 3] + a % 3];
                if (row < 0)
                    continue;
                rhs[row] += cell.force[a];
                for (unsigned b = 0; b < 24; ++b) {
                    const int col = free[3 * cell.nodes[b / 3] + b % 3];
                    if (col < 0)
                        continue;
                    if (triplets.size() >= limits.max_triplets)
                        throw ResourceLimit("Section mechanics assembly budget exceeded.");
                    triplets.emplace_back(row, col, cell.stiffness(a, b));
                }
            }
        }
        Eigen::SparseMatrix<double> matrix(free_count, free_count);
        matrix.setFromTriplets(triplets.begin(), triplets.end());
        result.assembly_triplets = triplets.size();
        result.matrix_entries    = std::size_t(matrix.nonZeros());
        // Sparse assembly owns its compressed entries now. Do not retain the much
        // larger duplicate cell contributions while allocating the factor.
        std::vector<Eigen::Triplet<double>>().swap(triplets);
        Eigen::VectorXd displacement = Eigen::VectorXd::Zero(free_count);
        if (free_count) {
            charge();
            // Congruence scaling changes coordinates, not the elastic problem:
            // B = S K S, g = S f, u = S y. It prevents a small, soft but anchored
            // component from being rejected solely against another component's
            // stiffness scale. Near-dependent cut bases must still pass the same
            // positive-pivot test; no stiffness floor or regularization is added.
            Eigen::VectorXd scale = matrix.diagonal();
            for (Eigen::Index i = 0; i < scale.size(); ++i) {
                charge();
                if (!std::isfinite(scale[i]) || scale[i] <= 0.)
                    throw std::runtime_error("Section mechanical diagonal is nonpositive or unresolved.");
                scale[i] = 1. / std::sqrt(scale[i]);
                if (!std::isfinite(scale[i]))
                    throw std::runtime_error("Section mechanical scaling is nonfinite.");
            }
            Eigen::SparseMatrix<double> equilibrated = matrix;
            for (int column = 0; column < equilibrated.outerSize(); ++column) {
                charge();
                for (Eigen::SparseMatrix<double>::InnerIterator entry(equilibrated, column); entry; ++entry) {
                    entry.valueRef() = (entry.value() * scale[entry.row()]) * scale[entry.col()];
                    if (!std::isfinite(entry.value()))
                        throw std::runtime_error("Section equilibrated stiffness is nonfinite.");
                }
            }
            const Eigen::VectorXd scaled_rhs = scale.cwiseProduct(rhs);
            if (!scaled_rhs.allFinite())
                throw std::runtime_error("Section equilibrated load is nonfinite.");
            // Select the ordering explicitly so the allocation preflight and the
            // numerical factor see the identical graph. NaturalOrdering below
            // prevents the factorizer from choosing a different fill pattern.
            Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> inverse_order;
            Eigen::AMDOrdering<int> ordering;
            Eigen::SparseMatrix<double> symmetric = equilibrated.selfadjointView<Eigen::Lower>();
            ordering(symmetric, inverse_order);
            const Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> permutation = inverse_order.inverse();
            Eigen::SparseMatrix<double> ordered(free_count, free_count);
            ordered.selfadjointView<Eigen::Upper>() = symmetric.selfadjointView<Eigen::Lower>().twistedBy(permutation);
            symmetric.resize(0, 0);
            std::optional<std::size_t> expected_entries;
            if (limits.sparse_factor)
                expected_entries = count_ldlt_factor_entries(ordered, limits.sparse_factor->max_entries, limits.sparse_factor->max_work,
                                                             checkpoint);
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>, Eigen::Upper, Eigen::NaturalOrdering<int>> solver;
            solver.compute(ordered);
            if (solver.info() != Eigen::Success || !solver.vectorD().allFinite())
                throw std::runtime_error("Section mechanical factorization failed or is nonfinite.");
            result.factor_entries = std::size_t(solver.matrixL().nestedExpression().nonZeros());
            if (expected_entries && result.factor_entries != *expected_entries)
                throw std::runtime_error("Section factor pattern differs from its allocation preflight.");
            if (solver.vectorD().minCoeff() <= 1e-12 * solver.vectorD().cwiseAbs().maxCoeff()) {
                Eigen::Index pivot = 0;
                solver.vectorD().minCoeff(&pivot);
                Eigen::VectorXd selector = Eigen::VectorXd::Zero(free_count);
                selector[pivot]          = 1.;
                selector                 = permutation.transpose() * selector;
                Eigen::Index row         = 0;
                selector.maxCoeff(&row);
                bool pivot_bed_connected = false;
                std::size_t pivot_node = result.nodes_m.size();
                for (std::size_t dof = 0; dof < free.size(); ++dof)
                    if (free[dof] == row) {
                        pivot_bed_connected = result.bed_connected_nodes[dof / 3];
                        pivot_node = dof / 3;
                        break;
                    }
                std::ostringstream message;
                message << "Section mechanical system is singular or numerically unresolved: scaled_min_pivot="
                        << solver.vectorD().minCoeff() << " scaled_max_pivot=" << solver.vectorD().cwiseAbs().maxCoeff()
                        << " nonpositive_pivots=" << (solver.vectorD().array() <= 0.).count()
                        << " free_components=" << result.free_components << " reference_dofs=" << result.reference_dofs
                        << " free_dofs=" << free_count << " cells=" << cells.size()
                        << " minimum_pivot_bed_connected=" << pivot_bed_connected;
                if (pivot_node < result.nodes_m.size()) {
                    message << " minimum_pivot_node_m=" << result.nodes_m[pivot_node].transpose();
                    for (const auto& [component, ids] : component_nodes) {
                        if (std::find(ids.begin(),ids.end(),pivot_node)==ids.end()) continue;
                        Eigen::Vector3d low=result.nodes_m[ids.front()], high=low;
                        for (const auto id : ids) {
                            charge();
                            low=low.cwiseMin(result.nodes_m[id]);
                            high=high.cwiseMax(result.nodes_m[id]);
                        }
                        message << " minimum_pivot_component=" << component << " component_nodes=" << ids.size()
                                << " component_basis_extent_m=" << (high-low).transpose();
                        break;
                    }
                }
                throw std::runtime_error(message.str());
            }
            const Eigen::VectorXd ordered_rhs         = permutation * scaled_rhs;
            Eigen::VectorXd scaled_displacement = permutation.transpose() * solver.solve(ordered_rhs);
            displacement                              = scale.cwiseProduct(scaled_displacement);
            charge();
            if (solver.info() != Eigen::Success || !displacement.allFinite())
                throw std::runtime_error("Section mechanical solve failed.");
            Eigen::VectorXd residual;
            // Bounded iterative refinement reuses the existing factor to solve
            // K*delta = K*u-f. It changes neither stiffness nor tolerances, and
            // targets row-wise backward error so large components cannot hide
            // inaccurate equilibrium in a weak cut cell.
            for (unsigned refinement = 0; ; ++refinement) {
                residual = matrix*displacement-rhs;
                Eigen::VectorXd denominator = rhs.cwiseAbs();
                for (int column = 0; column < matrix.outerSize(); ++column) {
                    charge();
                    for (Eigen::SparseMatrix<double>::InnerIterator entry(matrix,column); entry; ++entry)
                        denominator[entry.row()] += std::abs(entry.value())*std::abs(displacement[entry.col()]);
                }
                double error = 0.;
                for (Eigen::Index i = 0; i < denominator.size(); ++i) {
                    charge();
                    if (!std::isfinite(denominator[i]) || !std::isfinite(residual[i]))
                        throw std::runtime_error("Nonfinite section componentwise equilibrium.");
                    error = std::max(error,denominator[i] == 0. ? (residual[i] == 0. ? 0. : 1.) :
                        std::abs(residual[i])/denominator[i]);
                }
                result.componentwise_relative_residual = error;
                if (error <= 1e-8) break;
                if (refinement == 3) {
                    std::ostringstream message;
                    message << "Section componentwise equilibrium is unresolved after refinement: backward_error=" << error;
                    throw std::runtime_error(message.str());
                }
                const Eigen::VectorXd correction_rhs = permutation*scale.cwiseProduct(residual);
                const Eigen::VectorXd correction = permutation.transpose()*solver.solve(correction_rhs);
                if (solver.info() != Eigen::Success || !correction.allFinite())
                    throw std::runtime_error("Section residual correction failed.");
                scaled_displacement -= correction;
                displacement = scale.cwiseProduct(scaled_displacement);
                ++result.linear_refinement_steps;
            }
            const Eigen::VectorXd scaled_residual = equilibrated * scaled_displacement - scaled_rhs;
            if (!residual.allFinite() || !scaled_residual.allFinite())
                throw std::runtime_error("Section mechanical residual is nonfinite.");
            const double rhs_norm    = rhs.stableNorm();
            result.relative_residual = residual.stableNorm() / std::max(rhs_norm, 1e-30);
            if (!std::isfinite(rhs_norm) || !std::isfinite(result.relative_residual) || result.relative_residual > 1e-8)
                throw std::runtime_error("Section mechanical equilibrium is unresolved.");
            const double scaled_norm = scaled_rhs.stableNorm();
            if (!std::isfinite(scaled_norm) || (scaled_norm == 0. && scaled_residual.stableNorm() != 0.) ||
                (scaled_norm > 0. && scaled_residual.stableNorm() / scaled_norm > 1e-8))
                throw std::runtime_error("Section equilibrated equilibrium is unresolved.");
        }
        result.displacement_m = analytic_displacement;
        result.generalized_bed_reaction_n.assign(nodes.size(), Eigen::Vector3d::Zero());
        result.generalized_contact_reaction_n.assign(nodes.size(), 0.);
        for (std::size_t n = 0; n < nodes.size(); ++n)
            for (unsigned axis = 0; axis < 3; ++axis)
                if (free[3 * n + axis] >= 0)
                    result.displacement_m[n][axis] = displacement[free[3 * n + axis]];
        double energy_calculation_scale = 0.;
        std::vector<double> reference_force(nodes.size() * 3, 0.), reference_scale(nodes.size() * 3, 0.);
        for (const auto& [key, cell] : cells) {
            charge();
            // These isolated bodies use the exact compatible quadratic field,
            // not Q1 interpolation of its nodal samples. Their elastic strain
            // is identically zero; they transmit no load to another component.
            if (quadratic_centers.count(key[3])) continue;
            Vector24 u;
            for (unsigned n = 0; n < 8; ++n)
                u.segment<3>(3 * n) = result.displacement_m[cell.nodes[n]];
            const Vector24 force   = cell.stiffness * u - cell.force;
            const double quadratic = .5 * u.dot(cell.stiffness * u), load = u.dot(cell.force);
            result.elastic_energy_j += quadratic - load + cell.initial_energy;
            energy_calculation_scale += std::abs(quadratic) + std::abs(load) + std::abs(cell.initial_energy);
            for (unsigned n = 0; n < 8; ++n)
                if (fixed.count(cell.nodes[n]))
                    result.generalized_bed_reaction_n[cell.nodes[n]] += force.segment<3>(3 * n);
            for (unsigned n = 0; n < 8; ++n)
                if (contacting.count(cell.nodes[n])) {
                    const unsigned a = 3 * n + 2;
                    result.generalized_contact_reaction_n[cell.nodes[n]] += force[a];
                    contact_scale[cell.nodes[n]] += std::abs(cell.force[a]);
                    for (unsigned b = 0; b < 24; ++b)
                        contact_scale[cell.nodes[n]] += std::abs(cell.stiffness(a, b)) * std::abs(u[b]);
                }
            for (unsigned a = 0; a < 24; ++a) {
                const auto dof = 3 * cell.nodes[a / 3] + a % 3;
                if (!reference_dofs.count(dof) && !analytic_nodes.count(cell.nodes[a / 3]))
                    continue;
                charge(24);
                reference_force[dof] += force[a];
                reference_scale[dof] += std::abs(cell.force[a]);
                for (unsigned b = 0; b < 24; ++b)
                    reference_scale[dof] += std::abs(cell.stiffness(a, b)) * std::abs(u[b]);
            }
        }
        for (std::size_t dof = 0; dof < reference_force.size(); ++dof) {
            if (!reference_dofs.count(dof) && !analytic_nodes.count(dof / 3))
                continue;
            const double force = reference_force[dof], scale = reference_scale[dof];
            if (!std::isfinite(force) || !std::isfinite(scale) || (scale == 0. ? force != 0. : std::abs(force) > 1e-8 * scale))
                throw std::runtime_error("Free section reference or affine solution would carry a physical reaction.");
        }
        if (!std::isfinite(result.elastic_energy_j) || !std::isfinite(energy_calculation_scale) ||
            result.elastic_energy_j < -128. * std::numeric_limits<double>::epsilon() * energy_calculation_scale)
            throw std::runtime_error("Invalid section elastic energy.");
        // Cancellation of the three energy terms can leave a roundoff-negative
        // zero. The tolerance scales in joules with this calculation, not a fixed
        // absolute energy that could hide failure on a very small or soft body.
        result.elastic_energy_j = std::max(0., result.elastic_energy_j);
    };
    const auto inspect_trace = [&] {
        result.bed_lift_upper_bound_m = 0.;
        result.maximum_footprint_lift_m = 0.;
        result.maximum_footprint_penetration_m = 0.;
        bool penetrates = false;
        if (released_regions.empty()) return penetrates;
        for (const auto& [key,cell] : cells) {
            if (!(cell.bed_area > 0.)) continue;
            std::array<double,4> values;
            double displacement_scale = 0.;
            for (unsigned n = 0; n < 4; ++n) {
                const auto& u = result.displacement_m[cell.nodes[n]];
                values[n] = u.z();
                displacement_scale = std::max(displacement_scale,u.cwiseAbs().maxCoeff());
                result.bed_lift_upper_bound_m = std::max(result.bed_lift_upper_bound_m,u.z());
            }
            const std::array<coord_t,2> trace_origin{native_origin[0]+key[0]*native_spacing[0],native_origin[1]+key[1]*native_spacing[1]};
            const std::array<coord_t,2> trace_size{native_spacing[0],native_spacing[1]};
            const auto consider=[&](double minimum,double maximum) {
                penetrates = penetrates || minimum < -(1e-12+1e-8*displacement_scale);
                result.maximum_footprint_lift_m = std::max(result.maximum_footprint_lift_m,maximum);
                result.maximum_footprint_penetration_m = std::max(result.maximum_footprint_penetration_m,-minimum);
            };
            if (!cell.bed_trace.empty()) {
                const auto extrema=section_trace_extrema(cell.bed_trace,trace_origin,trace_size,values,[&] { charge(); });
                consider(extrema.minimum,extrema.maximum);
            }
            for (const auto& trace:cell.exact_bed_traces) {
                // Polygon extraction and edge visitation are temporary; the
                // retained handles have their separate aggregate bed budget.
                if (trace.vertex_count()>live_points/2)
                    throw ResourceLimit("Exact bed trace scratch point budget exceeded.");
                const auto consumed=integration ? assembly_work : work;
                if (consumed==limits.max_work) throw ResourceLimit("Exact bed trace work budget exhausted.");
                const auto extrema=exact_section_bilinear_extrema(trace,trace_origin,trace_size,values,
                    std::min<std::size_t>(4000000,limits.max_work-consumed),[&] { charge(0); });
                charge(extrema.work);
                // Use outward bounds for contact detection; rational corners
                // and edge-interior extrema never pass through native rounding.
                consider(extrema.minimum_bounds[0],extrema.maximum_bounds[1]);
            }
            if (cell.bed_trace.empty() && cell.exact_bed_traces.empty())
                throw std::invalid_argument("Missing occupied bed trace.");
        }
        return penetrates;
    };
    std::set<std::vector<std::size_t>> contact_history;
    for (unsigned iteration = 0;; ++iteration) {
        if (checkpoint) checkpoint();
        if (!contact_history.emplace(contacting.begin(),contacting.end()).second)
            throw std::runtime_error("Section cell-contact active set cycles.");
        solve_contact_set();
        const bool penetrates = inspect_trace();
        if (!cell_contact) {
            if (penetrates)
                throw std::runtime_error("Section bed release requires compressive re-contact on its occupied footprint.");
            break;
        }
        result.contact_iterations = iteration+1;
        result.conservative_cell_contact = !contacting.empty();
        // An unconstrained minimum that is feasible on the actual footprint
        // already solves the physical unilateral problem. Do not add contact
        // solely because its polynomial extension enters the bed in a void.
        if (contacting.empty() && !penetrates) break;
        auto next_contacting = contacting;
        for (const auto n : footprint_nodes) {
            charge();
            if (fixed.count(n)) continue;
            const double gap = result.displacement_m[n].z();
            if (!std::isfinite(result.generalized_contact_reaction_n[n]) || !std::isfinite(contact_scale[n]))
                throw std::runtime_error("Nonfinite section cell-contact reaction.");
            const double gap_tolerance = 1e-12+1e-8*result.displacement_m[n].cwiseAbs().maxCoeff();
            if (contacting.count(n)) {
                if (result.generalized_contact_reaction_n[n] < -1e-8*contact_scale[n])
                    next_contacting.erase(n);
            } else if (gap < -gap_tolerance) next_contacting.insert(n);
        }
        if (next_contacting == contacting) {
            if (penetrates) throw std::runtime_error("Section cell-contact footprint still penetrates the bed.");
            break;
        }
        if (iteration+1 >= cell_contact->max_iterations)
            throw ResourceLimit("Section cell-contact iteration budget exceeded.");
        contacting = std::move(next_contacting);
    }
    // Locate samples in certified occupied slabs, never in the background
    // grid alone. Half-open Z intervals give one material-side value at a
    // horizontal interface; XY boundary points deliberately remain unsampled.
    for (std::size_t request_index = 0; request_index < samples.size(); ++request_index) {
        charge();
        const auto& request = samples[request_index];
        const Point xy(static_cast<coord_t>(std::llround(request.position_m.x() * 1e9)),
                       static_cast<coord_t>(std::llround(request.position_m.y() * 1e9)));
        const coord_t z = std::llround(request.position_m.z() * 1e9);
        std::optional<std::size_t> component;
        for (const auto provisional : sample_components[request_index]) {
            charge();
            if (provisional >= component_mapping.size())
                throw std::invalid_argument("Invalid section sample source.");
            const auto mapped = component_mapping[provisional];
            if (component && *component != mapped)
                throw std::runtime_error("Ambiguous section sample component.");
            component = mapped;
        }
        if (!component) {
            ++result.unsampled_points;
            continue;
        }
        const Index key{index(xy.x(), 0), index(xy.y(), 1), index(z, 2), int(*component)};
        const auto cell = cells.find(key);
        if (cell == cells.end()) {
            ++result.unsampled_points;
            continue;
        }
        const auto quadratic = quadratic_centers.find(int(*component));
        if (quadratic != quadratic_centers.end()) {
            const auto displacement = compatible_displacement(component_material.at(int(*component)),quadratic->second,request.position_m);
            result.samples.push_back({request.id,request.position_m,displacement,
                Eigen::Matrix<double,6,1>::Zero(),0.,0.,false});
            continue;
        }
        Eigen::Vector3d local;
        for (unsigned axis = 0; axis < 3; ++axis)
            local[axis] = (request.position_m[axis] * 1e9 - double(native_origin[axis] + key[axis] * native_spacing[axis])) /
                          double(native_spacing[axis]);
        if (!native_z.empty())
            local.z() = (request.position_m.z() * 1e9 - double(z_coordinate(key[2]))) / double(z_size(key[2]));
        const auto local_basis=local_bases.find(key);
        if (local_basis!=local_bases.end())
            local=(request.position_m-local_basis->second.origin).cwiseQuotient(local_basis->second.size);
        const auto polynomial = local_basis==local_bases.end() ? polynomial_at(key[2]) : strain_polynomials(local_basis->second.size);
        Vector24 u;
        Eigen::Vector3d interpolated = Eigen::Vector3d::Zero();
        for (unsigned corner = 0; corner < 8; ++corner) {
            u.segment<3>(3 * corner) = result.displacement_m[cell->second.nodes[corner]];
            double weight            = 1.;
            for (unsigned axis = 0; axis < 3; ++axis)
                weight *= corner & (1u << axis) ? local[axis] : 1. - local[axis];
            interpolated += weight * u.segment<3>(3 * corner);
        }
        BMatrix b = BMatrix::Zero();
        for (unsigned monomial = 0; monomial < 8; ++monomial) {
            double weight = 1.;
            for (unsigned axis = 0; axis < 3; ++axis)
                if (monomial & (1u << axis))
                    weight *= local[axis];
            b += weight * polynomial[monomial];
        }
        auto material = constitutive[request.material];
        Eigen::Matrix<double,6,1> eigenstrain = material.eigenstrain + material.eigenstrain_gradient_per_m*request.position_m;
        if (load) {
            const auto& terms = spatial_loads.at({key,request.material});
            eigenstrain.setZero();
            for (unsigned term = 0; term < 8; ++term) {
                double weight = 1.;
                for (unsigned axis = 0; axis < 3; ++axis)
                    if (term & (1u << axis)) weight *= local[axis];
                eigenstrain += weight*terms[term];
            }
        }
        const Eigen::Matrix<double, 6, 1> stress = material.elasticity_pa *
            (b * u - eigenstrain);
        Eigen::Matrix3d tensor;
        tensor << stress[0], stress[3], stress[5], stress[3], stress[1], stress[4], stress[5], stress[4], stress[2];
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eig(tensor, Eigen::EigenvaluesOnly);
        const double mises = std::sqrt(
            .5 * (std::pow(stress[0] - stress[1], 2) + std::pow(stress[1] - stress[2], 2) + std::pow(stress[2] - stress[0], 2)) +
            3. * stress.tail<3>().squaredNorm());
        if (!stress.allFinite() || !interpolated.allFinite() || eig.info() != Eigen::Success || !std::isfinite(mises))
            throw std::runtime_error("Invalid recovered section field.");
        result.samples.push_back({request.id, request.position_m, interpolated, stress, eig.eigenvalues().maxCoeff(), mises,
                                  anchored_components.count(int(*component)) != 0});
    }
    return result;
}
} // namespace

SectionBedReactionProjection project_section_bed_reactions(const SectionMechanicalResult& field,
    const Eigen::Vector2d& origin, const Eigen::Vector2d& spacing, const std::function<void()>& checkpoint)
{
    const auto native = [](double value) -> coord_t {
        if (!std::isfinite(value) || std::abs(value) > 10. ||
            std::abs(value * 1e9 - std::round(value * 1e9)) > 1e-5)
            throw std::invalid_argument("Reaction reporting grid must use native coordinates.");
        return std::llround(value * 1e9);
    };
    if (field.nodes_m.empty() || field.nodes_m.size() > 96000 ||
        field.generalized_bed_reaction_n.size() != field.nodes_m.size() ||
        field.bed_connected_nodes.size() != field.nodes_m.size())
        throw std::invalid_argument("Invalid section reaction field.");
    std::array<coord_t,2> start, step, source_start, source_step;
    for (unsigned axis = 0; axis < 2; ++axis) {
        start[axis] = native(origin[axis]);
        step[axis] = native(spacing[axis]);
        source_start[axis] = native(field.basis_origin_m[axis]);
        source_step[axis] = native(field.basis_spacing_m[axis]);
        if (step[axis] <= 0 || source_step[axis] <= 0 || step[axis] % source_step[axis] ||
            (start[axis]-source_start[axis]) % source_step[axis])
            throw std::invalid_argument("Reaction reporting grid must be nested over the section basis.");
    }
    const coord_t bed = native(field.basis_origin_m.z());
    std::map<std::array<coord_t,2>,Eigen::Vector3d> reactions;
    for (std::size_t n = 0; n < field.nodes_m.size(); ++n) {
        if (checkpoint) checkpoint();
        const auto& force = field.generalized_bed_reaction_n[n];
        if (!force.allFinite()) throw std::invalid_argument("Nonfinite section reaction.");
        if (force.isZero(0.)) continue;
        if (!field.bed_connected_nodes[n] || native(field.nodes_m[n].z()) != bed)
            throw std::invalid_argument("Section reaction is not on the bonded bed.");
        std::array<coord_t,2> index;
        Eigen::Vector2d fraction;
        for (unsigned axis = 0; axis < 2; ++axis) {
            const coord_t position = native(field.nodes_m[n][axis]);
            if ((position-source_start[axis]) % source_step[axis])
                throw std::invalid_argument("Section reaction node is outside its declared basis.");
            const coord_t delta = position-start[axis];
            index[axis] = delta/step[axis] - (delta%step[axis] < 0 ? 1 : 0);
            fraction[axis] = double(delta-index[axis]*step[axis])/double(step[axis]);
        }
        for (unsigned corner = 0; corner < 4; ++corner) {
            auto key = index;
            double weight = 1.;
            for (unsigned axis = 0; axis < 2; ++axis) {
                const bool high = (corner & (1u << axis)) != 0;
                key[axis] += high;
                weight *= high ? fraction[axis] : 1.-fraction[axis];
            }
            if (weight == 0.) continue;
            auto entry = reactions.try_emplace(key,Eigen::Vector3d::Zero());
            entry.first->second += weight * force;
            if (!entry.first->second.allFinite()) throw std::runtime_error("Reaction projection overflow.");
        }
    }
    SectionBedReactionProjection result;
    for (const auto& [key,force] : reactions) {
        result.nodes_m.emplace_back(double(start[0]+key[0]*step[0])*1e-9,
                                    double(start[1]+key[1]*step[1])*1e-9, double(bed)*1e-9);
        result.generalized_reaction_n.push_back(force);
        result.tensile_generalized_force_n += std::max(0.,-force.z());
    }
    if (!std::isfinite(result.tensile_generalized_force_n)) throw std::runtime_error("Reaction projection overflow.");
    return result;
}

SectionMechanicalResult solve_deposited_sections(const SectionTopology& topology,
                                                 const std::vector<SectionMaterial>& materials,
                                                 const Eigen::Vector3d& origin,
                                                 const Eigen::Vector3d& spacing,
                                                 const SectionSolveLimits& limits,
                                                 const std::function<void()>& checkpoint,
                                                 const std::vector<SectionSampleRequest>& samples,
                                                 const std::vector<SectionBedReleaseRegion>& released_regions,
                                                 const std::optional<SectionCellContact>& cell_contact,
                                                 const SectionLoadField& load)
{
    if (topology.pieces().empty() || topology.pieces().size() > 10000 || topology.components() > 10000)
        throw std::invalid_argument("Invalid section topology.");
    return solve_section_source(
        topology.bed_z_m(),
        [&](const auto& accept) {
            for (const auto& piece : topology.pieces())
                accept(piece);
            std::vector<std::size_t> mapping(topology.components());
            for (std::size_t i = 0; i < mapping.size(); ++i)
                mapping[i] = i;
            return mapping;
        },
        materials, origin, spacing, limits, checkpoint, samples, nullptr, nullptr, {}, released_regions, cell_contact, load);
}

SectionMechanicalResult solve_local_section_tiles(const std::vector<LocalSectionTile>& tiles,
                                                  double bed,
                                                  const std::vector<SectionMaterial>& materials,
                                                  const Eigen::Vector3d& origin,
                                                  const Eigen::Vector3d& spacing,
                                                  const SectionSolveLimits& limits,
                                                  const SectionTopologyLimits& topology_limits,
                                                  const std::function<void()>& checkpoint,
                                                  const std::vector<SectionSampleRequest>& samples,
                                                  const std::vector<double>& z_planes,
                                                  const std::vector<SectionBedReleaseRegion>& released_regions,
                                                  const std::optional<SectionCellContact>& cell_contact,
                                                  const SectionLoadField& load)
{
    if (tiles.empty() || tiles.size() > 4096)
        throw std::invalid_argument("Invalid mechanical tile source.");
    std::size_t expected_free_components = 0;
    auto result = solve_section_source(
        bed,
        [&](const auto& accept) {
            LocalTopologyBuilder builder(bed, topology_limits, accept, checkpoint, {}, nullptr, LocalTileOrder::Arbitrary,
                limits.allow_free_components ? LocalAttachmentPolicy::AllowFreeComponents : LocalAttachmentPolicy::RequireBed);
            const auto source = builder.source();
            for (const auto& tile : tiles)
                builder.add_tile(tile);
            // No finish or field publication if any producer/callback failed.
            const auto certificate = builder.finish();
            if (certificate.source() != source || certificate.components() > 10000)
                throw std::invalid_argument("Mechanical topology source identity mismatch.");
            expected_free_components = certificate.free_components();
            return certificate.component_mapping();
        },
        materials, origin, spacing, limits, checkpoint, samples, nullptr, nullptr, z_planes, released_regions, cell_contact, load);
    if (result.free_components != expected_free_components)
        throw std::runtime_error("Section bed attachment disagrees with its topology certificate.");
    return result;
}
LocalSectionMechanicalResult solve_local_deposition_sections(const std::vector<DepositedSegment>& segments,
    const std::vector<SectionMaterial>& materials, const LocalSectionStreamSettings& settings,
    const std::function<std::optional<SectionOwnership>(double)>& ownership,
    const std::function<void()>& checkpoint, const std::vector<DepositedSegment>* partition_envelopes,
    const std::vector<SectionSampleRequest>& samples, StreamIntegrationStatistics* statistics,
    StreamTopologyStatistics* topology_statistics, LocalSectionStatistics* geometry_statistics)
{
    if (topology_statistics) *topology_statistics = {};
    if (geometry_statistics) *geometry_statistics = {};
    LocalSectionMechanicalResult result;
    std::size_t expected_free_components = 0;
    result.field = solve_section_source(settings.bed_z_m, [&](const auto& accept) {
        LocalTopologyBuilder builder(settings.bed_z_m, settings.topology, accept, checkpoint,
            settings.streaming_topology, topology_statistics, LocalTileOrder::Lexicographic,
            settings.mechanics.allow_free_components ? LocalAttachmentPolicy::AllowFreeComponents : LocalAttachmentPolicy::RequireBed);
        const auto source = builder.source();
        LocalSectionPacketVisitor visitor;
        visitor.begin_tile = [&](const LocalSectionTile& header) { builder.begin_tile(header); };
        visitor.band = [&](const LocalSectionTile& packet) { builder.add_band(packet); };
        visitor.end_tile = [&] { builder.end_tile(); };
        result.geometry = visit_local_deposition_packets(segments, settings.origin_m.head<2>(), settings.tile_size_m,
            settings.bed_z_m, settings.maximum_height_m, visitor,
            settings.geometry, ownership, checkpoint, settings.interval_subdivisions, geometry_statistics, partition_envelopes);
        // The traversal returned successfully; only now may connectivity be
        // finalized. A callback exception never reaches finish or publication.
        const auto certificate = builder.finish();
        if (certificate.source() != source || certificate.components() > 10000)
            throw std::invalid_argument("Mechanical topology source identity mismatch.");
        expected_free_components = certificate.free_components();
        return certificate.component_mapping();
    }, materials, settings.origin_m, settings.basis_spacing_m, settings.mechanics, checkpoint, samples,
        settings.integration ? &*settings.integration : nullptr, statistics, settings.basis_z_planes_m,
        settings.released_bed_regions, settings.cell_contact, settings.load);
    if (result.field.free_components != expected_free_components)
        throw std::runtime_error("Section bed attachment disagrees with its topology certificate.");
    return result;
}

ReferenceSectionMechanicalResult solve_reference_deposition_sections(
    const std::vector<DepositedSegment>& reference,const std::vector<DepositedSegment>& deposited,
    const std::vector<SectionMaterial>& materials,const LocalSectionStreamSettings& settings,
    const SectionGrowthStreamLimits& activation_limits,
    const std::function<std::optional<SectionOwnership>(double)>& ownership,const std::function<void()>& checkpoint)
{
    SectionGrowthStreamSettings activation;
    static_cast<SectionGrowthStreamLimits&>(activation)=activation_limits;
    activation.origin_m=settings.origin_m.head<2>(); activation.tile_size_m=settings.tile_size_m;
    activation.bed_z_m=settings.bed_z_m; activation.maximum_height_m=settings.maximum_height_m;
    activation.interval_subdivisions=settings.interval_subdivisions;
    activation.geometry=settings.geometry; activation.topology=settings.topology;
    activation.connectivity=settings.streaming_topology.value_or(StreamTopologyLimits{});
    ReferenceSectionMechanicalResult result;
    std::size_t expected_free=0;
    result.field=solve_section_source(settings.bed_z_m,[&](const auto& accept) {
        const auto proof=visit_reference_activation(reference,deposited,activation,accept,ownership,checkpoint);
        if (proof.active.components()>10000 || (!settings.mechanics.allow_free_components && proof.active.free_components()))
            throw std::invalid_argument("Reference activation has unsupported mechanical components: total=" +
                std::to_string(proof.active.components()) + ", free=" + std::to_string(proof.active.free_components()) + ".");
        expected_free=proof.active.free_components();
        result.excluded_nominal_volume_m3=proof.excluded_nominal_volume_m3;
        return proof.active.component_mapping();
    },materials,settings.origin_m,settings.basis_spacing_m,settings.mechanics,checkpoint,{},
        settings.integration ? &*settings.integration : nullptr,nullptr,settings.basis_z_planes_m,
        settings.released_bed_regions,settings.cell_contact,settings.load);
    if (result.field.free_components!=expected_free)
        throw std::runtime_error("Reference mechanical attachment disagrees with topology.");
    return result;
}

namespace {
SectionDirectionalProbe probe_refinement(const std::vector<SectionMaterial>& materials,
                                                 const Eigen::Vector3d& spacing,
                                                 const std::function<void(unsigned)>& progress,
                                                 const std::function<void()>& checkpoint,
                                                 const std::function<SectionMechanicalResult(const Eigen::Vector3d&,
                                                     const std::function<void()>&)>& solve)
{
    SectionDirectionalProbe probe;
    double energy_floor = 0., force_floor = 0.;
    bool callback_failed = false;
    const auto check     = [&] {
        try {
            if (checkpoint)
                checkpoint();
        } catch (...) {
            callback_failed = true;
            throw;
        }
    };
    try {
        for (unsigned run = 0; run < 4; ++run) {
            check();
            try {
                if (progress)
                    progress(run);
            } catch (...) {
                callback_failed = true;
                throw;
            }
            Eigen::Vector3d refined = spacing;
            if (run)
                refined[run - 1] *= .5;
            // Keep just one mechanical field alive; only scalar diagnostics
            // survive this iteration. Each solve keeps the existing hard caps.
            const auto field = solve(refined, check);
            if (run == 0) {
                if (!field.uses_spatial_load) for (const auto &material : materials) {
                    const auto stress = (material.elasticity_pa * material.eigenstrain).eval();
                    energy_floor = std::max(energy_floor, .5 * material.eigenstrain.dot(stress) * field.occupied_volume_m3);
                    force_floor = std::max(force_floor, stress.norm() * field.bed_area_m2);
                }
                energy_floor = std::max(energy_floor,field.prescribed_strain_energy_j);
                force_floor = std::max(force_floor,field.prescribed_stress_bound_pa*field.bed_area_m2);
                energy_floor *= 128. * std::numeric_limits<double>::epsilon();
                force_floor *= 128. * std::numeric_limits<double>::epsilon();
            }
            double tensile = 0.;
            for (const auto &reaction : field.generalized_bed_reaction_n)
                tensile += std::max(0., -reaction.z());
            if (!std::isfinite(tensile)) throw std::runtime_error("Nonfinite directional bed reaction.");
            SectionBasisSample sample{refined, field.elastic_energy_j, tensile, field.relative_residual,
                                      field.nodes_m.size(), field.occupied_axis_cells};
            if (run == 0) probe.reference = sample;
            else probe.directions[run - 1] = sample;
        }
    } catch (const ResourceLimit&) {
        if (callback_failed) throw;
        probe.state = SectionDirectionalProbe::State::ResourceLimited;
        return probe;
    }
    const auto relative_change = [](double a, double b, double floor) {
        const double scale = std::max(std::abs(a), std::abs(b));
        return scale > floor ? std::max(0., std::abs(a - b) - floor) / scale : 0.;
    };
    unsigned selected = 0;
    double largest_change = -1.;
    for (unsigned axis = 0; axis < 3; ++axis) {
        const auto &sample = *probe.directions[axis];
        probe.energy_change[axis] = relative_change(probe.reference->elastic_energy_j, sample.elastic_energy_j, energy_floor);
        probe.force_change[axis] = relative_change(probe.reference->tensile_generalized_bed_force_n,
                                                   sample.tensile_generalized_bed_force_n, force_floor);
        const double change = std::max(probe.energy_change[axis], probe.force_change[axis]);
        if (change > largest_change) { largest_change = change; selected = axis; }
    }
    if (largest_change > 0.) probe.recommended_spacing_m = probe.directions[selected]->spacing_m;
    probe.state = SectionDirectionalProbe::State::Complete;
    return probe;
}
} // namespace

SectionDirectionalProbe probe_section_refinement(const SectionTopology& topology,
    const std::vector<SectionMaterial>& materials, const Eigen::Vector3d& origin,
    const Eigen::Vector3d& spacing, const SectionSolveLimits& limits,
    const std::function<void(unsigned)>& progress, const std::function<void()>& checkpoint)
{
    return probe_refinement(materials, spacing, progress, checkpoint,
        [&](const Eigen::Vector3d& refined, const std::function<void()>& check) {
            return solve_deposited_sections(topology, materials, origin, refined, limits, check);
        });
}

SectionDirectionalProbe probe_local_deposition_refinement(const std::vector<DepositedSegment>& segments,
    const std::vector<SectionMaterial>& materials, const LocalSectionStreamSettings& settings,
    const std::function<std::optional<SectionOwnership>(double)>& ownership,
    const std::function<void(unsigned)>& progress, const std::function<void()>& checkpoint,
    const std::vector<DepositedSegment>* partition_envelopes)
{
    if (!settings.basis_z_planes_m.empty())
        throw std::invalid_argument("A graded section grid requires explicit plane refinement.");
    std::exception_ptr callback_error;
    const auto guard = [&](auto callback) {
        return [&, callback = std::move(callback)](const auto& argument) {
            try { return callback(argument); }
            catch (...) { callback_error = std::current_exception(); throw; }
        };
    };
    std::function<std::optional<SectionOwnership>(double)> owner;
    const auto guarded_source = guard(ownership);
    if (ownership) owner = [&](double z) {
        auto result = guarded_source(z);
        if (result) {
            if (result->resolve_conflict) result->resolve_conflict = guard(std::move(result->resolve_conflict));
            if (result->resolve_sampled_conflict) result->resolve_sampled_conflict = guard(std::move(result->resolve_sampled_conflict));
            if (result->resolve_partition) result->resolve_partition = guard(std::move(result->resolve_partition));
        }
        return result;
    };
    auto probe = probe_refinement(materials, settings.basis_spacing_m, progress, checkpoint,
        [&](const Eigen::Vector3d& refined, const std::function<void()>& check) {
            auto local = settings;
            local.basis_spacing_m = refined;
            if (settings.load) local.load = guard(settings.load);
            return solve_local_deposition_sections(segments, materials, local, owner, check, partition_envelopes).field;
        });
    if (callback_error) std::rethrow_exception(callback_error);
    return probe;
}

SectionBasisStudy refine_deposited_sections(const SectionTopology& topology,
                                            const std::vector<SectionMaterial>& materials,
                                            const Eigen::Vector3d& origin,
                                            const Eigen::Vector3d& initial_spacing,
                                            unsigned max_levels,
                                            double tolerance,
                                            const SectionSolveLimits& limits,
                                            const std::function<void(unsigned)>& progress,
                                            const std::function<void()>& checkpoint,
                                            const std::vector<SectionSampleRequest>& samples)
{
    if (max_levels == 0 || max_levels > 5 || !std::isfinite(tolerance) || tolerance <= 0. || tolerance > 1. || materials.empty() ||
        materials.size() > 256)
        throw std::invalid_argument("Invalid section basis refinement settings.");
    std::vector<double> energy_density, eigenstress;
    for (const auto& material : materials) {
        if (!material.eigenstrain.allFinite() || !material.eigenstrain_gradient_per_m.allFinite())
            throw std::invalid_argument("Invalid refinement eigenstrain.");
        const auto stress = (validated_elasticity(material.elasticity_pa) * material.eigenstrain).eval();
        energy_density.push_back(.5 * material.eigenstrain.dot(stress));
        eigenstress.push_back(stress.norm());
    }
    double energy_scale = 0., bed_area = 0., maximum_eigenstress = 0.;
    std::size_t exact_scale_work=0,exact_scale_points=0;
    for (const auto& piece : topology.pieces()) {
        if (checkpoint)
            checkpoint();
        const auto& slab = piece.slab;
        if (slab.material >= materials.size())
            throw std::invalid_argument("Missing refinement material.");
        double area = 0.;
        if (slab.footprint.is_exact()) {
            const auto count=slab.footprint.vertex_count();
            if (count>limits.max_points/2 || count>limits.max_points-exact_scale_points)
                throw ResourceLimit("Exact refinement scale point budget exceeded.");
            exact_scale_points+=count;
            if (const auto bounds=slab.footprint.bounds()) {
                if (exact_scale_work>=limits.max_work)
                    throw ResourceLimit("Exact refinement scale work budget exhausted.");
                const auto measure=exact_section_horizontal_area(slab.footprint.exact(),
                    {bounds->min.x(),bounds->min.y(),0},
                    {bounds->max.x()-bounds->min.x(),bounds->max.y()-bounds->min.y(),1},0,
                    std::min<std::size_t>(4000000,limits.max_work-exact_scale_work),checkpoint);
                exact_scale_work+=measure.work;
                area=measure.value;
                if (!(area>0.) || !std::isfinite(area))
                    throw std::invalid_argument("Invalid exact refinement scale area.");
            }
        } else for (const auto& polygon : slab.footprint.native())
            area += polygon.area() * 1e-18;
        energy_scale += energy_density[slab.material] * area * (slab.top_m - slab.bottom_m);
        maximum_eigenstress = std::max(maximum_eigenstress, eigenstress[slab.material]);
        if (std::llround(slab.bottom_m * 1e9) == std::llround(topology.bed_z_m() * 1e9))
            bed_area += area;
    }
    double force_scale = maximum_eigenstress * bed_area;
    if (!std::isfinite(energy_scale) || !std::isfinite(force_scale))
        throw std::invalid_argument("Nonfinite section refinement scales.");
    const auto stable = [tolerance](double a, double b, double scale) {
        return std::abs(a - b) <= tolerance * std::max(std::abs(a), std::abs(b)) + 128. * std::numeric_limits<double>::epsilon() * scale;
    };
    SectionBasisStudy study;
    unsigned stable_comparisons   = 0;
    Eigen::Vector3d spacing       = initial_spacing;
    bool callback_failed          = false;
    const auto guarded_checkpoint = [&] {
        try {
            if (checkpoint)
                checkpoint();
        } catch (...) {
            callback_failed = true;
            throw;
        }
    };
    try {
        for (unsigned level = 0; level < max_levels; ++level) {
            guarded_checkpoint();
            try {
                if (progress)
                    progress(level);
            } catch (...) {
                callback_failed = true;
                throw;
            }
            auto result          = solve_deposited_sections(topology, materials, origin, spacing, limits, guarded_checkpoint, samples);
            energy_scale = std::max(energy_scale,result.prescribed_strain_energy_j);
            force_scale = std::max(force_scale,result.prescribed_stress_bound_pa*result.bed_area_m2);
            double tensile_force = 0.;
            for (const auto& reaction : result.generalized_bed_reaction_n)
                tensile_force += std::max(0., -reaction.z());
            if (!std::isfinite(tensile_force))
                throw std::runtime_error("Nonfinite generalized tensile reaction.");
            if (!study.history.empty()) {
                const auto& previous = study.history.back();
                bool refined_all_axes = true;
                for (unsigned axis = 0; axis < 3; ++axis)
                    refined_all_axes = refined_all_axes && result.occupied_axis_cells[axis] > previous.occupied_axis_cells[axis];
                stable_comparisons   = refined_all_axes && result.nodes_m.size() > previous.nodes &&
                                               stable(previous.elastic_energy_j, result.elastic_energy_j, energy_scale) &&
                                               stable(previous.tensile_generalized_bed_force_n, tensile_force, force_scale) ?
                                           stable_comparisons + 1 :
                                           0;
            }
            study.history.push_back({spacing, result.elastic_energy_j, tensile_force, result.relative_residual,
                                     result.nodes_m.size(), result.occupied_axis_cells});
            if (stable_comparisons >= 2) {
                study.state = SectionBasisStudy::State::BasisConverged;
                study.field = std::move(result);
                break;
            }
            spacing *= .5;
        }
    } catch (const ResourceLimit&) {
        if (callback_failed)
            throw;
        study.state = SectionBasisStudy::State::ResourceLimited;
    }
    return study;
}
} // namespace Slic3r::Thermoelastic
