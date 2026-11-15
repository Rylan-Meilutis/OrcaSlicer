#pragma once
#include "ThermoelasticLayers.hpp"
#include "ThermoelasticMesh.hpp"
#include "ThermoelasticSections.hpp"
#include "ThermoelasticSectionTopology.hpp"
#include <Eigen/Core>
#include <cstddef>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace Slic3r::Thermoelastic {
struct SectionStudySettings
{
    double max_section_height_m  = 0.;
    unsigned max_geometry_levels = 3, max_basis_levels = 3;
    // Optional diagnostic-led seed selection before the unchanged full-axis
    // refinement study. Each probe uses at most four bounded solves.
    unsigned max_initial_basis_probes = 0;
    double relative_tolerance = 0.;
    // Explicit opt-in to bounded tile-local quadrature. This retains the
    // nominal paths/voids but avoids globally multiplying every scarf endpoint
    // height across unrelated regions. Numerical refinement is still required.
    std::optional<Eigen::Vector2d> local_tile_size_m;
    // Explicit quadrature-partition approximation, local path only. Endpoint
    // tops round up and bottoms down to choose slab boundaries on this grid,
    // independently halved at each geometry level. Samples evaluate ORIGINAL
    // beads: dimensions/materials/paths are never changed. Not an auto default.
    std::optional<double> initial_z_quantization_m;
    SectionSweepLimits sections;
    SectionTopologyLimits topology;
    SectionSolveLimits mechanics;
    std::size_t max_total_sections = 10000, max_total_segment_visits = 20000000;
    std::size_t max_total_boundary_points = 2000000, max_total_work = 100000000;
};
struct SectionGeometrySample
{
    double maximum_section_height_m, volume_m3, elastic_energy_j, tensile_generalized_bed_force_n;
    std::size_t sections, components;
    SectionBasisStudy::State basis_state;
    double z_quantization_m = 0., maximum_partition_shift_m = 0.;
    // Partition endpoint coalescing/enclosure only; no bead is removed or given
    // a finite cap. Original zero-height tips remain zero during evaluation.
    std::size_t collapsed_endpoint_intervals = 0;
    std::size_t expanded_zero_endpoints = 0;
    double maximum_boundary_uncertainty_m = 0.;
    // Upper bounds on changed total occupancy integrated across sampled slabs;
    // not a material-transfer volume or a mechanical error estimate.
    double added_representation_volume_bound_m3 = 0., removed_representation_volume_bound_m3 = 0.;
};
struct DepositionSectionStudy
{
    enum class State { Incomplete, NominalGeometryAndBasisConverged, ResourceLimited,
        ApproximateGeometryAndBasisConverged };
    State state = State::Incomplete;
    std::string failure_reason; // Bounded internal diagnostic, not a safety verdict.
    std::vector<SectionGeometrySample> history;
    // Scalar diagnostics from the last attempted geometry, including failed
    // basis convergence. Never retain its unqualified displacement field.
    std::vector<SectionBasisSample> last_basis_history;
    std::vector<SectionDirectionalProbe> last_initial_basis_probes;
    std::optional<SectionMechanicalResult> field;
    bool piecewise_constant_geometry = false;
    std::size_t section_count = 0, segment_visits = 0, boundary_points = 0, reserved_work = 0;
    std::size_t last_geometry_tiles = 0, last_geometry_slabs = 0;
};
// Sample the explicit nominal bead geometry, retaining every native-Z bead
// boundary. Nothing is inferred from filament names, display widths or CAD.
// Geometry and basis each require two stable comparisons, except horizontal
// rectangular envelopes without ownership are exactly constant between their
// included boundaries. This qualifies the specified nominal representation
// only, not physical deposition, thermal history, bonding or bed adhesion.
// Work is conservatively reserved from primitive per-call maxima. Ownership
// reconciliation reports consumed work, so unused per-plane reservations are
// released; other primitives retain their maximum charge. Exhaustion returns no field;
// caller callback exceptions (including ResourceLimit) always propagate.
// Progress arguments are (geometry-level index, basis-solve index), NOT
// completed/total work; neither is a progress percentage.
DepositionSectionStudy study_deposition_sections(const std::vector<DepositedSegment>& segments,
                                                 const std::vector<SectionMaterial>& materials,
                                                 const Eigen::Vector3d& origin_m,
                                                 const Eigen::Vector3d& initial_basis_spacing_m,
                                                 double bed_z_m,
                                                 const SectionStudySettings& settings,
                                                 const std::function<std::optional<SectionOwnership>(double)>& ownership = {},
                                                 const std::function<void(unsigned, unsigned)>& progress                 = {},
                                                 const std::function<void()>& checkpoint                                 = {},
                                                 const std::vector<SectionSampleRequest>& samples = {});
} // namespace Slic3r::Thermoelastic
