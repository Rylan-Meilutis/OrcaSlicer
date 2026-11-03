#pragma once

#include "ThermoelasticLayers.hpp"
#include "ThermoelasticMaterial.hpp"
#include "ThermoelasticSectionTopology.hpp"
#include "ThermoelasticThermalStrain.hpp"
#include <Eigen/Core>
#include <array>
#include <cstddef>
#include <functional>
#include <optional>
#include <variant>
#include <vector>

namespace Slic3r::Thermoelastic {
struct LocalSectionTile;
struct SectionMaterial
{
    ElasticityMatrix elasticity_pa;
    Eigen::Matrix<double, 6, 1> eigenstrain = Eigen::Matrix<double, 6, 1>::Zero();
    // Prescribed linear strain field in world coordinates (metres). The
    // constant term above is its value at world zero. This is a loading
    // scenario, not an inferred temperature or a heat-transfer solution.
    Eigen::Matrix<double, 6, 3> eigenstrain_gradient_per_m = Eigen::Matrix<double, 6, 3>::Zero();
};
struct SectionLoadQuery {
    std::size_t material, component;
    std::array<int, 3> cell;
    Eigen::Vector3d minimum_m, maximum_m;
    double occupied_volume_m3;
};
// Spatial loading on certified occupied cell/material groups. A supplied field
// replaces the material's prescribed strain, never its stiffness or geometry.
// Affine strain uses world coordinates; trilinear corner strain uses the
// requested cell box. Both must be valid throughout that box.
// Missing thermal coverage must throw; there is no zero-strain fallback.
// Called once per occupied group after topology completion; exceptions propagate.
using SectionSpatialStrain = std::variant<LinearThermalStrain,TrilinearThermalStrain>;
using SectionLoadField = std::function<SectionSpatialStrain(const SectionLoadQuery&)>;
struct SectionFactorLimits {
    // Strict factor off-diagonal entries, each a double plus a storage index.
    // The optional extended-grid path counts them BEFORE factor allocation.
    std::size_t max_entries = 8000000, max_work = 32000000;
};
struct SectionSolveLimits
{
    // Explicit offline ceilings: 4000 cells, 6000 nodes, 12000 solved DOFs.
    // At the DOF ceiling a dense triangular LDLT factor alone would use about
    // 0.81 GiB (double values + 32-bit indices); input/assembly are additional.
    // Defaults below are unchanged and deliberately smaller.
    std::size_t max_cells = 1000, max_nodes = 6000, max_free_dofs = 6000;
    std::size_t max_points = 200000, max_work = 2000000, max_triplets = 1000000;
    std::size_t max_bed_trace_points = 200000;
    // Explicit analysis of free bodies, not invented bed/interface attachment.
    // Their displacement has an arbitrary rigid-body reference frame.
    bool allow_free_components = false;
    // Explicit offline mode only: up to 16x the grid/assembly ceilings above,
    // with a symbolic factor-storage preflight (at most 64M entries / ~733 MiB).
    // Input, matrices, geometry and temporary storage have separate bounds.
    // No default budget increases: every larger allocation must be requested.
    std::optional<SectionFactorLimits> sparse_factor;
};
struct SectionSampleRequest {
    std::size_t id;
    Eigen::Vector3d position_m;
    std::size_t material;
};
struct SectionBedReleaseRegion {
    // Prescribed traction-free footprint rectangle, metres. Bounds must align
    // with XY basis cell faces, and remain fixed during refinement.
    Eigen::Vector2d minimum_m, maximum_m;
};
struct SectionCellContact {
    // Explicit conservative contact subspace: normal displacement is >= 0 at
    // every occupied bed-cell corner, including corners outside cut material.
    // This bounds penetration but can overconstrain a coarse cut-cell trace.
    // It requires basis refinement, not interpretation as local contact stress.
    unsigned max_iterations = 16;
};
struct SectionSpatialSample {
    std::size_t id;
    Eigen::Vector3d position_m, displacement_m;
    Eigen::Matrix<double, 6, 1> stress_pa;
    double maximum_principal_pa, von_mises_pa;
    bool bed_connected = true;
};
struct SectionMechanicalResult
{
    Eigen::Vector3d basis_origin_m = Eigen::Vector3d::Zero(), basis_spacing_m = Eigen::Vector3d::Zero();
    std::vector<Eigen::Vector3d> nodes_m, displacement_m, generalized_bed_reaction_n;
    // Detached one-cell components use support-scaled Q1 coordinates; their
    // nodes need not lie on the background lattice. Bed nodes remain on it.
    std::array<std::size_t, 3> occupied_axis_cells{};
    double occupied_volume_m3 = 0., bed_area_m2 = 0., elastic_energy_j = 0.;
    // Loading scales for numerical roundoff, not the relaxed field's stress.
    double prescribed_strain_energy_j = 0., prescribed_stress_bound_pa = 0.;
    // Maximum Q1 bed-cell corner lift bounds the cut footprint's displacement;
    // a background corner may lie outside deposited material.
    double released_bed_area_m2 = 0., bed_lift_upper_bound_m = 0.;
    double maximum_footprint_lift_m = 0., maximum_footprint_penetration_m = 0.;
    std::vector<double> generalized_contact_reaction_n;
    unsigned contact_iterations = 0;
    bool conservative_cell_contact = false;
    bool uses_spatial_load = false;
    double relative_residual = 0.; // Coarse Galerkin equilibrium only.
    double componentwise_relative_residual = 0.;
    unsigned linear_refinement_steps = 0; // Cumulative across contact iterations.
    std::size_t free_components = 0, reference_dofs = 0;
    std::size_t affine_free_components = 0;
    // Exact compatible linear-eigenstrain solution, not Q1 interpolation.
    std::size_t quadratic_free_components = 0;
    std::size_t solved_dofs = 0, factor_entries = 0, assembly_triplets = 0, matrix_entries = 0;
    std::vector<bool> bed_connected_nodes;
    // Point samples of this Q1 subspace (or the exact compatible free-body
    // solution), not independently converged local stresses. Points outside
    // certified material (including holes) are omitted.
    std::vector<SectionSpatialSample> samples;
    std::size_t unsampled_points = 0;
};
struct SectionBedReactionProjection {
    std::vector<Eigen::Vector3d> nodes_m, generalized_reaction_n;
    double tensile_generalized_force_n = 0.;
};
// I^T r on an explicit fixed XY Q1 reporting grid. Its planes must align with
// the source basis, with integer-multiple spacings. Preserves virtual work for
// reporting-grid displacement fields, total force and moment. This measures
// restraint at a FIXED spatial scale, not local traction or adhesive capacity;
// tensile/compressive contributions within a reporting basis can cancel.
// It does not replace raw-force convergence checks or qualify a field.
SectionBedReactionProjection project_section_bed_reactions(const SectionMechanicalResult& field,
    const Eigen::Vector2d& origin_m, const Eigen::Vector2d& spacing_m,
    const std::function<void()>& checkpoint = {});
// Static Q1 kinematic subspace on a regular background grid. Integrates the
// actual clipped slab polygons (including holes) and material tensors, never
// fills background cells or scales elasticity by a density fraction. Grid Z=0
// is the fully bonded rigid bed. Exact polynomial moments are used apart from
// native polygon coordinate quantization. Slabs must not overlap in volume.
// Positive-area bed traces are constrained exactly. Reactions are generalized
// basis forces, NOT local contact traction or available adhesion. No artificial
// physical supports, stiffness floor or disconnected-component stabilization
// is introduced. Explicit free-body analysis selects six reference coordinates
// per free component and verifies zero reaction there. Its absolute rigid-body
// motion is undetermined; bed_connected labels distinguish these results.
// Result remains UNQUALIFIED until independently refined in basis and geometry;
// topology factory verifies the sampled slab solid only, not original deposition.
// Separate physical components receive distinct basis nodes, even when coincident.
SectionMechanicalResult solve_deposited_sections(const SectionTopology& topology,
                                                 const std::vector<SectionMaterial>& materials,
                                                 const Eigen::Vector3d& origin_m,
                                                 const Eigen::Vector3d& spacing_m,
                                                 const SectionSolveLimits& limits        = {},
                                                 const std::function<void()>& checkpoint = {},
                                                 const std::vector<SectionSampleRequest>& samples = {},
                                                 const std::vector<SectionBedReleaseRegion>& released_regions = {},
                                                 const std::optional<SectionCellContact>& cell_contact = {},
                                                 const SectionLoadField& load = {});

// Closed tile producer -> connectivity -> moments -> solve path. Normalized
// pieces are integrated while borrowed; no historical integration polygons
// are retained. Only complete traversal and a matching final connectivity
// certificate permit assembly. Same unqualified static model and cumulative
// limits as solve_deposited_sections; no public raw moment input is accepted.
// Optional Z planes use the graded-grid contract of LocalSectionStreamSettings;
// empty preserves regular spacing, and incomplete solid coverage is rejected.
SectionMechanicalResult solve_local_section_tiles(const std::vector<LocalSectionTile>& tiles,
    double bed_z_m, const std::vector<SectionMaterial>& materials,
    const Eigen::Vector3d& origin_m, const Eigen::Vector3d& spacing_m,
    const SectionSolveLimits& limits = {}, const SectionTopologyLimits& topology_limits = {},
    const std::function<void()>& checkpoint = {}, const std::vector<SectionSampleRequest>& samples = {},
    const std::vector<double>& basis_z_planes_m = {},
    const std::vector<SectionBedReleaseRegion>& released_regions = {},
    const std::optional<SectionCellContact>& cell_contact = {},
    const SectionLoadField& load = {});

// Release regions remove only bed restraint, not deposited geometry. Shared
// nodes on retained cell traces stay fixed. Newly unrestrained components and
// bed penetration reject the default traction-free solve. Optional cell contact
// supplies a conservative normal-displacement subspace, not local plate pressure.
// Bed penetration/lift are checked on the occupied polygonal Q1 traces,
// including edge-interior extrema; a separate cell-corner bound is retained.
// Results remain prescribed-release studies,
// not an adhesion prediction, initiation criterion or qualified peel evolution.

struct SectionBasisSample
{
    Eigen::Vector3d spacing_m;
    double elastic_energy_j, tensile_generalized_bed_force_n, relative_residual;
    std::size_t nodes;
    std::array<std::size_t, 3> occupied_axis_cells{};
};
struct SectionBasisStudy
{
    enum class State { Incomplete, BasisConverged, ResourceLimited };
    State state = State::Incomplete;
    std::vector<SectionBasisSample> history;
    // No coarse field is returned after incomplete/budget-limited work.
    std::optional<SectionMechanicalResult> field;
};
struct SectionDirectionalProbe
{
    enum class State { Incomplete, Complete, ResourceLimited };
    State state = State::Incomplete;
    std::optional<SectionBasisSample> reference;
    std::array<std::optional<SectionBasisSample>, 3> directions;
    // Relative changes after halving one axis, excluding eigenstrain-scaled
    // roundoff. These are diagnostics, not a discretization error bound.
    std::array<double, 3> energy_change{}, force_change{};
    std::optional<Eigen::Vector3d> recommended_spacing_m;
};
// Four bounded solves (reference and independent X/Y/Z halvings). Recommends
// the axis with greatest energy/reaction change, but returns NO displacement
// field or convergence qualification. All three probes must succeed; a budget
// failure cannot select an axis using only the cheaper completed directions.
// No recommendation is made when changes are within roundoff in all axes.
// At most four per-solve work budgets are consumed. Cancellation propagates.
SectionDirectionalProbe probe_section_refinement(const SectionTopology& topology,
    const std::vector<SectionMaterial>& materials, const Eigen::Vector3d& origin_m,
    const Eigen::Vector3d& spacing_m, const SectionSolveLimits& limits = {},
    const std::function<void(unsigned)>& progress = {}, const std::function<void()>& checkpoint = {});
// Nested background-grid halvings on IDENTICAL certified slab geometry. Two
// successive stable energy/generalized tensile-force comparisons are required.
// Each comparison must increase occupied cell resolution in every axis; extra
// lateral nodes alone cannot qualify an unresolved through-thickness basis.
// This qualifies only this kinematic basis, not geometry, thermal history,
// material/interface data, local bed traction or physical adhesion. Roundoff
// floors use input eigenstrain energy and eigenstress times actual bed area.
SectionBasisStudy refine_deposited_sections(const SectionTopology& topology,
                                            const std::vector<SectionMaterial>& materials,
                                            const Eigen::Vector3d& origin_m,
                                            const Eigen::Vector3d& initial_spacing_m,
                                            unsigned max_levels,
                                            double relative_tolerance,
                                            const SectionSolveLimits& limits              = {},
                                            const std::function<void(unsigned)>& progress = {},
                                            const std::function<void()>& checkpoint       = {},
                                            const std::vector<SectionSampleRequest>& samples = {});
} // namespace Slic3r::Thermoelastic
