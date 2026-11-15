#include "ThermomechanicalSimulation.hpp"
#include "ThermalDeposition.hpp"
#include "ThermalHistory.hpp"
#include "ThermoelasticMaterial.hpp"
#include "ThermoelasticReduction.hpp"
#include "ThermoelasticMesh.hpp"
#include "ThermoelasticPeeling.hpp"
#include "ThermoelasticSolver.hpp"
#include "ThermoelasticThermalStrain.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {
SimulationResult simulate_deposition_grid(const DepositionGrid &grid,
    const std::vector<CellActivation> &activation, const std::vector<ThermalMechanicalMaterial> &materials,
    const std::vector<BoundaryScenario> &scenarios, double resistance, double end,
    const SimulationLimits &limits, const std::function<void()> &checkpoint,
    const std::vector<BedReleaseEvent> &bed_releases, const BedContactLimits &contact_limits)
{
    const auto valid_positive = [](double x) { return std::isfinite(x) && x > 0.; };
    if (limits.reduction) {
        const auto& settings = *limits.reduction;
        if (settings.axis > 2 || settings.strides.empty() || settings.strides.size() > 5 ||
            !valid_positive(settings.displacement_relative_tolerance) ||
            !valid_positive(settings.energy_relative_tolerance) ||
            settings.displacement_relative_tolerance > 1. || settings.energy_relative_tolerance > 1.)
            throw std::invalid_argument("Invalid reduction convergence settings.");
        unsigned previous = 0;
        for (unsigned stride : settings.strides) {
            if (stride == 0 || stride > 16 || (previous && (stride >= previous || previous % stride != 0)))
                throw std::invalid_argument("Reduction strides must be nested and decreasing.");
            previous = stride;
        }
    }
    if (activation.size() != grid.voxels.size() || grid.voxels.empty() || scenarios.empty() ||
        !grid.origin_m.allFinite() || !grid.spacing_m.allFinite() || (grid.spacing_m.array() <= 0.).any() ||
        !std::isfinite(end) || end < 0. || !std::isfinite(resistance) || resistance < 0.)
        throw std::invalid_argument("Invalid deposition-grid simulation inputs.");
    if (grid.voxels.size() > 200000 || materials.size()>256 || limits.max_schedule_entries==0 ||
        limits.max_schedule_entries > 10000000 || limits.max_mechanical_nodes==0 || limits.max_mechanical_nodes>12000 ||
        grid.voxels.size() > limits.thermal.solve.max_cells || scenarios.size() > limits.thermal.max_events)
        throw ResourceLimit("Deposition-grid simulation exceeds its cell/event budget.");
    if (bed_releases.size() > grid.voxels.size() || bed_releases.size() > limits.thermal.max_events)
        throw ResourceLimit("Bed release schedule exceeds its event budget.");
    if (!bed_releases.empty() && limits.reduction)
        throw std::invalid_argument("Reduced mechanics does not support bed release histories.");
    const double volume = grid.spacing_m.prod();
    if (!valid_positive(volume)) throw std::invalid_argument("Invalid thermal cell volume.");
    for (const auto &material : materials) {
        if (!valid_positive(material.conductivity_w_m_k) || !valid_positive(material.volumetric_heat_capacity_j_m3_k))
            throw std::invalid_argument("Invalid thermal material properties.");
        thermal_eigenstrain(material.expansion, material.local_to_world,
                            material.stress_free_temperature_k, material.stress_free_temperature_k);
        if (material.elasticity_pa) validated_elasticity(*material.elasticity_pa);
    }
    std::map<std::array<int,3>,std::size_t> occupied;
    std::vector<std::size_t> order(grid.voxels.size());
    std::iota(order.begin(),order.end(),std::size_t(0));
    for (std::size_t i=0;i<grid.voxels.size();++i) {
        if (checkpoint && i%128==0) checkpoint();
        const auto &voxel=grid.voxels[i];
        if (voxel.material >= materials.size() || voxel.index[2]<0 ||
            !std::isfinite(activation[i].time_s) || activation[i].time_s<0. || !valid_positive(activation[i].temperature_k) ||
            !occupied.emplace(voxel.index,i).second)
            throw std::invalid_argument("Invalid or duplicated deposition cell.");
        for (const int coordinate : voxel.index)
            if (coordinate==std::numeric_limits<int>::min() || coordinate==std::numeric_limits<int>::max())
                throw std::invalid_argument("Deposition index cannot form neighbors.");
    }
    std::stable_sort(order.begin(),order.end(),[&](std::size_t a,std::size_t b){return activation[a].time_s<activation[b].time_s;});
    std::vector<bool> scheduled_release(grid.voxels.size(), false);
    double previous_release = -1.;
    for (const auto &release : bed_releases) {
        if (checkpoint) checkpoint();
        if (!std::isfinite(release.time_s) || release.time_s < 0. || release.time_s < previous_release ||
            release.voxel >= grid.voxels.size() || scheduled_release[release.voxel] ||
            grid.voxels[release.voxel].index[2] != 0 || release.time_s < activation[release.voxel].time_s ||
            !std::isfinite(release.coefficient_w_m2_k) || release.coefficient_w_m2_k < 0. ||
            !valid_positive(release.temperature_k))
            throw std::invalid_argument("Invalid prescribed bed release event.");
        scheduled_release[release.voxel] = true;
        previous_release = release.time_s;
    }
    double previous=-1.;
    for (const auto &scenario : scenarios) {
        if (!std::isfinite(scenario.time_s) || scenario.time_s<0. || scenario.time_s<=previous ||
            !std::isfinite(scenario.bed_coefficient_w_m2_k) || scenario.bed_coefficient_w_m2_k<0. || !valid_positive(scenario.bed_temperature_k))
            throw std::invalid_argument("Invalid thermal boundary scenario.");
        previous=scenario.time_s;
        if (scenario.exposed.size()>16) throw ResourceLimit("Too many thermal surface exchanges.");
        for (const auto &surface : scenario.exposed)
            if (!std::isfinite(surface.coefficient_w_m2_k) || surface.coefficient_w_m2_k<0. ||
                !valid_positive(surface.temperature_k) || surface.kind==ThermalHistory::BoundaryKind::Bed ||
                surface.kind<ThermalHistory::BoundaryKind::Ambient || surface.kind>ThermalHistory::BoundaryKind::Exhaust ||
                surface.face_mask==0 || (surface.face_mask & ~0x3fu)!=0 ||
                (surface.bounds && (!surface.bounds->minimum_m.allFinite() || !surface.bounds->maximum_m.allFinite() ||
                    (surface.bounds->minimum_m.array()>surface.bounds->maximum_m.array()).any())))
                throw std::invalid_argument("Invalid exposed surface exchange.");
    }
    if (scenarios.front().time_s>activation[order.front()].time_s)
        throw std::invalid_argument("Initial boundary conditions are missing.");
    const std::size_t missing=std::numeric_limits<std::size_t>::max();
    std::vector<std::size_t> active_id(grid.voxels.size(),missing), active;
    std::vector<std::size_t> released_by(grid.voxels.size(), missing);
    std::size_t position=0, boundary=0, release_position=0, entries=0, event_count=0;
    const auto charge=[&](std::size_t count) {
        if (count>limits.max_schedule_entries-entries) throw ResourceLimit("Thermal schedule exceeds its storage/work budget.");
        entries+=count;
    };
    const auto next_event = [&](ThermalHistory::DepositionEvent &event) {
        if (position==order.size() && boundary==scenarios.size() && release_position==bed_releases.size()) return false;
        if (checkpoint) checkpoint();
        const double deposition=position<order.size()?activation[order[position]].time_s:std::numeric_limits<double>::infinity();
        const double change=boundary<scenarios.size()?scenarios[boundary].time_s:std::numeric_limits<double>::infinity();
        const double release_time=release_position<bed_releases.size()?bed_releases[release_position].time_s:std::numeric_limits<double>::infinity();
        const double time=std::min({deposition,change,release_time});
        if (time>end) return false;
        if (event_count++>=limits.thermal.max_events) throw ResourceLimit("Thermal event budget exceeded.");
        event = {};
        event.time_s=time;
        while (boundary<scenarios.size() && scenarios[boundary].time_s<=time) ++boundary;
        const auto &scenario=scenarios[boundary-1];
        while (release_position<bed_releases.size() && bed_releases[release_position].time_s<=time) {
            charge(1);
            released_by[bed_releases[release_position].voxel] = release_position;
            ++release_position;
        }
        const std::size_t old_count=active.size();
        while (position<order.size() && activation[order[position]].time_s==time) {
            const auto index=order[position++];
            const auto &material=materials[grid.voxels[index].material];
            charge(1);
            active_id[index]=active.size(); active.push_back(index);
            event.deposited_cells.push_back({volume,material.volumetric_heat_capacity_j_m3_k,activation[index].temperature_k,0.});
        }
        // Add each new shared face once. Empty cells remain true thermal gaps.
        for (std::size_t id=old_count;id<active.size();++id) {
            const auto &voxel=grid.voxels[active[id]];
            for (int axis=0;axis<3;++axis) for (int sign : {-1,1}) {
                auto neighbor=voxel.index; neighbor[axis]+=sign;
                const auto found=occupied.find(neighbor);
                if (found==occupied.end() || active_id[found->second]>=id) continue;
                const auto other=active_id[found->second];
                const auto material=voxel.material, adjacent=grid.voxels[found->second].material;
                const double length=grid.spacing_m[axis], area=volume/length;
                const double denominator=.5*length/materials[material].conductivity_w_m_k+
                    .5*length/materials[adjacent].conductivity_w_m_k+(material==adjacent?0.:resistance);
                charge(1);
                event.added_connections.push_back({other,id,area/denominator});
            }
        }
        // Rebuild exposure at each activation/boundary event, removing buried faces.
        charge(active.size());
        for (std::size_t id=0;id<active.size();++id) {
            if (checkpoint && id%128==0) checkpoint();
            const auto &voxel=grid.voxels[active[id]];
            for (int axis=0;axis<3;++axis) for (int sign : {-1,1}) {
                auto neighbor=voxel.index; neighbor[axis]+=sign;
                const auto found=occupied.find(neighbor);
                if (found!=occupied.end() && active_id[found->second]!=missing) continue;
                const double area=volume/grid.spacing_m[axis];
                const double half_cell_resistance=.5*grid.spacing_m[axis]/materials[voxel.material].conductivity_w_m_k;
                const auto conductance=[&](double h) {
                    return h==0. ? 0. : area/(half_cell_resistance+1./h);
                };
                if (axis==2 && sign==-1 && voxel.index[2]==0) {
                    charge(1);
                    const auto released = released_by[active[id]];
                    if (released == missing)
                        event.exposed_boundaries.push_back({id,conductance(scenario.bed_coefficient_w_m2_k),scenario.bed_temperature_k,ThermalHistory::BoundaryKind::Bed});
                    else {
                        const auto &gap = bed_releases[released];
                        event.exposed_boundaries.push_back({id,conductance(gap.coefficient_w_m2_k),gap.temperature_k,ThermalHistory::BoundaryKind::Ambient});
                    }
                } else {
                    Eigen::Vector3d center = grid.origin_m + grid.spacing_m.cwiseProduct(
                        Eigen::Vector3d(double(voxel.index[0])+.5,double(voxel.index[1])+.5,double(voxel.index[2])+.5));
                    center[axis] += .5*sign*grid.spacing_m[axis];
                    if (!center.allFinite()) throw std::invalid_argument("Thermal face coordinates overflow.");
                    double combined_h=0., combined_temperature=0.;
                    std::size_t matched=0;
                    auto kind=ThermalHistory::BoundaryKind::Ambient;
                    // Parallel surface exchanges share ONE solid half-cell
                    // resistance. Selection never changes geometry or mass.
                    charge(scenario.exposed.size());
                    for (const auto &surface : scenario.exposed) {
                        if ((surface.face_mask & (1u << (2*axis+(sign>0 ? 1 : 0))))==0 ||
                            (surface.bounds && ((center.array()<surface.bounds->minimum_m.array()).any() ||
                                                (center.array()>surface.bounds->maximum_m.array()).any())) ||
                            surface.coefficient_w_m2_k==0.) continue;
                        const double next_h=combined_h+surface.coefficient_w_m2_k;
                        if (!std::isfinite(next_h)) throw std::invalid_argument("Surface exchange coefficient overflow.");
                        combined_temperature=combined_temperature*(combined_h/next_h)+
                            surface.temperature_k*(surface.coefficient_w_m2_k/next_h);
                        combined_h=next_h;
                        ++matched;
                        kind=surface.kind;
                    }
                    if (!std::isfinite(combined_temperature)) throw std::invalid_argument("Surface reservoir temperature overflow.");
                    if (combined_h==0.) continue;
                    charge(1);
                    // Mixed surface reservoirs share a solid resistance; do not
                    // mislabel their net heat as natural ambient convection.
                    event.exposed_boundaries.push_back({id,conductance(combined_h),combined_temperature,
                        matched==1?kind:ThermalHistory::BoundaryKind::Mixed});
                }
            }
        }
        return true;
    };
    SimulationResult result;
    result.thermal=ThermalHistory::evolve_generated(next_event,end,limits.thermal,checkpoint);
    if (result.thermal.state!=ThermalHistory::State::Converged || active.empty()) return result;
    result.active_voxel_indices=active;
    std::vector<Voxel> mechanical_voxels;
    mechanical_voxels.reserve(active.size());
    for (std::size_t id=0;id<active.size();++id) {
        auto voxel=grid.voxels[active[id]];
        const auto &material=materials[voxel.material];
        voxel.eigenstrain=thermal_eigenstrain(material.expansion,material.local_to_world,
            material.stress_free_temperature_k,result.thermal.temperatures_k[id]);
        mechanical_voxels.push_back(voxel);
    }
    result.mesh=mesh_voxels(mechanical_voxels,grid.origin_m,grid.spacing_m,limits.max_mechanical_nodes,checkpoint);
    // mesh_voxels emits six tetrahedra per input voxel, in input order.
    if (result.mesh.elements.size()!=6*active.size()) throw std::logic_error("Unexpected voxel tetrahedralization.");
    for (std::size_t id=0;id<active.size();++id)
        for (std::size_t tet=0;tet<6;++tet)
            result.mesh.elements[6*id+tet].elasticity_pa=materials[grid.voxels[active[id]].material].elasticity_pa;
    if (release_position != 0) {
        // mesh_voxels emits two bed triangles per z=0 voxel, in input order.
        // Use stable voxel identity, never a nearest-coordinate match across
        // distinct material islands or a guessed correspondence after refinement.
        std::vector<std::size_t> retained;
        std::size_t face = 0;
        for (const auto index : active) {
            if (grid.voxels[index].index[2] != 0) continue;
            if (released_by[index] == missing) {
                retained.push_back(face);
                retained.push_back(face + 1);
            }
            face += 2;
        }
        if (face != result.mesh.bed_faces.size()) throw std::logic_error("Unexpected voxel bed-face ordering.");
        result.partial_mechanical = solve_partially_bonded_contact(result.mesh.nodes_m, result.mesh.elements,
            result.mesh.bed_faces, retained, grid.origin_m.z(), contact_limits, checkpoint);
        return result;
    }
    if (!limits.reduction) {
        result.mechanical=solve_bonded_bed(result.mesh.nodes_m,result.mesh.elements,result.mesh.bed_faces,grid.origin_m.z(),checkpoint);
        return result;
    }
    result.reduction.emplace();
    auto& study = *result.reduction;
    const auto& settings = *limits.reduction;
    std::vector<std::size_t> pinned;
    std::vector<Constraint> constraints;
    for (std::size_t node = 0; node < result.mesh.nodes_m.size(); ++node)
        if (result.mesh.nodes_m[node].z() == grid.origin_m.z()) {
            pinned.push_back(node);
            for (unsigned axis = 0; axis < 3; ++axis) constraints.push_back({node, axis});
        }
    unsigned stable_comparisons = 0;
    try {
        for (unsigned stride : settings.strides) {
            const auto basis = occupied_line_basis(result.mesh.nodes_m, result.mesh.elements, pinned,
                settings.axis, stride, settings.solve, checkpoint);
            auto next = solve_reduced(result.mesh.nodes_m, result.mesh.elements, constraints,
                basis.nodes_m, basis.rows, settings.solve, checkpoint);
            if (study.latest) {
                double difference = 0., scale = 0.;
                for (std::size_t node = 0; node < result.mesh.nodes_m.size(); ++node) {
                    difference += (next.recovered.displacement_m[node] - study.latest->recovered.displacement_m[node]).squaredNorm();
                    scale += next.recovered.displacement_m[node].squaredNorm();
                }
                study.displacement_relative_change = scale > 0. ? std::sqrt(difference / scale) : (difference == 0. ? 0. : std::numeric_limits<double>::infinity());
                const double energy = next.recovered.elastic_energy_j;
                const double energy_difference = std::abs(energy - study.latest->recovered.elastic_energy_j);
                study.energy_relative_change = energy > 0. ? energy_difference / energy : (energy_difference == 0. ? 0. : std::numeric_limits<double>::infinity());
                // Smaller strides can yield the same basis on short lines. Such
                // duplicate solves are not independent refinement evidence.
                if (next.free_reduced_dofs > study.latest->free_reduced_dofs)
                    stable_comparisons = study.displacement_relative_change <= settings.displacement_relative_tolerance &&
                        study.energy_relative_change <= settings.energy_relative_tolerance ? stable_comparisons + 1 : 0;
            }
            study.latest = std::move(next);
            study.solved_strides.push_back(stride);
            if (stride == 1) {
                study.state = ReductionStudy::State::FineReference;
                break;
            }
            if (stable_comparisons >= 2) {
                study.state = ReductionStudy::State::Stabilized;
                break;
            }
        }
    } catch (const ReductionBudget&) {
        study.state = ReductionStudy::State::BudgetExceeded;
    }
    return result;
}
PeelingHistoryResult simulate_peeling_history(const DepositionGrid &grid,
    const std::vector<CellActivation> &activation,
    const std::vector<ThermalMechanicalMaterial> &materials,
    const std::vector<BoundaryScenario> &scenarios, double resistance,
    const std::vector<double> &times, const std::vector<GridReleasePatch> &patches,
    const PeelingHistoryLimits &limits, const std::function<void()> &checkpoint)
{
    if (limits.max_observations == 0 || limits.max_observations > 16 ||
        limits.max_equilibrium_solves == 0 || limits.max_equilibrium_solves > 256 ||
        times.empty() || patches.empty() || activation.size() != grid.voxels.size() || limits.simulation.reduction)
        throw std::invalid_argument("Invalid peeling history inputs or limits.");
    if (times.size() > limits.max_observations || patches.size() > 32 || grid.voxels.size() > 200000)
        throw ResourceLimit("Peeling history exceeds its observation or geometry budget.");
    double previous = -1.;
    for (const auto time : times) {
        if (!std::isfinite(time) || time < 0. || time <= previous)
            throw std::invalid_argument("Peeling observation times must be finite and strictly increasing.");
        previous = time;
    }
    std::vector<std::size_t> owner(grid.voxels.size(), patches.size());
    for (std::size_t p = 0; p < patches.size(); ++p) {
        const auto &patch = patches[p];
        if (patch.voxels.empty() || !std::isfinite(patch.fracture_energy_j_m2) || patch.fracture_energy_j_m2 < 0. ||
            !std::isfinite(patch.gap_coefficient_w_m2_k) || patch.gap_coefficient_w_m2_k < 0. ||
            !std::isfinite(patch.gap_temperature_k) || patch.gap_temperature_k <= 0.)
            throw std::invalid_argument("Invalid peeling history patch properties.");
        for (const auto voxel : patch.voxels) {
            if (checkpoint) checkpoint();
            if (voxel >= owner.size() || owner[voxel] != patches.size() || grid.voxels[voxel].index[2] != 0 ||
                !std::isfinite(activation[voxel].time_s) || activation[voxel].time_s > times.front())
                throw std::invalid_argument("Peeling patches must be disjoint bottom cells active at the first observation.");
            owner[voxel] = p;
        }
    }
    PeelingHistoryResult result;
    std::vector<bool> released(patches.size(), false);
    const auto remaining_solves = [&] {
        if (result.equilibrium_solves >= limits.max_equilibrium_solves)
            throw ResourceLimit("Peeling history exhausted its global equilibrium budget.");
        return limits.max_equilibrium_solves - result.equilibrium_solves;
    };
    for (const auto time : times) {
        if (checkpoint) checkpoint();
        auto contact_limits = limits.release.contact;
        contact_limits.max_iterations = std::min(contact_limits.max_iterations, remaining_solves());
        auto snapshot = simulate_deposition_grid(grid, activation, materials, scenarios, resistance, time,
            limits.simulation, checkpoint, result.releases, contact_limits);
        if (snapshot.thermal.state != ThermalHistory::State::Converged ||
            (!snapshot.mechanical && !snapshot.partial_mechanical))
            throw std::runtime_error("Peeling history thermal or mechanical snapshot is unresolved.");
        result.equilibrium_solves += snapshot.partial_mechanical ? snapshot.partial_mechanical->iterations : 1;
        result.maximum_accepted_thermal_interval_s = std::max(result.maximum_accepted_thermal_interval_s,
            snapshot.thermal.maximum_accepted_interval_s);
        std::vector<BedReleasePatch> face_patches(patches.size());
        for (std::size_t p = 0; p < patches.size(); ++p)
            face_patches[p].fracture_energy_j_m2 = patches[p].fracture_energy_j_m2;
        std::vector<std::size_t> retained;
        std::size_t face = 0;
        for (const auto voxel : snapshot.active_voxel_indices) {
            if (grid.voxels[voxel].index[2] != 0) continue;
            const auto p = owner[voxel];
            if (p != patches.size()) {
                face_patches[p].faces.push_back(face);
                face_patches[p].faces.push_back(face + 1);
            }
            if (p == patches.size() || !released[p]) {
                retained.push_back(face);
                retained.push_back(face + 1);
            }
            face += 2;
        }
        if (face != snapshot.mesh.bed_faces.size()) throw std::logic_error("Unexpected peeling bed-face mapping.");
        auto release_limits = limits.release;
        release_limits.max_equilibrium_solves = std::min(release_limits.max_equilibrium_solves, remaining_solves());
        auto evolution = evolve_bed_release(snapshot.mesh.nodes_m, snapshot.mesh.elements, snapshot.mesh.bed_faces,
            retained, face_patches, grid.origin_m.z(), release_limits, checkpoint);
        result.equilibrium_solves += evolution.equilibrium_solves;
        for (const auto &step : evolution.steps) {
            const auto &patch = patches[step.patch];
            released[step.patch] = true;
            for (const auto voxel : patch.voxels)
                result.releases.push_back({time, voxel, patch.gap_coefficient_w_m2_k, patch.gap_temperature_k});
        }
        const auto &before_release = snapshot.mechanical ? *snapshot.mechanical : snapshot.partial_mechanical->bed.bonded;
        result.samples.push_back({time, evolution.equilibrium.released_area_total_m2,
            evolution.equilibrium.maximum_lift_m, evolution.equilibrium.bonded.equilibrium.elastic_energy_j,
            evolution.steps.size(), evolution.equilibrium.bonded.tensile_force_n,
            before_release.tensile_force_n, before_release.equilibrium.elastic_energy_j});
        // Do not retain a stale fully bonded field after accepting new damage.
        snapshot.mechanical.reset();
        snapshot.partial_mechanical.reset();
        if (time == times.back()) {
            result.final_snapshot = std::move(snapshot);
            result.final_release = std::move(evolution);
        }
    }
    return result;
}
namespace {
bool peeling_measures_close(const PeelingHistoryResult &a, const PeelingHistoryResult &b,
    const PeelingTimeStudyLimits &limits)
{
    const auto close = [&](double x, double y, double absolute) {
        return std::isfinite(x) && std::isfinite(y) &&
            std::abs(x-y) <= absolute + limits.relative_tolerance * std::max(std::abs(x),std::abs(y));
    };
    std::size_t next = 0;
    for (const auto &sample : a.samples) {
        while (next < b.samples.size() && b.samples[next].time_s < sample.time_s) ++next;
        if (next == b.samples.size() || b.samples[next].time_s != sample.time_s ||
            !close(sample.maximum_lift_m,b.samples[next].maximum_lift_m,limits.lift_tolerance_m) ||
            !close(sample.released_area_m2,b.samples[next].released_area_m2,limits.area_tolerance_m2) ||
            !close(sample.tensile_force_n,b.samples[next].tensile_force_n,limits.force_tolerance_n) ||
            !close(sample.pre_release_force_n,b.samples[next].pre_release_force_n,limits.force_tolerance_n) ||
            !close(sample.pre_release_energy_j,b.samples[next].pre_release_energy_j,limits.energy_tolerance_j) ||
            !close(sample.elastic_energy_j,b.samples[next].elastic_energy_j,limits.energy_tolerance_j)) return false;
    }
    const auto peaks = [](const PeelingHistoryResult &run) {
        std::array<double,3> peak{};
        for (const auto &sample : run.samples) {
            peak[0] = std::max(peak[0],sample.maximum_lift_m);
            peak[1] = std::max({peak[1],sample.elastic_energy_j,sample.pre_release_energy_j});
            peak[2] = std::max({peak[2],sample.tensile_force_n,sample.pre_release_force_n});
        }
        return peak;
    };
    const auto peak_a = peaks(a), peak_b = peaks(b);
    // Net heat can stabilize while bed heating and fan cooling errors cancel.
    // Require each physical boundary channel to stabilize independently.
    for (std::size_t kind = 0; kind < ThermalHistory::boundary_kind_count; ++kind)
        if (!close(a.final_snapshot.thermal.boundary_energy_by_kind_j[kind],
                   b.final_snapshot.thermal.boundary_energy_by_kind_j[kind], limits.energy_tolerance_j)) return false;
    return close(peak_a[0],peak_b[0],limits.lift_tolerance_m) &&
        close(peak_a[1],peak_b[1],limits.energy_tolerance_j) &&
        close(peak_a[2],peak_b[2],limits.force_tolerance_n) &&
        close(a.final_snapshot.thermal.final_energy_j,b.final_snapshot.thermal.final_energy_j,limits.energy_tolerance_j) &&
        close(a.final_snapshot.thermal.boundary_energy_j,b.final_snapshot.thermal.boundary_energy_j,limits.energy_tolerance_j);
}
} // namespace

PeelingTimeStudy study_peeling_time(const DepositionGrid &grid,
    const std::vector<CellActivation> &activation,
    const std::vector<ThermalMechanicalMaterial> &materials,
    const std::vector<BoundaryScenario> &scenarios, double resistance,
    const std::vector<double> &times, const std::vector<GridReleasePatch> &patches,
    const PeelingTimeStudyLimits &limits, const std::function<void()> &checkpoint)
{
    const auto positive = [](double x) { return std::isfinite(x) && x > 0.; };
    if (times.size() < 2 || times.size() > 4 || limits.max_history_runs == 0 || limits.max_history_runs > 6 ||
        limits.max_total_observations == 0 || limits.max_total_observations > 64 ||
        limits.max_total_equilibrium_solves == 0 || limits.max_total_equilibrium_solves > 256 ||
        !positive(limits.relative_tolerance) || limits.relative_tolerance >= 1. ||
        !positive(limits.temperature_tolerance_k) || !positive(limits.release_time_tolerance_s) ||
        !positive(limits.lift_tolerance_m) || !positive(limits.area_tolerance_m2) || !positive(limits.energy_tolerance_j) ||
        !positive(limits.force_tolerance_n) ||
        !positive(limits.run.simulation.thermal.max_step_temperature_error_k))
        throw std::invalid_argument("Invalid peeling time-convergence criteria.");
    const auto histories_close = [&](const PeelingHistoryResult &a, const PeelingHistoryResult &b) {
        if (checkpoint) checkpoint();
        if (a.final_snapshot.active_voxel_indices != b.final_snapshot.active_voxel_indices ||
            a.final_snapshot.thermal.temperatures_k.size() != b.final_snapshot.thermal.temperatures_k.size()) return false;
        for (std::size_t i=0; i<a.final_snapshot.thermal.temperatures_k.size(); ++i)
            if (!std::isfinite(a.final_snapshot.thermal.temperatures_k[i]) ||
                !std::isfinite(b.final_snapshot.thermal.temperatures_k[i]) ||
                std::abs(a.final_snapshot.thermal.temperatures_k[i]-b.final_snapshot.thermal.temperatures_k[i]) > limits.temperature_tolerance_k)
                return false;
        const auto &displacement_a = a.final_release.equilibrium.bonded.equilibrium.displacement_m;
        const auto &displacement_b = b.final_release.equilibrium.bonded.equilibrium.displacement_m;
        if (displacement_a.size() != displacement_b.size()) return false;
        for (std::size_t i=0; i<displacement_a.size(); ++i)
            if (!displacement_a[i].allFinite() || !displacement_b[i].allFinite() ||
                (displacement_a[i]-displacement_b[i]).norm() > limits.lift_tolerance_m +
                    limits.relative_tolerance * std::max(displacement_a[i].norm(),displacement_b[i].norm())) return false;
        std::vector<double> release_a(grid.voxels.size(),-1.), release_b(grid.voxels.size(),-1.);
        for (const auto &event : a.releases) release_a[event.voxel] = event.time_s;
        for (const auto &event : b.releases) release_b[event.voxel] = event.time_s;
        for (std::size_t i=0; i<release_a.size(); ++i)
            if ((release_a[i] < 0.) != (release_b[i] < 0.) ||
                std::abs(release_a[i]-release_b[i]) > limits.release_time_tolerance_s) return false;
        return peeling_measures_close(a,b,limits);
    };
    PeelingTimeStudy study;
    std::vector<double> schedule = times;
    std::optional<PeelingHistoryResult> previous;
    for (unsigned level=0; level<3; ++level) {
        auto run_limits = limits.run;
        const auto run = [&] {
            if (checkpoint) checkpoint();
            if (study.completed_runs >= limits.max_history_runs ||
                schedule.size() > limits.max_total_observations-study.completed_observations ||
                study.equilibrium_solves >= limits.max_total_equilibrium_solves)
                throw ResourceLimit("Peeling time study exhausted its global budget.");
            run_limits.max_equilibrium_solves = std::min(run_limits.max_equilibrium_solves,
                limits.max_total_equilibrium_solves-study.equilibrium_solves);
            auto result = simulate_peeling_history(grid,activation,materials,scenarios,resistance,schedule,patches,run_limits,checkpoint);
            ++study.completed_runs;
            study.completed_observations += static_cast<unsigned>(schedule.size());
            study.equilibrium_solves += result.equilibrium_solves;
            return result;
        };
        auto coarse = run();
        run_limits.simulation.thermal.max_step_s *= .5;
        if (coarse.maximum_accepted_thermal_interval_s > 0.)
            run_limits.simulation.thermal.max_step_s = std::min(run_limits.simulation.thermal.max_step_s,
                .5*coarse.maximum_accepted_thermal_interval_s);
        run_limits.simulation.thermal.max_step_temperature_error_k *= .25;
        auto fine = run();
        if (!histories_close(coarse,fine) || (previous && !histories_close(*previous,fine))) return study;
        if (level == 2) {
            study.temporally_qualified = std::move(fine);
            return study;
        }
        previous = std::move(fine);
        std::vector<double> refined;
        for (std::size_t i=0; i<schedule.size(); ++i) {
            if (i != 0) {
                const double midpoint = schedule[i-1] + .5*(schedule[i]-schedule[i-1]);
                if (!(midpoint > schedule[i-1] && midpoint < schedule[i]))
                    throw std::invalid_argument("Peeling observation interval cannot be refined.");
                refined.push_back(midpoint);
            }
            refined.push_back(schedule[i]);
        }
        schedule = std::move(refined);
    }
    return study;
}
PeelingSpaceStudy study_peeling_space(const DepositionGrid &grid,
    const std::vector<CellActivation> &activation,
    const std::vector<ThermalMechanicalMaterial> &materials,
    const std::vector<BoundaryScenario> &scenarios, double resistance,
    const std::vector<double> &times, const std::vector<GridReleasePatch> &patches,
    const PeelingSpaceStudyLimits &limits, const std::function<void()> &checkpoint)
{
    if (limits.max_history_runs == 0 || limits.max_history_runs > 18 ||
        limits.max_total_observations == 0 || limits.max_total_observations > 192 ||
        limits.max_total_equilibrium_solves == 0 || limits.max_total_equilibrium_solves > 1024 ||
        activation.size() != grid.voxels.size())
        throw std::invalid_argument("Invalid peeling spatial-study limits.");
    if (grid.voxels.size() > 200000) throw ResourceLimit("Peeling spatial-study grid exceeds its budget.");
    PeelingSpaceStudy study;
    auto current_grid = grid;
    auto current_activation = activation;
    auto current_patches = patches;
    std::optional<PeelingHistoryResult> previous;
    std::vector<double> previous_release_times;
    unsigned subdivisions = 1;
    for (unsigned level=0; level<3; ++level) {
        if (checkpoint) checkpoint();
        if (study.completed_runs >= limits.max_history_runs || study.completed_observations >= limits.max_total_observations ||
            study.equilibrium_solves >= limits.max_total_equilibrium_solves)
            throw ResourceLimit("Peeling spatial study exhausted its global budget.");
        auto time_limits = limits.time;
        time_limits.max_history_runs = std::min(time_limits.max_history_runs,limits.max_history_runs-study.completed_runs);
        time_limits.max_total_observations = std::min(time_limits.max_total_observations,
            limits.max_total_observations-study.completed_observations);
        time_limits.max_total_equilibrium_solves = std::min(time_limits.max_total_equilibrium_solves,
            limits.max_total_equilibrium_solves-study.equilibrium_solves);
        auto temporal = study_peeling_time(current_grid,current_activation,materials,scenarios,resistance,times,
            current_patches,time_limits,checkpoint);
        study.completed_runs += temporal.completed_runs;
        study.completed_observations += temporal.completed_observations;
        study.equilibrium_solves += temporal.equilibrium_solves;
        if (!temporal.temporally_qualified) return study;
        ++study.completed_levels;
        auto current = std::move(*temporal.temporally_qualified);
        std::vector<double> cell_release_times(current_grid.voxels.size(),-1.);
        for (const auto &event : current.releases) cell_release_times[event.voxel] = event.time_s;
        std::vector<double> release_times;
        for (const auto &patch : current_patches) {
            const double time = cell_release_times[patch.voxels.front()];
            for (const auto voxel : patch.voxels)
                if (cell_release_times[voxel] != time) throw std::logic_error("Peeling patch released inconsistently.");
            release_times.push_back(time);
        }
        if (previous) {
            if (!peeling_measures_close(*previous,current,limits.time)) return study;
            for (std::size_t p=0; p<release_times.size(); ++p)
                if ((release_times[p] < 0.) != (previous_release_times[p] < 0.) ||
                    std::abs(release_times[p]-previous_release_times[p]) > limits.time.release_time_tolerance_s) return study;
            // Each parent cell has eight equal-volume children with the same
            // activation. Compare cell averages, not unrelated node indices.
            std::vector<double> average(current_grid.voxels.size()/8,0.);
            std::vector<unsigned> count(average.size(),0);
            for (std::size_t i=0; i<current.final_snapshot.active_voxel_indices.size(); ++i) {
                const auto parent = current.final_snapshot.active_voxel_indices[i]/8;
                average[parent] += current.final_snapshot.thermal.temperatures_k[i]/8.;
                ++count[parent];
            }
            for (std::size_t i=0; i<previous->final_snapshot.active_voxel_indices.size(); ++i) {
                const auto parent = previous->final_snapshot.active_voxel_indices[i];
                if (count[parent] != 8 || !std::isfinite(average[parent]) ||
                    std::abs(average[parent]-previous->final_snapshot.thermal.temperatures_k[i]) > limits.time.temperature_tolerance_k)
                    return study;
            }
        }
        if (level == 2) {
            study.discretization_qualified = QualifiedPeelingGrid{std::move(current_grid),std::move(current),subdivisions};
            return study;
        }
        previous = std::move(current);
        previous_release_times = std::move(release_times);
        const auto cell_limit = std::min(std::size_t(200000),limits.time.run.simulation.thermal.solve.max_cells);
        if (current_grid.voxels.size() > cell_limit/8)
            throw ResourceLimit("Refined peeling grid exceeds the cell budget.");
        DepositionGrid refined{current_grid.origin_m,current_grid.spacing_m*.5,{}};
        std::vector<CellActivation> refined_activation;
        refined.voxels.reserve(8*current_grid.voxels.size());
        refined_activation.reserve(8*current_activation.size());
        for (std::size_t i=0; i<current_grid.voxels.size(); ++i) {
            if (checkpoint && i%128==0) checkpoint();
            for (int z=0; z<2; ++z) for (int y=0; y<2; ++y) for (int x=0; x<2; ++x) {
                auto child = current_grid.voxels[i];
                const std::array<int,3> offset{x,y,z};
                for (unsigned axis=0; axis<3; ++axis) {
                    const long long coordinate = 2LL*child.index[axis]+offset[axis];
                    if (coordinate < std::numeric_limits<int>::min() || coordinate > std::numeric_limits<int>::max())
                        throw std::invalid_argument("Refined peeling coordinates overflow.");
                    child.index[axis] = static_cast<int>(coordinate);
                }
                refined.voxels.push_back(std::move(child));
                refined_activation.push_back(current_activation[i]);
            }
        }
        for (auto &patch : current_patches) {
            std::vector<std::size_t> children;
            children.reserve(4*patch.voxels.size());
            for (const auto parent : patch.voxels)
                for (std::size_t child=0; child<4; ++child) children.push_back(8*parent+child);
            patch.voxels = std::move(children);
        }
        current_grid = std::move(refined);
        current_activation = std::move(refined_activation);
        subdivisions *= 2;
    }
    return study;
}
} // namespace Slic3r::Thermoelastic
