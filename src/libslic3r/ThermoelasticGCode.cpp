#include "ThermoelasticGCode.hpp"

#include "ExtrusionEntity.hpp"
#include "GCode/GCodeProcessor.hpp"
#include "GCodeReader.hpp"
#include "ThermoelasticField.hpp"
#include "ThermoelasticMesh.hpp"
#include "ThermalAirflow.hpp"
#include "Config.hpp"
#include "PrintConfig.hpp"
#include "WarpMitigation.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <ios>
#include <istream>
#include <limits>
#include <new>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Slic3r::Thermoelastic {

namespace {
AirflowTraceSummary trace_airflow(const AirflowTraceInput &input,
    const std::function<void()> &checkpoint)
{
    if (input.max_bytes == 0 || input.max_bytes > 256 * 1024 * 1024 ||
        input.max_events == 0 || input.max_events > 100000)
        throw std::invalid_argument("Invalid airflow trace budget.");
    ThermalHistory::AirflowTimeline timeline(input.channels, input.policy, input.max_events);
    AirflowTraceSummary result;
    GCodeReader reader;
    std::array<char, 65537> buffer;
    std::size_t bytes = 0;
    result.state = AirflowTraceState::Complete;
    while (true) {
        if (checkpoint && result.source_lines % 256 == 0) checkpoint();
        try {
            input.source.getline(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        } catch (const std::ios_base::failure &) {
            // Inspect stream flags below, including exception-enabled EOF.
        }
        const auto count = static_cast<std::size_t>(input.source.gcount());
        if (input.source.bad()) { result.state = AirflowTraceState::ReadError; break; }
        if (input.source.fail() && !input.source.eof()) {
            result.state = count == buffer.size() - 1 ? AirflowTraceState::BudgetExceeded : AirflowTraceState::ReadError;
            break;
        }
        if (count == 0 && input.source.eof()) break;
        if (count > input.max_bytes - bytes) { result.state = AirflowTraceState::BudgetExceeded; break; }
        bytes += count;
        ++result.source_lines;
        try {
            reader.parse_line(std::string(buffer.data()), [&](auto &, const auto &line) {
                timeline.consume(line, result.source_lines);
                if ((line.cmd_is("G1") || line.cmd_is("G2") || line.cmd_is("G3")) && line.has('E'))
                    for (std::size_t role = 0; role < result.unresolved_during_motion.size(); ++role)
                        result.unresolved_during_motion[role] |= !timeline.state().pwm[role].has_value();
            });
        } catch (const std::runtime_error &) {
            result.state = AirflowTraceState::BudgetExceeded;
            break;
        }
        if (input.source.eof()) break;
    }
    result.final_state = timeline.state();
    result.events = timeline.events().size();
    if (result.state == AirflowTraceState::Complete)
        result.decoded_events = timeline.events();
    return result;
}

bool deposited(const GCodeProcessorResult::MoveVertex &move)
{
    return move.type == EMoveType::Extrude && move.extrusion_role != erWipeTower &&
        move.extrusion_role != erCustom && move.extrusion_role != erNone;
}

// Static analysis does not require timing, but uses the same physical endpoint
// convention as timed_deposition. Display markers must not alter bead geometry
// or stress sample locations.
template<class Visitor>
void visit_deposited_moves(const GCodeProcessorResult &gcode, Visitor visit,
    const std::function<void()> &checkpoint)
{
    std::optional<std::size_t> previous;
    for (std::size_t i = 0; i < gcode.moves.size(); ++i) {
        if (checkpoint && i % 1024 == 0) checkpoint();
        const auto &move = gcode.moves[i];
        if (move.internal_only || move.type == EMoveType::Seam) continue;
        if (deposited(move)) {
            if (!previous)
                throw std::invalid_argument("Mechanical deposition has no physical starting point.");
            visit(i, gcode.moves[*previous].position, move);
        }
        previous = i;
    }
}

std::vector<float> sample_moves(const GCodeProcessorResult &gcode, const RefinementResult &analysis,
    const std::function<void()> &checkpoint)
{
    const Field field(analysis.mesh.nodes_m, analysis.mesh.elements, analysis.bed.equilibrium);
    std::vector<float> samples(gcode.moves.size(), -1.f);
    visit_deposited_moves(gcode, [&](std::size_t i, const Eigen::Vector3f &start,
                                   const GCodeProcessorResult::MoveVertex &move) {
        Eigen::Vector3d center = (start.cast<double>() + move.position.cast<double>()) * .0005;
        center.z() -= double(move.height) * .0005;
        const auto sample = field.sample(center);
        if (!sample) return;
        const double stress = std::max(0., sample->maximum_principal_pa) * 1e-6;
        if (!std::isfinite(stress) || stress > std::numeric_limits<float>::max())
            throw std::runtime_error("Mechanical preview stress is not representable.");
        samples[i] = float(stress);
    }, checkpoint);
    return samples;
}
} // namespace

ConfiguredAnalysis analyze_configured_gcode(GCodeProcessorResult &gcode,
    const DynamicPrintConfig &config, const AnalysisLimits &limits,
    const std::function<void(unsigned, unsigned)> &progress,
    const std::function<void()> &checkpoint,
    const std::vector<std::optional<MaterialScenario>> &material_scenarios,
    const AirflowTraceInput *airflow_input)
{
    gcode.thermal_stress_mpa.clear();
    gcode.thermal_analysis = {};
    ConfiguredAnalysis result;
    if (airflow_input) result.airflow = trace_airflow(*airflow_input, checkpoint);
    // Bound factorization as well as rasterization. Raising this limit without
    // checking sparse-factor memory use can exhaust a desktop's swap space.
    if (limits.max_nodes < 4 || limits.max_nodes > 60000 || limits.max_segments == 0 ||
        limits.max_subdivisions == 0 || limits.max_subdivisions > 16 ||
        !std::isfinite(limits.relative_tolerance) || limits.relative_tolerance <= 0. ||
        limits.relative_tolerance >= 1.)
        throw std::invalid_argument("Invalid configured mechanical analysis limits.");
    if (material_scenarios.size() > 256)
        throw std::invalid_argument("Too many explicit mechanical material slots.");
    for (const auto &entry : material_scenarios) if (entry) {
        const auto &m = *entry;
        if (m.filament_type.empty() || m.provenance.empty() || m.provenance.size() > 4096 ||
            !std::isfinite(m.young_pa + m.poisson + m.cte_min + m.cte_max +
                           m.stress_free_temperature_c + m.maximum_chamber_temperature_c) ||
            m.young_pa <= 0. || m.poisson <= -1. || m.poisson >= .5 ||
            m.cte_min < 0. || m.cte_max < m.cte_min || m.maximum_chamber_temperature_c <= 0. ||
            m.stress_free_temperature_c < m.maximum_chamber_temperature_c)
            throw std::invalid_argument("Invalid explicit mechanical material scenario.");
    }
    bool cancelled = false;
    const auto check = [&] {
        if (checkpoint) try { checkpoint(); }
        catch (...) { cancelled = true; throw; }
    };
    try {
        check();
        const auto *types = config.option<ConfigOptionStrings>("filament_type");
        const auto *chambers = config.option<ConfigOptionInts>("chamber_temperature");
        const auto *nozzles = config.option<ConfigOptionInts>("nozzle_temperature");
        const auto *plate = config.option<ConfigOptionEnum<BedType>>("curr_bed_type");
        const auto *beds = plate ? config.option<ConfigOptionInts>(get_bed_temp_key(plate->value)) : nullptr;
        std::set<std::size_t> used;
        double width = std::numeric_limits<double>::infinity(), height = width, top = 0.;
        std::size_t count = 0;
        visit_deposited_moves(gcode, [&](std::size_t, const Eigen::Vector3f &,
                                       const GCodeProcessorResult::MoveVertex &move) {
            if (++count > limits.max_segments) throw ResourceLimit("Mechanical path budget exceeded.");
            used.insert(move.extruder_id);
            width = std::min(width, double(move.width));
            height = std::min(height, double(move.height));
            top = std::max(top, double(move.position.z()));
        }, check);
        if (used.empty() || !(width > 0. && height > 0. && top > 0.) || !std::isfinite(width + height + top))
            return result;
        const auto resin = *warp_resin_material("ASA");
        std::vector<MaterialScenario> properties(*used.rbegin() + 1);
        for (const auto id : used) {
            // The cited ASA constitutive study explicitly uses nu=0.35. That
            // modelling assumption must not silently become measured data for
            // TPU, filled resins, arbitrary blends or other polymer families.
            if (types && id < types->values.size() && id < material_scenarios.size() && material_scenarios[id] &&
                material_scenarios[id]->filament_type == types->values[id])
                properties[id] = *material_scenarios[id];
            else if (types && id < types->values.size() && types->values[id] == "ASA" &&
                     !(id < material_scenarios.size() && material_scenarios[id]))
                properties[id] = {"ASA", "Published ASA screening scenario", resin.modulus_max, .35,
                    resin.cte_min, resin.cte_max, resin.tg_max, resin.tg_min};
            else
                result.unsupported_filaments.push_back(id);
        }
        if (!result.unsupported_filaments.empty()) {
            result.state = AnalysisState::UnsupportedMaterial;
            return result;
        }
        if (!chambers || !nozzles || !beds || chambers->values.empty() ||
            nozzles->values.empty() || beds->values.empty()) {
            result.state = AnalysisState::InvalidThermalInputs;
            return result;
        }
        std::vector<Voxel> materials(properties.size(), Voxel{{0,0,0}, 0., 0.});
        for (const auto id : used) {
            const auto &property = properties[id];
            materials[id].young_pa = property.young_pa;
            materials[id].poisson = property.poisson;
            if (chambers->get_at(id) <= 0 || chambers->get_at(id) >= property.maximum_chamber_temperature_c ||
                nozzles->get_at(id) <= property.stress_free_temperature_c || beds->get_at(id) <= 0) {
                result.state = AnalysisState::InvalidThermalInputs;
                return result;
            }
        }
        std::vector<std::size_t> canonical(materials.size());
        result.scenario_parameters.push_back(1.); // Correlated CTE/temperature scenario version.
        for (const auto id : used) {
            const auto &property = properties[id];
            result.scenario_parameters.insert(result.scenario_parameters.end(), {
                double(id), property.young_pa, property.poisson, property.cte_min, property.cte_max,
                property.stress_free_temperature_c, property.maximum_chamber_temperature_c,
                double(chambers->get_at(id)), double(nozzles->get_at(id)), double(beds->get_at(id))});
            result.scenario_sources.push_back(property.filament_type);
            result.scenario_sources.push_back(property.provenance);
            canonical[id] = id;
            for (const auto previous : used) {
                if (previous >= id) break;
                const auto &a = properties[id], &b = properties[previous];
                if (a.filament_type == b.filament_type && a.provenance == b.provenance &&
                    a.young_pa == b.young_pa && a.poisson == b.poisson &&
                    a.cte_min == b.cte_min && a.cte_max == b.cte_max &&
                    a.stress_free_temperature_c == b.stress_free_temperature_c &&
                    a.maximum_chamber_temperature_c == b.maximum_chamber_temperature_c &&
                    chambers->get_at(id) == chambers->get_at(previous) &&
                    nozzles->get_at(id) == nozzles->get_at(previous) &&
                    beds->get_at(id) == beds->get_at(previous)) {
                    canonical[id] = canonical[previous];
                    break;
                }
            }
        }
        std::vector<DepositedSegment> segments;
        segments.reserve(count);
        visit_deposited_moves(gcode, [&](std::size_t, const Eigen::Vector3f &start,
                                       const GCodeProcessorResult::MoveVertex &move) {
            if (move.extrusion_role == erMixed)
                throw std::invalid_argument("Mixed extrusion roles cannot be resolved for mechanical analysis.");
            segments.push_back({start.cast<double>() * .001,
                move.position.cast<double>() * .001, double(move.width) * .001,
                double(move.height) * .001, canonical[move.extruder_id]});
        }, check);
        auto grid = voxelize_deposition(segments, materials, {.0005 * width, .0005 * width, .0005 * height},
            0., limits.max_nodes, check);
        const auto component_tops = component_top_z(grid, check);
        std::vector<float> envelope;
        // Uniform chamber temperature and a linear bed-to-chamber gradient,
        // each at the published CTE endpoints. These are prescribed scenarios,
        // not a solved transient temperature history or a universal worst case.
        for (unsigned scenario = 0; scenario < 4; ++scenario) {
            check();
            for (std::size_t i = 0; i < grid.voxels.size(); ++i) {
                auto &voxel = grid.voxels[i];
                const auto id = voxel.material;
                const auto &property = properties[id];
                const double reference = property.stress_free_temperature_c;
                const double chamber = double(chambers->get_at(id));
                const double fraction = std::clamp((voxel.index[2] + .5) * grid.spacing_m.z() / component_tops[i], 0., 1.);
                const double temperature = scenario < 2 ? chamber :
                    chamber + (std::min(double(beds->get_at(id)), reference) - chamber) * (1. - fraction);
                const double alpha = scenario % 2 == 0 ? property.cte_min : property.cte_max;
                voxel.eigenstrain.head<3>().setConstant(-alpha * std::max(0., reference - temperature));
            }
            const auto analysis = solve_refined_voxels(grid.voxels, grid.origin_m, grid.spacing_m,
                limits.max_nodes, limits.max_subdivisions, limits.relative_tolerance,
                [&](unsigned subdivision) { if (progress) progress(scenario, subdivision); }, check);
            if (!analysis.converged) {
                result.state = analysis.resource_limited ? AnalysisState::BudgetExceeded : AnalysisState::Unconverged;
                result.scenarios.clear();
                return result;
            }
            result.scenarios.push_back({analysis.bed.tensile_force_n, analysis.bed.equilibrium.elastic_energy_j});
            auto samples = sample_moves(gcode, analysis, check);
            if (envelope.empty()) envelope = std::move(samples);
            else for (std::size_t i = 0; i < envelope.size(); ++i)
                envelope[i] = envelope[i] < 0.f || samples[i] < 0.f ? -1.f : std::max(envelope[i], samples[i]);
        }
        check();
        gcode.thermal_stress_mpa = std::move(envelope);
        result.state = AnalysisState::Converged;
    } catch (const std::bad_alloc &) {
        if (cancelled) throw;
        result.state = AnalysisState::BudgetExceeded;
        result.scenarios.clear();
    } catch (const ResourceLimit &) {
        if (cancelled) throw;
        result.state = AnalysisState::BudgetExceeded;
        result.scenarios.clear();
    } catch (const std::runtime_error &) {
        if (cancelled) throw;
        result.state = AnalysisState::InvalidGeometry;
        result.scenarios.clear();
    } catch (const std::invalid_argument &) {
        if (cancelled) throw;
        result.state = AnalysisState::InvalidGeometry;
        result.scenarios.clear();
    }
    return result;
}

bool mechanically_improves(const ConfiguredAnalysis &baseline,
    const ConfiguredAnalysis &candidate, double relative_tolerance)
{
    if (!std::isfinite(relative_tolerance) || relative_tolerance <= 0. || relative_tolerance >= 1. ||
        baseline.state != AnalysisState::Converged || candidate.state != AnalysisState::Converged ||
        baseline.scenario_parameters.empty() || baseline.scenario_sources.empty() ||
        baseline.scenario_parameters != candidate.scenario_parameters ||
        baseline.scenario_sources != candidate.scenario_sources ||
        baseline.scenarios.empty() || baseline.scenarios.size() != candidate.scenarios.size()) return false;
    for (std::size_t i = 0; i < baseline.scenarios.size(); ++i) {
        const auto &a = baseline.scenarios[i], &b = candidate.scenarios[i];
        if (!std::isfinite(a.tensile_force_n + a.elastic_energy_j + b.tensile_force_n + b.elastic_energy_j) ||
            b.tensile_force_n < 0. || b.elastic_energy_j < 0. ||
            !(b.tensile_force_n * (1. + relative_tolerance) < a.tensile_force_n * (1. - relative_tolerance)) ||
            !(b.elastic_energy_j * (1. + relative_tolerance) < a.elastic_energy_j * (1. - relative_tolerance))) return false;
    }
    return true;
}

RefinementResult analyze_gcode(GCodeProcessorResult &gcode, const std::vector<Voxel> &materials,
    const Eigen::Vector3d &spacing_m, double bed_z_m, std::size_t max_nodes,
    unsigned max_subdivisions, double relative_tolerance,
    const std::function<void(unsigned)> &progress, const std::function<void()> &checkpoint)
{
    gcode.thermal_stress_mpa.clear();
    gcode.thermal_analysis = {};
    std::vector<DepositedSegment> segments;
    visit_deposited_moves(gcode, [&](std::size_t, const Eigen::Vector3f &start,
                                   const GCodeProcessorResult::MoveVertex &move) {
        if (move.extrusion_role == erMixed)
            throw std::invalid_argument("Mixed extrusion roles cannot be resolved for mechanical analysis.");
        segments.push_back({start.cast<double>() * .001,
            move.position.cast<double>() * .001, double(move.width) * .001,
            double(move.height) * .001, move.extruder_id});
    }, checkpoint);
    const auto grid = voxelize_deposition(segments, materials, spacing_m, bed_z_m, max_nodes, checkpoint);
    auto analysis = solve_refined_voxels(grid.voxels, grid.origin_m, grid.spacing_m, max_nodes,
        max_subdivisions, relative_tolerance, progress, checkpoint);
    if (!analysis.converged) return analysis;

    auto samples = sample_moves(gcode, analysis, checkpoint);
    if (checkpoint) checkpoint();
    gcode.thermal_stress_mpa = std::move(samples);
    return analysis;
}

} // namespace Slic3r::Thermoelastic
