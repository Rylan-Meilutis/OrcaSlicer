#include <catch2/catch_message.hpp>
#include "libslic3r/Format/ProjectStrategy.hpp"
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "libslic3r/Config.hpp"
#include "libslic3r/Format/bbs_3mf.hpp"
#include "libslic3r/GCode.hpp"
#include "libslic3r/GCode/DepositionProvenance.hpp"
#include "libslic3r/GCode/DepositionSection.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/Model.hpp"
#include "libslic3r/Preset.hpp"
#include "libslic3r/PresetBundle.hpp"
#include "libslic3r/Print.hpp"
#include "libslic3r/PrintConfig.hpp"
#include "libslic3r/Semver.hpp"
#include "libslic3r/ThermoelasticSourceDiagnostics.hpp"
#include "libslic3r/ThermoelasticMesh.hpp"
#include "libslic3r/ThermoelasticToolpath.hpp"
#include "libslic3r/ThermoelasticLayers.hpp"
#include "libslic3r/ThermoelasticLocalSections.hpp"
#include "libslic3r/ThermoelasticLocalTopology.hpp"
#include "libslic3r/ThermoelasticSectionPreview.hpp"
#include "libslic3r/ThermoelasticSections.hpp"
#include "libslic3r/ThermoelasticGCodeSections.hpp"
#include "libslic3r/ThermoelasticMaterial.hpp"
#include "libslic3r/ThermoelasticConflictOwnership.hpp"
#include "libslic3r/SectionGCodeHistory.hpp"
#include "libslic3r/SectionThermalAssembly.hpp"
#include "libslic3r/SectionIntegration.hpp"
#include "libslic3r/ThermalHistory.hpp"
#include "libslic3r/ThermalDeposition.hpp"
#include "libslic3r/ThermoelasticExactSections.hpp"
#include "libslic3r/ThermoelasticSectionFootprint.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Polygon.hpp"
#include "libslic3r/Point.hpp"
#include "libslic3r/Line.hpp"
#include "libslic3r/BoundingBox.hpp"
#include <Eigen/Core>
#include "libslic3r/ExtrusionEntity.hpp"
#include "libslic3r/Utils.hpp"
#include "libslic3r/libslic3r.h"
#include "test_utils.hpp"

#include <boost/filesystem/operations.hpp>
#include <boost/nowide/fstream.hpp>
#include <boost/nowide/cstdio.hpp>
#include <nlohmann/json.hpp>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <exception>
#include <iostream>
#include <initializer_list>
#include <ios>
#include <iosfwd>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include <utility>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <limits>
#include <sstream>
#include <stdexcept>

using namespace Slic3r;

namespace {
class ProjectSectionAudit : public GCode
{
public:
    using GCode::m_deposition_sections;
};

void audit_collapsed_envelopes(const GCodeProcessorResult& source, const Thermoelastic::TimedToolpath& path)
{
    std::size_t count = 0, unknown_volume = 0, neighbor_candidates = 0;
    double total_mm3 = 0., maximum_mm3 = 0., maximum_span_mm = 0.;
    for (const auto& point : path.unresolved_deposition) {
        if (point.reason != Thermoelastic::UnresolvedSpatialDeposition::Reason::CollapsedEnvelope) continue;
        ++count;
        double volume_mm3 = 0.;
        if (point.material < source.filament_diameters.size() &&
            std::isfinite(source.filament_diameters[point.material]) && source.filament_diameters[point.material] > 0.) {
            const double diameter = source.filament_diameters[point.material];
            volume_mm3 = std::acos(-1.) * .25 * diameter * diameter * point.extruded_filament_mm;
            total_mm3 += volume_mm3;
            maximum_mm3 = std::max(maximum_mm3, volume_mm3);
        } else ++unknown_volume;
        maximum_span_mm = std::max(maximum_span_mm, (point.end_m - point.start_m).norm() * 1000.);
        const auto next = std::lower_bound(path.deposition.begin(), path.deposition.end(), point.move_index,
            [](const auto& bead, std::size_t index) { return bead.move_index < index; });
        // Report only immediate same-source endpoint candidates. Their commanded
        // volume is not an omission allowance or a mechanical error bound.
        const auto candidate = [&](const auto& bead, bool before) {
            const auto& move = source.moves.at(bead.move_index);
            const bool touches = before ?
                bead.end_s == point.start_s && (bead.geometry.end_m - point.start_m).squaredNorm() == 0. :
                bead.start_s == point.end_s && (bead.geometry.start_m - point.end_m).squaredNorm() == 0.;
            if (!touches || move.deposition_source_id != point.source_id ||
                bead.geometry.material != point.material || !bead.commanded_area_m2) return 0.;
            ++neighbor_candidates;
            return *bead.commanded_area_m2 * (bead.geometry.end_m - bead.geometry.start_m).norm() * 1e9;
        };
        const double previous_mm3 = next == path.deposition.begin() ? 0. : candidate(*(next - 1), true);
        const double next_mm3 = next == path.deposition.end() ? 0. : candidate(*next, false);
        if (count <= 64)
            std::cout << "collapsed_envelope move=" << point.move_index << " source=" << point.source_id
                      << " command=" << point.command_id << " material=" << point.material
                      << " commanded_mm3=" << volume_mm3 << " span_mm=" << (point.end_m - point.start_m).norm() * 1000.
                      << " previous_candidate_mm3=" << previous_mm3 << " next_candidate_mm3=" << next_mm3 << '\n';
    }
    std::cout << "collapsed_envelope_summary count=" << count << " commanded_mm3=" << total_mm3
              << " maximum_commanded_mm3=" << maximum_mm3 << " maximum_span_mm=" << maximum_span_mm
              << " unknown_volume=" << unknown_volume << " neighbor_candidates=" << neighbor_candidates
              << "; no omission admitted; detail limited to64 spans\n" << std::flush;
    std::size_t vertical_count = 0, vertical_unknown_volume = 0;
    double vertical_filament_mm = 0., vertical_volume_mm3 = 0.;
    for (const auto& point : path.unresolved_deposition) {
        if (point.reason != Thermoelastic::UnresolvedSpatialDeposition::Reason::UnresolvedXYFootprint ||
            point.move_index >= source.moves.size() || !source.moves[point.move_index].deposition_vertical_extrusion)
            continue;
        ++vertical_count;
        vertical_filament_mm += point.extruded_filament_mm;
        std::optional<double> volume;
        if (point.material < source.filament_diameters.size() &&
            std::isfinite(source.filament_diameters[point.material]) && source.filament_diameters[point.material] > 0.) {
            const double diameter = source.filament_diameters[point.material];
            volume = std::acos(-1.) * .25 * diameter * diameter * point.extruded_filament_mm;
            vertical_volume_mm3 += *volume;
        } else ++vertical_unknown_volume;
        if (vertical_count <= 8)
            std::cout << "unresolved_pure_z move=" << point.move_index << " source=" << point.source_id
                      << " command=" << point.command_id << " material=" << point.material
                      << " filament_mm=" << point.extruded_filament_mm
                      << " commanded_mm3=" << (volume ? std::to_string(*volume) : "unknown")
                      << " start_z_mm=" << point.start_m.z() * 1000.
                      << " end_z_mm=" << point.end_m.z() * 1000.
                      << " start_s=" << point.start_s << " end_s=" << point.end_s << '\n';
    }
    std::cout << "unresolved_pure_z_summary count=" << vertical_count << " filament_mm=" << vertical_filament_mm
              << " commanded_mm3=" << vertical_volume_mm3 << " unknown_volume=" << vertical_unknown_volume
              << "; commanded E volume only, no spatial model or omission admitted; detail limited to8 spans\n" << std::flush;
    struct StationaryAudit {
        std::size_t count = 0, unknown_volume = 0;
        double filament_mm = 0., volume_mm3 = 0., maximum_mm3 = 0.;
    };
    std::array<StationaryAudit,2> stationary{};
    std::size_t details = 0;
    for (const auto& point : path.unresolved_deposition) {
        if (point.reason != Thermoelastic::UnresolvedSpatialDeposition::Reason::UnresolvedXYFootprint ||
            (point.end_m-point.start_m).squaredNorm()!=0.) continue;
        // Match the history's explicit role filter, without treating the
        // remaining entries as admitted or assigning them any spatial extent.
        const bool excluded = point.role==erWipeTower || point.role==erCustom || point.role==erNone;
        auto& summary = stationary[excluded ? 1 : 0];
        ++summary.count;
        const bool valid_e = std::isfinite(point.extruded_filament_mm) && point.extruded_filament_mm>0.;
        if (valid_e) summary.filament_mm += point.extruded_filament_mm;
        std::optional<double> volume;
        if (valid_e && point.material<source.filament_diameters.size() &&
            std::isfinite(source.filament_diameters[point.material]) && source.filament_diameters[point.material]>0.) {
            const double diameter = source.filament_diameters[point.material];
            const double candidate = std::acos(-1.)*.25*diameter*diameter*point.extruded_filament_mm;
            if (std::isfinite(candidate) && candidate>0.) volume = candidate;
        }
        if (volume) {
            summary.volume_mm3 += *volume;
            summary.maximum_mm3 = std::max(summary.maximum_mm3,*volume);
        } else ++summary.unknown_volume;
        if (details++<8) {
            std::cout << "unresolved_stationary move=" << point.move_index << " source=" << point.source_id
                      << " command=" << point.command_id << " material=" << point.material
                      << " role=" << unsigned(point.role) << " excluded_role=" << excluded
                      << " filament_mm=" << point.extruded_filament_mm << " commanded_mm3=";
            if (volume) std::cout << *volume;
            else std::cout << "unknown";
            std::cout << " start_s=" << point.start_s << " end_s=" << point.end_s << '\n';
        }
    }
    for (std::size_t group=0;group<stationary.size();++group) {
        const auto& summary = stationary[group];
        std::cout << "unresolved_stationary_summary excluded_role=" << group << " count=" << summary.count
                  << " filament_mm=" << summary.filament_mm << " commanded_mm3=" << summary.volume_mm3
                  << " commanded_m3=" << summary.volume_mm3*1e-9 << " maximum_commanded_mm3=" << summary.maximum_mm3
                  << " unknown_volume=" << summary.unknown_volume
                  << "; zero3D displacement only; no omission admitted; detail limited to8 spans total\n";
    }
    std::cout << std::flush;
}

Thermoelastic::SectionBedReleaseRegion mechanical_release_region(const std::string& text)
{
    if (text.empty() || text.size() > 256 || text.back() == ',')
        throw std::invalid_argument("Expected a bounded release rectangle in millimetres.");
    std::istringstream input(text);
    std::array<double,4> values;
    std::string token;
    std::size_t count = 0;
    while (std::getline(input,token,',')) {
        if (count == values.size()) throw std::invalid_argument("Too many release coordinates.");
        std::size_t consumed = 0;
        const double value = std::stod(token,&consumed)*.001;
        if (consumed != token.size() || !std::isfinite(value) || std::abs(value) > 10.)
            throw std::invalid_argument("Invalid release coordinate.");
        values[count++] = value;
    }
    if (count != 4 || values[0] >= values[2] || values[1] >= values[3])
        throw std::invalid_argument("Empty or incomplete release rectangle.");
    return {{values[0],values[1]},{values[2],values[3]}};
}
std::vector<double> mechanical_z_planes(const std::string& text)
{
    if (text.empty() || text.size() > 4096 || text.back() == ',')
        throw std::invalid_argument("Invalid graded plane list.");
    std::istringstream input(text);
    std::vector<double> planes;
    std::string token;
    while (std::getline(input,token,',')) {
        std::size_t consumed = 0;
        const double plane = std::stod(token,&consumed) * .001;
        if (consumed != token.size() || !std::isfinite(plane) || std::abs(plane) > 10. ||
            planes.size() == 128 || (!planes.empty() && plane <= planes.back()))
            throw std::invalid_argument("Invalid graded plane coordinate.");
        planes.push_back(plane);
    }
    if (planes.size() < 2) throw std::invalid_argument("Graded planes must enclose a solid.");
    return planes;
}
std::vector<Eigen::Vector3d> mechanical_bases(const Eigen::Vector3d& base, const std::string& mode)
{
    if (mode != "1" && mode != "x" && mode != "y" && mode != "z" && mode != "xy")
        throw std::invalid_argument("Expected independent axes (1), x, y, z or combined lateral (xy) basis probe.");
    std::vector<Eigen::Vector3d> result(4, base);
    if (mode == "xy") {
        result[1].x() *= .5;
        result[2].y() *= .5;
        result[3].head<2>() *= .5;
        return result;
    }
    for (unsigned run = 1; run < result.size(); ++run) {
        const unsigned axis = mode == "1" ? run-1 : mode == "x" ? 0 : mode == "y" ? 1 : 2;
        result[run][axis] *= mode == "1" ? .5 : std::ldexp(1., -int(run));
    }
    return result;
}
std::vector<double> mechanical_quanta(const std::string& text)
{
    if (text.empty() || text.size() > 128) throw std::invalid_argument("Invalid mechanical quantum list.");
    std::vector<double> result;
    std::size_t begin = 0;
    long long previous_nm = 0;
    for (;;) {
        const auto comma = text.find(',', begin);
        const auto value = text.substr(begin, comma == std::string::npos ? comma : comma - begin);
        if (value.empty() || result.size() == 4) throw std::invalid_argument("Expected two to four mechanical quanta.");
        std::size_t consumed = 0;
        const double quantum = std::stod(value, &consumed) * 1e-6;
        if (consumed != value.size() || !std::isfinite(quantum) || quantum < 2e-9 || quantum > .001)
            throw std::invalid_argument("Invalid mechanical quantum.");
        const auto native = std::llround(quantum * 1e9);
        if (std::abs(quantum * 1e9 - double(native)) > 1e-6 || (previous_nm && 2 * native != previous_nm))
            throw std::invalid_argument("Mechanical quanta must halve on the native grid.");
        result.push_back(quantum);
        previous_nm = native;
        if (comma == std::string::npos) break;
        begin = comma + 1;
    }
    if (result.size() < 2) throw std::invalid_argument("Expected two to four mechanical quanta.");
    return result;
}
}

TEST_CASE("Captured native ownership failures replay without project history", "[.][ThermoelasticOwnershipReplay]")
{
    using namespace Thermoelastic;
    const char* path=std::getenv("ORCA_FIXED_HISTORY_REPLAY_FIXTURE");
    REQUIRE(path!=nullptr);
    boost::nowide::ifstream input(path,std::ios::binary);
    REQUIRE(input.good());
    input.seekg(0,std::ios::end);
    const auto bytes=input.tellg();
    REQUIRE(bytes>std::streampos(0));
    REQUIRE(bytes<=std::streampos(1024*1024));
    input.seekg(0);
    const auto fixture=nlohmann::json::parse(input);
    REQUIRE(fixture.at("schema").get<unsigned>()==1);
    REQUIRE(fixture.at("native_fill").get<std::string>()=="nonzero-winding");
    REQUIRE(fixture.at("allowance_nm").get<unsigned>()==7);
    const auto max_work=fixture.at("max_work").get<std::size_t>();
    REQUIRE(max_work>0);
    REQUIRE(max_work<=4000000);
    const auto choice=fixture.at("choice").get<std::string>();
    const bool supported_choice=choice=="first" || choice=="second";
    REQUIRE(supported_choice);
    std::size_t vertices=0;
    const auto ring=[&](const nlohmann::json& points) {
        REQUIRE(points.is_array());
        REQUIRE(points.size()>=3);
        REQUIRE(points.size()<=2048);
        Polygon polygon;
        for (const auto& point:points) {
            REQUIRE(point.is_array());
            REQUIRE(point.size()==2);
            REQUIRE(point[0].is_number_integer());
            REQUIRE(point[1].is_number_integer());
            const auto x=point[0].get<coord_t>(),y=point[1].get<coord_t>();
            REQUIRE(std::abs(double(x))<=1e10);
            REQUIRE(std::abs(double(y))<=1e10);
            REQUIRE(++vertices<=8192);
            polygon.points.emplace_back(x,y);
        }
        return polygon;
    };
    const auto polygons=[&](const nlohmann::json& operand) {
        REQUIRE(operand.is_array());
        REQUIRE(operand.size()<=2048);
        const auto begin=vertices;
        ExPolygons output;
        for (const auto& value:operand) {
            ExPolygon polygon(ring(value.at("contour")));
            const auto& holes=value.at("holes");
            REQUIRE(holes.is_array());
            REQUIRE(holes.size()<=2048);
            for (const auto& hole:holes) polygon.holes.push_back(ring(hole));
            output.push_back(std::move(polygon));
        }
        REQUIRE(vertices-begin<=2048);
        return output;
    };
    REQUIRE(fixture.at("native").size()==2);
    REQUIRE(fixture.at("envelopes").size()==2);
    const std::array<ExPolygons,2> masks{polygons(fixture.at("native")[0]),polygons(fixture.at("native")[1])};
    const std::array<ExPolygons,2> envelopes{polygons(fixture.at("envelopes")[0]),polygons(fixture.at("envelopes")[1])};
    const bool single_sided=fixture.at("single_sided").get<bool>();
    const bool recorded_dual=fixture.value("dual_boundary",false);
    if (const char* override_dual=std::getenv("ORCA_FIXED_HISTORY_REPLAY_DUAL_BOUNDARY")) {
        REQUIRE(std::string(override_dual)=="1");
        REQUIRE(fixture.at("exact").get<bool>());
        const auto first=make_exact_section(envelopes[0]).section,second=make_exact_section(envelopes[1]).section;
        const auto deposited=boolean_exact_sections(first,second,ExactSectionOperation::Union).section;
        for (const auto priority:{ExactInterfaceChoice::FavorFirst,ExactInterfaceChoice::FavorSecond}) {
            const auto partition=exact_partition_uncertain_interface(masks[0],masks[1],envelopes,priority,
                7,max_work,{},NativeWitnessFill::NonzeroWinding,single_sided,true);
            CHECK(boolean_exact_sections(partition.materials[0],partition.materials[1],ExactSectionOperation::Intersection).section.empty());
            const auto occupied=boolean_exact_sections(partition.materials[0],partition.materials[1],ExactSectionOperation::Union).section;
            CHECK(boolean_exact_sections(occupied,deposited,ExactSectionOperation::SymmetricDifference).section.empty());
            std::cout<<"Conditional dual-boundary sensitivity, priority="<<(priority==ExactInterfaceChoice::FavorFirst ? "first" : "second")
                <<", ambiguity_area_m2="<<partition.dual_native_boundary_area_m2<<", work="<<partition.work
                <<"; unchanged deposited union, NOT a physical uncertainty bound\n"<<std::flush;
        }
        return;
    }
    try {
        if (fixture.at("exact").get<bool>())
            exact_partition_uncertain_interface(masks[0],masks[1],envelopes,
                choice=="first" ? ExactInterfaceChoice::FavorFirst : ExactInterfaceChoice::FavorSecond,
                7,max_work,{},NativeWitnessFill::NonzeroWinding,single_sided,recorded_dual);
        else
            partition_uncertain_interface({}, {},masks[0],masks[1],7,max_work,{},&envelopes,
                NativeWitnessFill::NonzeroWinding,single_sided);
        FAIL("Captured rejected ownership unexpectedly passed its original proof.");
    } catch (const std::invalid_argument& error) {
        const auto recorded=fixture.at("failure").get<std::string>();
        const auto category=recorded.substr(0,recorded.find(':'));
        CHECK(std::string(error.what()).find(category)==0);
        std::cout<<"Replayed ownership failure: "<<error.what()<<'\n'<<std::flush;
    }
}

TEST_CASE("Shared-prefix mechanical basis probes isolate the requested axis", "[DepositionProvenance]")
{
    const Eigen::Vector3d base(.02,.02,.005);
    const auto independent = mechanical_bases(base,"1");
    CHECK_THAT((independent[2]-Eigen::Vector3d(.02,.01,.005)).norm(),Catch::Matchers::WithinAbs(0.,1e-15));
    const auto y = mechanical_bases(base,"y");
    REQUIRE(y.size() == 4);
    CHECK_THAT((y[3]-Eigen::Vector3d(.02,.0025,.005)).norm(),Catch::Matchers::WithinAbs(0.,1e-15));
    const auto xy = mechanical_bases(base,"xy");
    REQUIRE(xy.size() == 4);
    CHECK_THAT((xy[1]-Eigen::Vector3d(.01,.02,.005)).norm(),Catch::Matchers::WithinAbs(0.,1e-15));
    CHECK_THAT((xy[2]-Eigen::Vector3d(.02,.01,.005)).norm(),Catch::Matchers::WithinAbs(0.,1e-15));
    CHECK_THAT((xy[3]-Eigen::Vector3d(.01,.01,.005)).norm(),Catch::Matchers::WithinAbs(0.,1e-15));
    CHECK_THROWS_AS(mechanical_bases(base,"all"),std::invalid_argument);
}

TEST_CASE("Diagnostic graded plane lists are bounded and increasing", "[DepositionProvenance]")
{
    const auto planes = mechanical_z_planes("0,0.625,1.25,2.5");
    REQUIRE(planes.size() == 4);
    CHECK_THAT(planes.back(),Catch::Matchers::WithinAbs(.0025,1e-15));
    for (const char* invalid : {"","0","0,","0,,1","0,1,1","0,2,1","0,nan","0,2junk"})
        CHECK_THROWS_AS(mechanical_z_planes(invalid),std::invalid_argument);
}

TEST_CASE("Diagnostic bed release rectangles reject incomplete or nonfinite bounds", "[DepositionProvenance]")
{
    const auto region = mechanical_release_region("0,150,80,210");
    CHECK_THAT(region.minimum_m.y(),Catch::Matchers::WithinAbs(.15,1e-15));
    CHECK_THAT(region.maximum_m.x(),Catch::Matchers::WithinAbs(.08,1e-15));
    for (const char* invalid : {"","0,1,2","0,1,2,3,4","0,1,2,3,","0,nan,2,3","0,0,0,1","0,0,1,2junk"})
        CHECK_THROWS_AS(mechanical_release_region(invalid),std::invalid_argument);
}

TEST_CASE("Shared-prefix mechanical diagnostic quanta are bounded nested halvings", "[DepositionProvenance]")
{
    const auto quanta = mechanical_quanta("10,5,2.5");
    REQUIRE(quanta.size() == 3);
    CHECK_THAT(quanta[0],Catch::Matchers::WithinAbs(1e-5,1e-18));
    CHECK_THAT(quanta[2],Catch::Matchers::WithinAbs(2.5e-6,1e-18));
    for (const char* invalid : {"", "10", "10,", "10,6", "10,5,2.5,1.25,0.625", "10junk,5", "nan,5", "0,0", "1001,500.5"})
        CHECK_THROWS_AS(mechanical_quanta(invalid),std::invalid_argument);
}

TEST_CASE("Saved project exports audit local deposition provenance without changing settings",
          "[DepositionProvenance][.DepositionProvenanceProject]")
{
    const char* path = std::getenv("ORCA_PROVENANCE_PROJECT");
    if (!path)
        SKIP("Set ORCA_PROVENANCE_PROJECT to a saved single-plate 3MF.");
    const auto previous_resources = resources_dir();
    set_resources_dir(std::string(TEST_DATA_DIR) + "/../../resources");
    const ScopeGuard restore_resources([&] { set_resources_dir(previous_resources); });
    ScopedSlic3rTemporaryDir temporary_dir("deposition_provenance_project");
    constexpr std::uintmax_t reserve = std::uintmax_t(30) * 1024 * 1024 * 1024;
    REQUIRE(boost::filesystem::space(temporary_dir.path()).available >= reserve);

    DynamicPrintConfig imported;
    ConfigSubstitutionContext substitutions{ForwardCompatibilitySubstitutionRule::Enable};
    PlateDataPtrs plates;
    std::vector<Preset*> presets;
    const ScopeGuard release_import([&] {
        release_PlateData_list(plates);
        for (Preset* preset : presets)
            delete preset;
    });
    bool is_bbl = false;
    Semver version;
    auto model = Model::read_from_file(path, &imported, &substitutions,
                                       LoadStrategy::LoadModel | LoadStrategy::LoadConfig | LoadStrategy::AddDefaultInstances, &plates,
                                       &presets, &is_bbl, &version);
    REQUIRE_FALSE(model.objects.empty());
    REQUIRE(plates.size() <= 1);
    Preset::normalize(imported);
    PresetBundle bundle;
    bundle.load_config_model(path, imported, version);
    auto config = bundle.full_config(false);
    if (!plates.empty())
        config.apply(plates.front()->config);
    // Read-only isolation experiments: only this temporary in-memory copy is
    // changed. Never save these diagnostic settings into the user's project.
    const bool disable_beams = std::getenv("ORCA_PROVENANCE_DISABLE_BEAMS") != nullptr;
    const bool disable_roots = std::getenv("ORCA_PROVENANCE_DISABLE_ROOTS") != nullptr;
    const char* standard_walls = std::getenv("ORCA_PROVENANCE_STANDARD_WALLS");
    const char* planar_surfaces = std::getenv("ORCA_PROVENANCE_PLANAR_SURFACES");
    if (standard_walls) {
        REQUIRE(std::string(standard_walls) == "1");
        REQUIRE_FALSE(disable_beams);
        REQUIRE_FALSE(disable_roots);
        DynamicPrintConfig candidate;
        candidate.set_deserialize_strict("perimeter_layering", "standard");
        config.apply(candidate);
        for (auto* object : model.objects) {
            object->config.apply(candidate);
            for (auto* volume : object->volumes) volume->config.apply(candidate);
        }
        std::cout << "Explicit numerical candidate ONLY: perimeter_layering=standard; temporary config/object/volume overrides; "
                     "all other settings, roots, beams and model geometry retained; saved project unchanged; "
                     "NOT validation of the original wall geometry and NOT automatic adoption\n" << std::flush;
    }
    if (planar_surfaces) {
        REQUIRE(std::string(planar_surfaces)=="1");
        REQUIRE(standard_walls!=nullptr);
        DynamicPrintConfig candidate;
        candidate.set_deserialize_strict("top_surface_z_mode","disabled");
        candidate.set_deserialize_strict("zaa_enabled","0");
        candidate.set_deserialize_strict("nonplanar_top_surface","0");
        config.apply(candidate);
        for (auto* object:model.objects) {
            object->config.apply(candidate);
            for (auto* volume:object->volumes) volume->config.apply(candidate);
        }
        std::cout << "Additional numerical candidate ONLY: top-surface Z methods disabled in temporary config/object/volume copies; "
                     "roots, beams and model mesh retained; saved project unchanged; "
                     "NOT validation of the original nonplanar toolpaths and NOT automatic adoption\n" << std::flush;
    }
    for (const auto& key : {"interlocking_beam", "rooting", "rooting_xy"}) {
        if ((std::string(key) == "interlocking_beam") ? !disable_beams : !disable_roots) continue;
        config.set_deserialize_strict(key, "0");
        for (auto* object : model.objects) object->config.set(key, false);
    }
    if (disable_beams || disable_roots)
        std::cout << "Diagnostic configuration ONLY: beams_disabled=" << disable_beams
                  << " roots_disabled=" << disable_roots << "; saved project unchanged\n";
    Print print;
    print.set_status_silent();
    print.apply(model, config);
    for (const auto* object : print.objects()) {
        if (disable_beams) REQUIRE_FALSE(object->config().interlocking_beam.value);
        if (disable_roots) {
            REQUIRE_FALSE(object->config().rooting.value);
            REQUIRE_FALSE(object->config().rooting_xy.value);
        }
    }
    // This audit uses the Cartesian generator, not the belt subclass.
    REQUIRE_FALSE(print.config().belt_printer.value);
    print.process();
    if (standard_walls) {
        REQUIRE(print.num_print_regions() > 0);
        for (std::size_t i = 0; i < print.num_print_regions(); ++i)
            REQUIRE(print.get_print_region(i).config().perimeter_layering.value == PerimeterLayeringMode::Standard);
    }
    if (planar_surfaces)
        for (std::size_t i=0;i<print.num_print_regions();++i) {
            const auto& region=print.get_print_region(i).config();
            REQUIRE(region.top_surface_z_mode.value==TopSurfaceZMode::Disabled);
            REQUIRE_FALSE(region.zaa_enabled.value);
            REQUIRE_FALSE(region.nonplanar_top_surface.value);
        }
    REQUIRE(boost::filesystem::space(temporary_dir.path()).available >= reserve);

    ScopedTemporaryFile output(".gcode");
    ProjectSectionAudit generator;
    print.update_imex_slice_offset();
    const auto origin = print.get_plate_origin();
    const auto imex = print.get_imex_slice_offset();
    generator.set_gcode_offset_with_imex_shift(origin.x(),origin.y(),imex.x(),imex.y());
    generator.enable_deposition_provenance();
    GCodeProcessorResult result;
    generator.do_export(&print, output.string().c_str(), &result);
    if (std::getenv("ORCA_PROVENANCE_PARENT_COMMAND_AUDIT") && generator.m_deposition_sections) {
        // Bounded diagnostic for the first missing-child failure. This uses
        // the generator registry even when final verification rejects it.
        const std::string prefix = "Missing/reordered section command: expected ";
        if (result.deposition_sections_failure.compare(0, prefix.size(), prefix) == 0) {
            const auto first = std::stoul(result.deposition_sections_failure.substr(prefix.size()));
            for (std::uint32_t id = first > 1 ? first - 1 : 1; id <= first + 3; ++id) {
                const auto* command = generator.m_deposition_sections->section(id);
                if (!command) break;
                std::cout.precision(17);
                std::cout << "parent_command_audit id=" << id << " transform=" << unsigned(command->transform)
                    << " start=" << command->start[0] << ',' << command->start[1] << ',' << command->start[2]
                    << " end=" << command->end[0] << ',' << command->end[1] << ',' << command->end[2]
                    << " height=" << command->start_height << ',' << command->end_height
                    << " width=" << command->width << " filament=" << command->filament << '\n';
            }
        }
    }
    const bool gated = print.config().spiral_mode.value || !print.config().time_lapse_gcode.value.empty();
    if (gated)
        CHECK_FALSE(result.deposition_provenance);
    else
        REQUIRE(result.deposition_provenance);

    // Four counts per final logical filament: unowned, unknown-height,
    // native-planar, transformed. None of these qualify physical ownership.
    std::map<unsigned, std::array<std::size_t, 4>> counts;
    std::size_t extrusions = 0;
    for (const auto& move : result.moves) {
        if (move.type != EMoveType::Extrude)
            continue;
        ++extrusions;
        auto& row = counts[move.extruder_id];
        if (move.deposition_source_id == 0) {
            ++row[0];
            continue;
        }
        REQUIRE(result.deposition_provenance);
        const auto* source = result.deposition_provenance->source(move.deposition_source_id);
        REQUIRE(source != nullptr);
        REQUIRE(source->logical_filament == move.extruder_id);
        const auto column = source->height == DepositionSource::Height::NativePlanar ? 2 :
                            source->height == DepositionSource::Height::Transformed  ? 3 :
                                                                                       1;
        ++row[column];
    }
    REQUIRE(extrusions > 0);

    boost::nowide::ifstream file(output.string(), std::ios::binary);
    REQUIRE(file.is_open());
    // Fixed-size streaming search includes chunk-boundary matches without
    // retaining the G-code or relying on its maximum line length.
    constexpr const char* marker = "_ORCA_INTERNAL_SOURCE";
    std::array<char, 65536> buffer{};
    std::string tail;
    bool leaked_marker = false;
    while (file.read(buffer.data(), buffer.size()) || file.gcount() != 0) {
        tail.append(buffer.data(), static_cast<std::size_t>(file.gcount()));
        leaked_marker |= tail.find(marker) != std::string::npos;
        leaked_marker |= tail.find("_ORCA_INTERNAL_SECTION") != std::string::npos;
        if (tail.size() > 32)
            tail.erase(0, tail.size() - 32);
    }
    CHECK(file.eof());
    CHECK_FALSE(leaked_marker);
    std::cout << "Provenance project audit: bytes=" << boost::filesystem::file_size(output.path()) << " extrusion_moves=" << extrusions
              << " registry=" << bool(result.deposition_provenance) << " ownership_snapshot=" << bool(result.native_ownership)
              << " unsupported_gate=" << gated << '\n';
    std::size_t section_moves = 0, variable_sections = 0;
    for (const auto& move : result.moves) {
        if (move.internal_only || !move.deposition_design_heights) continue;
        REQUIRE(result.deposition_sections);
        REQUIRE(move.deposition_design_width);
        ++section_moves;
        variable_sections += (*move.deposition_design_heights)[0] != (*move.deposition_design_heights)[1];
    }
    std::cout << "command_sections registry=" << bool(result.deposition_sections)
              << " failure=" << result.deposition_sections_failure
              << " parent_commands=" << (result.deposition_sections ? result.deposition_sections->size() : 0)
              << " verified_moves=" << section_moves << " variable_moves=" << variable_sections << '\n';
    for (const auto& [filament, row] : counts)
        std::cout << "logical_filament=" << filament << " unowned=" << row[0] << " unknown_height=" << row[1] << " native_planar=" << row[2]
                  << " transformed=" << row[3] << '\n';
    if (std::getenv("ORCA_PROVENANCE_MISSING_COMMANDS")) {
        REQUIRE(result.deposition_provenance);
        std::vector<unsigned int> lines;
        std::size_t missing = 0;
        for (std::size_t i = 0; i < result.moves.size(); ++i) {
            const auto& move = result.moves[i];
            if (move.internal_only || move.type != EMoveType::Extrude || move.deposition_command_section_id ||
                (move.extrusion_role != erPerimeter && move.extrusion_role != erExternalPerimeter))
                continue;
            const auto* source = result.deposition_provenance->source(move.deposition_source_id);
            if (!source || source->height != DepositionSource::Height::Transformed)
                continue;
            ++missing;
            if (lines.size() < 8) {
                lines.push_back(move.gcode_id);
                std::cout << "Missing transformed command: move=" << i << " source=" << move.deposition_source_id
                          << " line=" << move.gcode_id << " width=" << move.width << " height=" << move.height << '\n';
            }
        }
        std::cout << "Missing transformed perimeter commands=" << missing << '\n';
        boost::nowide::ifstream commands(output.string());
        std::string line;
        unsigned int number = 0;
        while (std::getline(commands,line)) {
            ++number;
            if (std::any_of(lines.begin(),lines.end(),[&](unsigned int id) {
                    return number >= id && number - id <= 1;
                }))
                std::cout << "Emitted line " << number << ": " << line.substr(0,512) << '\n';
        }
        return;
    }
    if (std::getenv("ORCA_FIXED_HISTORY_ADMISSION")) {
        using namespace Thermoelastic;
        // Full-project numerical admission only. Explicit uniform/isothermal
        // synthetic properties, not ASA/TPU calibration or an adhesion claim.
        const auto timeline = timed_deposition(result, 2000000, 2000000);
        audit_collapsed_envelopes(result, timeline);
        REQUIRE(timeline.state == ToolpathState::Complete);
        REQUIRE_FALSE(counts.empty());
        const std::size_t material_count = std::size_t(counts.rbegin()->first) + 1;
        REQUIRE(material_count <= 256);
        SectionGCodeScenario scenario(BeadSection::RoundedRectangle);
        scenario.acknowledge_modeled_timing    = true;
        scenario.use_native_ownership          = true;
        scenario.use_fixed_interface_reference = true;
        scenario.spacing_m                     = {.005, .005, .001};
        // Independent basis refinement must not change the deposition envelope
        // or the tile partition used to integrate it.
        const std::array<const char*,3> fixed_basis_variables{
            "ORCA_FIXED_HISTORY_BASIS_X_MM", "ORCA_FIXED_HISTORY_BASIS_Y_MM", "ORCA_FIXED_HISTORY_BASIS_Z_MM"};
        for (unsigned axis=0;axis<fixed_basis_variables.size();++axis) {
            if (const char* setting=std::getenv(fixed_basis_variables[axis])) {
                std::size_t consumed=0;
                const std::string value(setting);
                const double millimetres=std::stod(value,&consumed);
                REQUIRE(consumed==value.size());
                REQUIRE(std::isfinite(millimetres));
                REQUIRE(millimetres>=.2);
                REQUIRE(millimetres<=20.);
                scenario.spacing_m[axis]=millimetres*.001;
            }
        }
        std::cout << "Fixed history basis_m=" << scenario.spacing_m.transpose() << '\n';
        scenario.stream_tile_size_m            = Eigen::Vector2d{.005, .005};
        if (const char* tile_setting = std::getenv("ORCA_FIXED_HISTORY_TILE_MM")) {
            std::size_t consumed = 0;
            const std::string value(tile_setting);
            const double tile_mm = std::stod(value, &consumed);
            REQUIRE(consumed == value.size());
            REQUIRE(std::isfinite(tile_mm));
            REQUIRE(tile_mm >= 1.);
            REQUIRE(tile_mm <= 10.);
            scenario.stream_tile_size_m = Eigen::Vector2d::Constant(tile_mm * .001);
        }
        scenario.stream_section_height_m       = .0002;
        scenario.stream_vertical_normal_bins   = 32;
        if (const char* bins=std::getenv("ORCA_FIXED_HISTORY_NORMAL_BINS"))
            scenario.stream_vertical_normal_bins=unsigned(std::stoul(bins));
        std::size_t activation_count = 1;
        if (const char* setting=std::getenv("ORCA_FIXED_HISTORY_ACTIVATIONS")) {
            const std::string value(setting);
            REQUIRE(value.size()==1);
            REQUIRE(value.front()>='1');
            REQUIRE(value.front()<='8');
            activation_count=std::size_t(value.front()-'0');
        }
        REQUIRE(std::isfinite(timeline.total_s));
        REQUIRE(timeline.total_s>0.);
        // Independent temporal batching control, not continuous deposition or
        // physical calibration. Preserve the exact final processor timestamp.
        for (std::size_t i=1;i<activation_count;++i)
            scenario.activation_times_s.push_back(timeline.total_s*(double(i)/double(activation_count)));
        scenario.activation_times_s.push_back(timeline.total_s);
        for (std::size_t i = 0; i < material_count; ++i) {
            scenario.filament_to_material.push_back(i);
            scenario.thermal_materials.push_back({Eigen::Matrix3d::Identity() * .2, 2e6, 400.});
            scenario.mechanical_materials.push_back({Eigen::Matrix<double, 6, 6>::Identity() * 2e9});
            scenario.expansion.push_back({{{250., 600., Eigen::Matrix3d::Identity() * 1e-5}}, Eigen::Matrix3d::Identity(), 400.});
        }
        const bool homogeneous = std::getenv("ORCA_FIXED_HISTORY_HOMOGENEOUS") != nullptr;
        if (homogeneous) {
            REQUIRE(std::string(std::getenv("ORCA_FIXED_HISTORY_HOMOGENEOUS")) == "1");
            // A synthetic homogeneous control isolates whole-model streaming
            // and solver admission from native multimaterial interfaces.
            // It is NOT the physical ASA/TPU case or a bound on its curl.
            scenario.filament_to_material.assign(material_count, 0);
            scenario.use_native_ownership = false;
        }
        scenario.deposition_temperatures_k = [material_count](double) { return std::vector<double>(material_count, 400.); };
        scenario.exposure                  = [](double, const SectionThermalAssembly&) { return std::vector<ThermalHistory::Element>{}; };
        const bool cooling_control=std::getenv("ORCA_FIXED_HISTORY_COOLING_CONTROL") != nullptr;
        std::vector<SectionThermalNode> cooling_nodes;
        if (cooling_control) {
            // Synthetic coupled-load control, not inferred printer fan physics.
            scenario.cooldown_s=10.;
            scenario.exposure=[&cooling_nodes](double,const SectionThermalAssembly& assembly) {
                cooling_nodes=assembly.nodes;
                for (std::size_t i=0;i<assembly.initial.deposited_cells.size();++i) {
                    const auto& cell=assembly.initial.deposited_cells[i];
                    if (!(cell.volume_m3>0.) || !(cell.volumetric_heat_capacity_j_m3_k>0.) || !(cell.temperature_k>0.))
                        throw std::invalid_argument("Synthetic thermal cell invalid: id="+std::to_string(i)+
                            ", volume="+std::to_string(cell.volume_m3)+", temperature="+std::to_string(cell.temperature_k));
                }
                for (std::size_t i=0;i<assembly.initial.added_elements.size();++i) {
                    const auto& element=assembly.initial.added_elements[i];
                    if (!ThermalHistory::valid_element(element,assembly.nodes.size())) {
                        std::ostringstream diagnostic;
                        diagnostic.precision(std::numeric_limits<double>::max_digits10);
                        diagnostic << "Synthetic conduction control invalid: id=" << i << ", matrix=";
                        for (unsigned row=0;row<8;++row)
                            for (unsigned column=0;column<8;++column) diagnostic << element.conductance(row,column) << ',';
                        throw std::invalid_argument(diagnostic.str());
                    }
                }
                std::vector<ThermalHistory::Element> boundaries;
                const auto add=[&](const auto& face,bool bed) {
                    const auto robin=SectionIntegration::thermal_boundary_element(face.area_moments,bed ? 100. : 20.,bed ? 400. : 300.,
                        assembly.quadrature);
                    ThermalHistory::Element element;
                    element.nodes=face.nodes;
                    element.conductance=robin.conductance;
                    element.load=robin.load;
                    element.boundary=true;
                    element.kind=bed ? ThermalHistory::BoundaryKind::Bed : ThermalHistory::BoundaryKind::Ambient;
                    if (!ThermalHistory::valid_element(element,assembly.nodes.size())) {
                        std::ostringstream diagnostic;
                        diagnostic.precision(std::numeric_limits<double>::max_digits10);
                        diagnostic << "Synthetic boundary control invalid: bed=" << bed << ", moments=";
                        for (const auto value : face.area_moments) diagnostic << value << ',';
                        diagnostic << " nodes=";
                        for (const auto node : face.nodes) diagnostic << node << ',';
                        throw std::invalid_argument(diagnostic.str());
                    }
                    boundaries.push_back(std::move(element));
                };
                for (const auto& face : assembly.horizontal_faces) add(face,face.kind==SectionHorizontalFaceKind::Bed);
                for (const auto& face : assembly.vertical_faces) add(face,false);
                if (std::getenv("ORCA_FIXED_HISTORY_THERMAL_STUDY")) {
                    // Reuse the identical assembled geometry for a temporal
                    // sensitivity study; no disk cache or modified boundaries.
                    const auto study_started=std::chrono::steady_clock::now();
                    for (const unsigned quadrature : {0u,1u,2u}) for (const double dt : {1.,.5,.25}) {
                        auto event=assembly.initial;
                        event.time_s=0.;
                        event.exposed_elements=boundaries;
                        // Diagnostic quadrature comparison only. Row-sum
                        // boundary lumping preserves reservoir loads and heat
                        // transfer for a constant field; it is not enabled in
                        // the actual scenario or accepted as a curl model.
                        if (quadrature) for (auto& element : event.exposed_elements) {
                            const auto diagonal=element.conductance.rowwise().sum().eval();
                            element.conductance=diagonal.asDiagonal();
                        }
                        if (quadrature==2) for (auto& element : event.added_elements) {
                            // Synthetic diagonal-conductivity experiment only:
                            // transverse gradient lumping yields twelve positive
                            // edge conductances. Check preservation of the exact
                            // affine patch action; mixed modes change and still
                            // require independent spatial refinement.
                            const auto original=element.conductance;
                            element.conductance.setZero();
                            for (unsigned axis=0;axis<3;++axis) {
                                Eigen::Matrix<double,8,1> coordinate;
                                for (unsigned node=0;node<8;++node) coordinate[node]=(node>>axis)&1u;
                                const auto force=(original*coordinate).eval();
                                for (unsigned low=0;low<8;++low) if (!(low&(1u<<axis))) {
                                    const unsigned high=low|(1u<<axis);
                                    const double conductance=.5*(force[high]-force[low]);
                                    if (!std::isfinite(conductance) || conductance<0.)
                                        throw std::invalid_argument("Synthetic directional quadrature has a negative conductance.");
                                    element.conductance(low,low)+=conductance;
                                    element.conductance(high,high)+=conductance;
                                    element.conductance(low,high)-=conductance;
                                    element.conductance(high,low)-=conductance;
                                }
                            }
                            for (unsigned axis=0;axis<3;++axis) {
                                Eigen::Matrix<double,8,1> coordinate;
                                for (unsigned node=0;node<8;++node) coordinate[node]=(node>>axis)&1u;
                                if (((element.conductance-original)*coordinate).norm()>1e-10*original.norm())
                                    throw std::invalid_argument("Synthetic directional quadrature fails the affine patch check.");
                            }
                        }
                        std::vector<ThermalHistory::DepositionEvent> events;
                        events.push_back(std::move(event));
                        ThermalHistory::HistoryLimits trial_limits;
                        trial_limits.max_step_s=dt;
                        std::cout << "Synthetic thermal time-step study dt_s=" << dt
                                  << " quadrature=" << quadrature << '\n' << std::flush;
                        const auto trial=ThermalHistory::evolve(events,10.,trial_limits,[&] {
                            if (std::chrono::steady_clock::now()-study_started>std::chrono::seconds(60))
                                throw ResourceLimit("Synthetic thermal time-step study exceeded sixty seconds.");
                        });
                        std::cout << "Synthetic thermal study state=" << int(trial.state)
                                  << " accepted_steps=" << trial.steps;
                        if (!trial.temperatures_k.empty()) {
                            const auto range=std::minmax_element(trial.temperatures_k.begin(),trial.temperatures_k.end());
                            std::cout << " minimum_k=" << *range.first << " maximum_k=" << *range.second
                                      << " balance_error_j=" << trial.energy_balance_error_j;
                        }
                        if (trial.rejected_temperature) {
                            const auto& failed=*trial.rejected_temperature;
                            std::cout << " rejected_kelvin=" << failed.kelvin << " node=" << failed.node;
                            if (failed.node<assembly.nodes.size())
                                std::cout << " location_m=" << assembly.nodes[failed.node].position_m.transpose()
                                          << " component=" << assembly.nodes[failed.node].component;
                        }
                        std::cout << '\n' << std::flush;
                    }
                }
                return boundaries;
            };
        }
        scenario.stationary_allowance      = StationaryExtrusionAllowance{1e-12, 1e-6, 1e-3};
        if (const char* collapsed = std::getenv("ORCA_FIXED_HISTORY_COLLAPSED_ENVELOPE")) {
            REQUIRE(std::string(collapsed) == "1");
            scenario.collapsed_envelope_allowance = CollapsedEnvelopeAllowance{1e-12, 1e-11, 1e-6, .1};
            std::cout << "Explicit collapsed-envelope omission: per_span_m3=1e-12 total_m3=1e-11"
                         " total_fraction=1e-6 neighbor_fraction=0.1; NOT a stress-error bound\n" << std::flush;
        }
        SectionGCodeLimits limits;
        if (std::getenv("ORCA_FIXED_HISTORY_LUMPED_THERMAL"))
            limits.assembly.quadrature=SectionIntegration::ThermalQuadrature::LumpedDiagonal;
        limits.max_segments      = 2000000;
        const auto prefix_paths=timeline.deposition.size()+timeline.unresolved_deposition.size();
        constexpr std::size_t prefix_visit_ceiling=200000000;
        // One reference traversal plus one traversal per selected snapshot.
        // Keep the original allowance for the default single-frame control.
        REQUIRE(prefix_paths<=prefix_visit_ceiling/(activation_count+1));
        limits.max_prefix_visits=std::max<std::size_t>(4000000,prefix_paths*(activation_count+1));
        // Explicit full-project offline work envelope; the ten-minute
        // checkpoint deadline and independent live-storage caps remain active.
        // These allowances are not normal slicing defaults.
        limits.streaming_thermal_work  = 2000000000;
        // Smaller spatial tiles reduce local arrangement complexity without
        // changing the nominal beads. Their additional index buckets remain
        // subject to the independent reference, point and work ceilings.
        // Include the 5 mm control: the shared planner and material cursors
        // consume one aggregate index allowance, not a fresh budget each.
        if (scenario.stream_tile_size_m->x()<=.005)
            limits.streaming_geometry.max_tiles=16384;
        if (const char* increment=std::getenv("ORCA_FIXED_HISTORY_PARTITION_UM"))
            limits.streaming_geometry.partition_increment_m=std::stod(increment)*1e-6;
        // Whole-project traversal processes more than8M boundary points, but
        // releases each tile packet. These cumulative counters are not live
        // storage limits.
        limits.streaming_geometry.max_boundary_points = 32000000;
        limits.streaming_geometry.max_index_work = 100000000;
        limits.streaming_thermal_points = 32000000;
        // The exact surface traversal has its own bounded CPU allowance.
        limits.streaming_skin.max_work = 2000000000;
        // Exact full-case frontiers exceed 100k traces after ownership-aware
        // compaction (~1.51 GiB measured process RSS). Use the existing 250k
        // offline ceiling: at most ~16 MiB of Seam payload, not unbounded growth.
        // This hidden diagnostic does not alter normal slicing defaults.
        limits.streaming_skin.max_seam_patches = 250000;
        // Provisional thermal moments retain source labels until the complete
        // certificate merges them. Bound this offline store at the existing
        // 400k ceiling (~83 MiB of 27-double moment payload per map); final
        // physical cell/group limits remain unchanged.
        limits.streaming_thermal_groups = 400000;
        limits.streaming_growth.max_work = 400000000;
        limits.streaming_growth.max_output_points = 128000000;
        // This is cumulative processed slab count, not retained geometry.
        // The topology's live-point, label, and frontier limits remain unchanged.
        limits.streaming_topology.max_slabs = 250000;
        limits.streaming_topology.max_processed_points = 128000000;
        limits.streaming_topology.max_work = 400000000;
        if (std::getenv("ORCA_FIXED_HISTORY_LONG_TOPOLOGY")) {
            // Explicit offline CPU budget only. The measured 400M-work run
            // takes about 26 seconds; retain all live-memory caps and the
            // ten-minute cancellation deadline for this longer validation.
            limits.streaming_topology.max_work = 1600000000;
        }
        if (std::getenv("ORCA_FIXED_HISTORY_EXTENDED_MECHANICS")) {
            // Explicit offline sparse study: bounded grid/assembly storage,
            // with the existing factor preflight capped near 96 MiB of entries.
            // The 10-minute deadline and live polygon caps remain unchanged.
            limits.mechanics.sparse_factor=SectionFactorLimits{};
            limits.mechanics.max_cells=24000;
            limits.mechanics.max_nodes=32000;
            limits.mechanics.max_free_dofs=96000;
            limits.mechanics.max_work=32000000;
            limits.mechanics.max_triplets=12000000;
            limits.streaming_integration.max_work=2000000000;
            limits.streaming_integration.max_total_points=32000000;
            limits.streaming_integration.max_groups=125000;
            // Cumulative traversal count, not retained polygons. Use the
            // existing sparse-study ceiling for the finer partition control;
            // live points, component storage, factor entries and time stay capped.
            limits.streaming_integration.max_pieces=500000;
            limits.streaming_integration.max_provisional_components=250000;
        }
        limits.mechanics.allow_free_components=std::getenv("ORCA_FIXED_HISTORY_ALLOW_FREE") != nullptr;
        // Preserve one-nanometre source intervals rather than merging them or
        // discarding their material; use their exact half-native midplanes.
        limits.streaming_geometry.geometry.half_native_section_planes = true;
        // Explicit 4nm Boolean representation sensitivity, including exact
        // predicates. Material-priority variants require the explicit test
        // environment option; they are not a physical ownership decision.
        limits.streaming_geometry.geometry.ownership_boundary_grid_units = 4;
        limits.ownership.geometry.ownership_boundary_grid_units          = 4;
        const bool exact_sections=std::getenv("ORCA_FIXED_HISTORY_EXACT_SECTIONS")!=nullptr;
        limits.streaming_geometry.geometry.retain_exact_interface_sections=exact_sections;
        limits.ownership.geometry.retain_exact_interface_sections=exact_sections;
        // Match the existing exact interface operation's4M ceiling, without
        // raising its vertex/storage bounds or the prepared total-work cap.
        limits.ownership.max_polygon_work = 4000000;
        limits.ownership.max_source_work = 4000000;
        if (const char* value=std::getenv("ORCA_FIXED_HISTORY_OWNERSHIP_WORK")) {
            std::size_t consumed=0;
            const auto budget=std::stoull(value,&consumed);
            REQUIRE(consumed==std::string(value).size());
            REQUIRE(budget>=80000000);
            REQUIRE(budget<=1000000000);
            // Existing prepared-resolver CPU ceiling, explicitly selected for
            // offline full-project work. Local geometry, live storage and the
            // independent cancellation deadline remain unchanged.
            limits.ownership.max_total_work=std::size_t(budget);
            std::cout << "Offline prepared ownership work=" << budget << '\n';
        }
        // Also budget a separate exact occupied-union topology certificate
        // when integer Boolean output has an ambiguous sub-grid channel.
        limits.streaming_geometry.geometry.max_slab_group_visits = 8000000;
        const char* interface_setting                                    = std::getenv("ORCA_FIXED_HISTORY_INTERFACE");
        const std::string interface_variant                              = interface_setting ? interface_setting : "strict";
        limits.ownership.allow_single_sided_native_boundary=std::getenv("ORCA_FIXED_HISTORY_SINGLE_SIDED_BOUNDARY")!=nullptr;
        if (const char* dual=std::getenv("ORCA_FIXED_HISTORY_DUAL_BOUNDARY")) {
            REQUIRE(std::string(dual)=="1");
            REQUIRE(exact_sections);
            limits.ownership.allow_dual_native_boundary=true;
        }
        if (const char* failure_path=std::getenv("ORCA_FIXED_HISTORY_FAILURE_FIXTURE")) {
            limits.ownership.observe_partition_failure=[path=std::string(failure_path),captured=std::make_shared<bool>(false)]
                (const NativeOwnershipPartitionFailure& failure) {
                if (*captured) return;
                *captured=true;
                try {
                    std::size_t vertices=0;
                    const auto ring=[&](const Polygon& polygon) {
                        nlohmann::json points=nlohmann::json::array();
                        for (const auto& point:polygon.points) {
                            if (++vertices>8192) throw std::runtime_error("Failure fixture vertex limit exceeded.");
                            points.push_back({point.x(),point.y()});
                        }
                        return points;
                    };
                    const auto polygons=[&](const ExPolygons& input) {
                        const auto begin=vertices;
                        nlohmann::json output=nlohmann::json::array();
                        for (const auto& polygon:input) {
                            nlohmann::json holes=nlohmann::json::array();
                            for (const auto& hole:polygon.holes) holes.push_back(ring(hole));
                            output.push_back({{"contour",ring(polygon.contour)},{"holes",std::move(holes)}});
                        }
                        if (vertices-begin>2048) throw std::runtime_error("Failure fixture operand limit exceeded.");
                        return output;
                    };
                    nlohmann::json fixture{{"schema",1},{"z_m",failure.sampled.conflict().bottom_m},
                        {"bed_z_m",failure.sampled.bed_z_m()},{"allowance_nm",7},
                        {"choice",failure.variant==NativeInterfaceVariant::FavorFirst ? "first" : "second"},
                        {"native_fill","nonzero-winding"},{"max_work",failure.max_work},
                        {"single_sided",failure.allow_single_sided_native_boundary},
                        {"dual_boundary",failure.allow_dual_native_boundary},
                        {"exact",failure.sampled.geometry_limits().retain_exact_interface_sections},
                        {"half_native_planes",failure.sampled.geometry_limits().half_native_section_planes},
                        {"materials",{failure.sampled.original_regions()[0].material,failure.sampled.original_regions()[1].material}},
                        {"failure",failure.message}};
                    fixture["native"]={polygons(failure.native_masks[0]),polygons(failure.native_masks[1])};
                    fixture["envelopes"]={polygons(failure.envelopes[0]),polygons(failure.envelopes[1])};
                    if (failure.sampled.clip()) {
                        const auto bounds=get_extents(failure.sampled.clip()->contour.points);
                        fixture["clip_bounds_nm"]={{bounds.min.x(),bounds.min.y()},{bounds.max.x(),bounds.max.y()}};
                    }
                    const auto payload=fixture.dump();
                    if (payload.size()>1024*1024) throw std::runtime_error("Failure fixture byte limit exceeded.");
                    // Exclusive creation also protects a concurrently created
                    // destination or symlink; a preflight exists() cannot.
                    std::unique_ptr<FILE,decltype(&std::fclose)> output(boost::nowide::fopen(path.c_str(),"wbx"),&std::fclose);
                    if (!output) throw std::runtime_error("Failure fixture destination exists or cannot be created.");
                    if (std::fwrite(payload.data(),1,payload.size(),output.get())!=payload.size())
                        throw std::runtime_error("Failure fixture write failed.");
                    if (std::fclose(output.release())!=0) throw std::runtime_error("Failure fixture close failed.");
                    std::cout<<"Ownership failure fixture saved: "<<path<<", vertices="<<vertices<<", bytes="<<payload.size()<<'\n'<<std::flush;
                } catch (const std::exception& error) {
                    std::cerr<<"Ownership failure fixture not saved: "<<error.what()<<'\n'<<std::flush;
                }
            };
        }
        if (interface_variant == "first")
            limits.ownership.interface_variant = NativeInterfaceVariant::FavorFirst;
        else if (interface_variant == "second")
            limits.ownership.interface_variant = NativeInterfaceVariant::FavorSecond;
        else if (interface_variant != "both")
            REQUIRE(interface_variant == "strict");
        if (const char* deposited = std::getenv("ORCA_FIXED_HISTORY_FIRST_DEPOSITION")) {
            REQUIRE(std::string(deposited) == "1");
            REQUIRE_FALSE(homogeneous);
            REQUIRE(interface_variant == "strict");
            REQUIRE(exact_sections);
            REQUIRE_FALSE(limits.ownership.allow_single_sided_native_boundary);
            REQUIRE_FALSE(limits.ownership.allow_dual_native_boundary);
            REQUIRE_FALSE(scenario.collapsed_envelope_allowance);
            scenario.use_native_ownership = false;
            scenario.use_deposited_ownership = true;
            limits.ownership.source_centerline_boundary_allowance_nm = 0;
            limits.ownership.deposited_overlap_policy = DepositedOverlapPolicy::FirstDepositionWins;
            std::cout << "Explicit first-deposition envelope ownership: authenticated command order, no native boundary assumption, "
                         "existing explicit bounded stationary allowance retained, no moving-vertical omission; "
                         "exact union refers to retained envelopes only, NOT physical displacement or bonding validation\n" << std::flush;
        }
        if (const char* probe=std::getenv("ORCA_FIXED_HISTORY_ROLE_PROBE")) {
            REQUIRE(std::string(probe)=="1");
            REQUIRE_FALSE(homogeneous);
            REQUIRE_FALSE(scenario.collapsed_envelope_allowance);
            REQUIRE_FALSE(std::getenv("ORCA_FIXED_HISTORY_ORDER_PROBE"));
            REQUIRE(scenario.origin_m.head<2>().isZero());
            REQUIRE_THAT(scenario.stream_tile_size_m->x(),Catch::Matchers::WithinAbs(.005,1e-12));
            REQUIRE_THAT(scenario.stream_tile_size_m->y(),Catch::Matchers::WithinAbs(.005,1e-12));
            const auto began=std::chrono::steady_clock::now();
            const auto check_probe=[&] {
                if (std::chrono::steady_clock::now()-began>std::chrono::seconds(60))
                    throw ResourceLimit("Role probe exceeded sixty seconds.");
            };
            const auto nominal=resolve_nominal_deposition(result,timeline,timeline.total_s,
                scenario.section,scenario.stationary_allowance,check_probe);
            REQUIRE(nominal.segments.size()==timeline.deposition.size());
            REQUIRE(nominal.retained_fractions.size()==nominal.segments.size());
            std::vector<DepositedSegment> selected;
            std::vector<ExtrusionRole> roles;
            for (std::size_t i=0;i<nominal.segments.size();++i) {
                check_probe();
                REQUIRE_THAT(nominal.retained_fractions[i],Catch::Matchers::WithinAbs(1.,0.));
                const auto& bead=nominal.segments[i];
                REQUIRE(bead.start_m.allFinite());
                REQUIRE(bead.end_m.allFinite());
                REQUIRE(std::isfinite(bead.width_m));
                REQUIRE(bead.width_m>0.);
                REQUIRE(bead.width_m<=10.);
                const auto width=std::llround(bead.width_m*1e9);
                const double radius=(std::max(double((width+1)/2),std::ceil(double(float(.5*double(width)))))+2.)*1e-9;
                if (std::max(bead.start_m.x(),bead.end_m.x())+radius<.050 ||
                    std::min(bead.start_m.x(),bead.end_m.x())-radius>.055 ||
                    std::max(bead.start_m.y(),bead.end_m.y())+radius<.110 ||
                    std::min(bead.start_m.y(),bead.end_m.y())-radius>.115) continue;
                REQUIRE(selected.size()<100000);
                const auto move=timeline.deposition[i].move_index;
                REQUIRE(move<result.moves.size());
                REQUIRE(result.moves[move].extrusion_role==timeline.deposition[i].role);
                selected.push_back(bead);
                roles.push_back(timeline.deposition[i].role);
            }
            REQUIRE_FALSE(selected.empty());
            // The complete-reference traversal uses two identical source
            // lists. Keep its shared interval planner and all thin endpoints.
            LocalSectionSamplingSchedule schedule({&selected,&selected},scenario.origin_m.head<2>(),
                *scenario.stream_tile_size_m,scenario.bed_z_m,scenario.stream_section_height_m,
                limits.streaming_geometry,check_probe);
            const auto& intervals=schedule.intervals({10,22});
            const ExPolygon tile(Polygon{{50000000,110000000},{55000000,110000000},
                {55000000,115000000},{50000000,115000000}});
            auto geometry=limits.streaming_geometry.geometry;
            geometry.exact_per_command_composition=true;
            geometry.retain_exact_interface_sections=true;
            geometry.max_boundary_points=4096;
            geometry.max_slab_group_visits=4000000;
            struct Measurement { SectionSamplerWork work; std::size_t samples=0, commands=0; };
            // Kind -2: complete material; -1: without internal infill;
            // nonnegative: that role alone. These are NOT additive costs.
            std::map<std::pair<std::size_t,int>,Measurement> measurements;
            std::size_t spent=0, visits=0, completed=0;
            coord_t last_top=0;
            bool truncated=false;
            std::cout << "Role counterfactual probe tile=10,22 selected=" << selected.size()
                      << " planned_bands=" << intervals.size()
                      << "; kind -2=complete, -1=without internal infill, others=ExtrusionRole enum"
                      << "; internal_infill_role=" << int(erInternalInfill)
                      << "; role-only composition is NOT additive attribution; no ownership or physics qualification\n";
            try {
                for (const auto& interval:intervals) {
                    if (interval.first>=8190000) break;
                    if (completed==512) { truncated=true; break; }
                    const double plane=(double(interval.first)+.5*double(interval.second-interval.first))*1e-9;
                    std::map<std::pair<std::size_t,int>,std::vector<DepositedSegment>> groups;
                    std::size_t active=0;
                    for (std::size_t i=0;i<selected.size();++i) {
                        check_probe();
                        if (++visits>4000000) throw ResourceLimit("Role probe source-visit budget exceeded.");
                        const auto& bead=selected[i];
                        if (plane+2e-9<std::min(bead.start_m.z()-bead.height_m,
                                bead.end_m.z()-bead.end_height_m.value_or(bead.height_m)) ||
                            plane-2e-9>=std::max(bead.start_m.z(),bead.end_m.z())) continue;
                        if (++active>4096) throw ResourceLimit("Role probe active command budget exceeded.");
                        groups[{bead.material,-2}].push_back(bead);
                        groups[{bead.material,int(roles[i])}].push_back(bead);
                        groups.try_emplace({bead.material,-1});
                        if (roles[i]!=erInternalInfill) groups[{bead.material,-1}].push_back(bead);
                    }
                    // Commit only a complete band, so reported variants cover
                    // exactly the same completed interval prefix on failure.
                    std::map<std::pair<std::size_t,int>,Measurement> band;
                    for (const auto& group:groups) {
                        if (group.second.empty()) { band[group.first]={{},1,0}; continue; }
                        if (spent>100000000-geometry.max_slab_group_visits)
                            throw ResourceLimit("Role probe cumulative work reservation exceeded.");
                        spent+=geometry.max_slab_group_visits;
                        const auto section=section_deposition(group.second,plane,scenario.bed_z_m,
                            geometry,check_probe,nullptr,&tile);
                        spent-=geometry.max_slab_group_visits-section.ownership_work;
                        band[group.first]={section.sampler_work,1,group.second.size()};
                    }
                    for (const auto& entry:band) {
                        auto& total=measurements[entry.first];
                        total.work.exact_construction+=entry.second.work.exact_construction;
                        total.work.exact_composition+=entry.second.work.exact_composition;
                        total.work.other+=entry.second.work.other;
                        total.samples+=entry.second.samples;
                        total.commands+=entry.second.commands;
                    }
                    ++completed;
                    last_top=interval.second;
                }
            } catch (const ResourceLimit& error) {
                truncated=true;
                std::cout << "Role probe stopped=" << error.what() << '\n';
            }
            for (const auto& entry:measurements)
                std::cout << "Role counterfactual material=" << entry.first.first << " kind=" << entry.first.second
                          << " samples=" << entry.second.samples << " command_visits=" << entry.second.commands
                          << " construction=" << entry.second.work.exact_construction
                          << " composition=" << entry.second.work.exact_composition << " other=" << entry.second.work.other << '\n';
            std::cout << "Role probe completed_bands=" << completed << " last_top_nm=" << last_top
                      << " truncated=" << truncated << " charged_or_reserved_work=" << spent
                      << " source_visits=" << visits << " seconds="
                      << std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count() << '\n' << std::flush;
            REQUIRE(completed>0);
            return;
        }
        const char* grouped_probe = std::getenv("ORCA_FIXED_HISTORY_GROUPED_WALL_PROBE");
        if (const char* probe=grouped_probe ? grouped_probe : std::getenv("ORCA_FIXED_HISTORY_ORDER_PROBE")) {
            REQUIRE(std::string(probe)=="1");
            if (grouped_probe) REQUIRE_FALSE(std::getenv("ORCA_FIXED_HISTORY_ORDER_PROBE"));
            REQUIRE_FALSE(homogeneous);
            REQUIRE_FALSE(scenario.collapsed_envelope_allowance);
            const auto probe_started=std::chrono::steady_clock::now();
            const auto probe_check=[&] {
                if (std::chrono::steady_clock::now()-probe_started>std::chrono::seconds(60))
                    throw ResourceLimit("Tile ordering probe exceeded sixty seconds.");
            };
            // Same full generated source and nominal resolver as history, but
            // no source-ownership qualification or mechanics in this probe.
            const auto nominal=resolve_nominal_deposition(result,timeline,timeline.total_s,
                scenario.section,scenario.stationary_allowance,probe_check);
            if (grouped_probe) {
                REQUIRE(nominal.segments.size() == timeline.deposition.size());
                REQUIRE(nominal.retained_fractions.size() == nominal.segments.size());
                for (std::size_t i = 0; i < nominal.segments.size(); ++i) {
                    probe_check();
                    REQUIRE_THAT(nominal.retained_fractions[i], Catch::Matchers::WithinAbs(1., 0.));
                    const auto& bead = nominal.segments[i];
                    const auto& entry = timeline.deposition[i];
                    REQUIRE_THAT((bead.start_m - entry.geometry.start_m).norm(), Catch::Matchers::WithinAbs(0., 0.));
                    REQUIRE_THAT((bead.end_m - entry.geometry.end_m).norm(), Catch::Matchers::WithinAbs(0., 0.));
                    REQUIRE_THAT(bead.width_m, Catch::Matchers::WithinAbs(entry.geometry.width_m, 0.));
                    REQUIRE_THAT(bead.height_m, Catch::Matchers::WithinAbs(entry.geometry.height_m, 0.));
                    REQUIRE(bead.end_height_m.has_value() == entry.geometry.end_height_m.has_value());
                    if (bead.end_height_m)
                        REQUIRE_THAT(*bead.end_height_m, Catch::Matchers::WithinAbs(*entry.geometry.end_height_m, 0.));
                    REQUIRE(bead.material == entry.geometry.material);
                    REQUIRE(bead.section == scenario.section);
                    REQUIRE(entry.move_index < result.moves.size());
                    REQUIRE(result.moves[entry.move_index].extrusion_role == entry.role);
                    REQUIRE(result.moves[entry.move_index].extruder_id == bead.material);
                }
            }
            const std::vector<double> planes = grouped_probe ? std::vector<double>{.0041, .00415, .007182184} :
                                                               std::vector<double>{.007182184};
            const ExPolygon tile(Polygon{{50000000,110000000},{55000000,110000000},
                {55000000,115000000},{50000000,115000000}});
            for (const double plane : planes) {
            std::map<std::size_t,std::vector<DepositedSegment>> groups;
            std::map<std::size_t,std::vector<std::size_t>> source_indices;
            std::size_t selected=0;
            for (std::size_t source_index = 0; source_index < nominal.segments.size(); ++source_index) {
                const auto& bead = nominal.segments[source_index];
                probe_check();
                REQUIRE(bead.start_m.allFinite());
                REQUIRE(bead.end_m.allFinite());
                REQUIRE(std::isfinite(bead.width_m));
                REQUIRE(bead.width_m>0.);
                REQUIRE(bead.width_m<=10.);
                const double end_height=bead.end_height_m.value_or(bead.height_m);
                REQUIRE(std::isfinite(bead.height_m));
                REQUIRE(std::isfinite(end_height));
                // Outward padding includes native rounding near XY/Z bounds;
                // the ordinary sampler still makes the actual clipped cut.
                constexpr double pad=2e-9;
                if (plane+pad<std::min(bead.start_m.z()-bead.height_m,bead.end_m.z()-end_height) ||
                    plane-pad>=std::max(bead.start_m.z(),bead.end_m.z())) continue;
                // Match LocalSections' outward filter, including the native
                // width rounding and float radius actually passed to offset.
                const auto width=std::llround(bead.width_m*1e9);
                const double radius=(std::max(double((width+1)/2),
                    std::ceil(double(float(.5*double(width)))))+2.)*1e-9;
                if (std::max(bead.start_m.x(),bead.end_m.x())+radius<.050 ||
                    std::min(bead.start_m.x(),bead.end_m.x())-radius>.055 ||
                    std::max(bead.start_m.y(),bead.end_m.y())+radius<.110 ||
                    std::min(bead.start_m.y(),bead.end_m.y())-radius>.115) continue;
                REQUIRE(++selected<=4096);
                groups[bead.material].push_back(bead);
                source_indices[bead.material].push_back(source_index);
            }
            if (grouped_probe && selected == 0) {
                std::cout << "Grouped wall plane_m=" << plane << " candidates=0\n";
                continue;
            }
            REQUIRE(selected>0);
            auto geometry=limits.streaming_geometry.geometry;
            geometry.exact_per_command_composition=true;
            geometry.retain_exact_interface_sections=true;
            geometry.max_slab_group_visits=4000000;
            geometry.max_boundary_points=4096;
            std::cout << "Tile order probe tile=10,22 plane_m=" << plane << " candidates=" << selected
                      << " stationary_omitted_m3=" << nominal.stationary_extrusion.commanded_volume_m3
                      << " geometry-only; no source qualification\n";
            for (const auto& group:groups) {
                const auto& original=group.second;
                if (grouped_probe) {
                    const auto& indices = source_indices.at(group.first);
                    const auto native = [](double value) { return std::llround(value * 1e9); };
                    // Nonexclusive counters explain every rejected condition;
                    // they do not count calls made while extending a chain.
                    std::array<std::size_t, 6> rejected{};
                    double maximum_z_change_m = 0., maximum_height_change_m = 0.;
                    for (std::size_t i = 0; i < original.size(); ++i) {
                        const auto& bead = original[i];
                        const auto& entry = timeline.deposition.at(indices[i]);
                        const auto& move = result.moves.at(entry.move_index);
                        rejected[0] += entry.role != erPerimeter && entry.role != erExternalPerimeter;
                        rejected[1] += move.deposition_source_id == 0;
                        rejected[2] += native(bead.start_m.z()) != native(bead.end_m.z());
                        rejected[3] += native(bead.height_m) != native(bead.end_height_m.value_or(bead.height_m));
                        rejected[4] += bead.section != BeadSection::RectangularEnvelope && bead.section != BeadSection::RoundedRectangle;
                        rejected[5] += native(plane) < native(bead.start_m.z()) - native(bead.height_m) ||
                            native(plane) >= native(bead.start_m.z());
                        maximum_z_change_m = std::max(maximum_z_change_m, std::abs(bead.start_m.z() - bead.end_m.z()));
                        maximum_height_change_m = std::max(maximum_height_change_m,
                            std::abs(bead.height_m - bead.end_height_m.value_or(bead.height_m)));
                    }
                    const auto eligible = [&](std::size_t i) {
                        const auto& bead = original[i];
                        const auto& entry = timeline.deposition.at(indices[i]);
                        const auto& move = result.moves.at(entry.move_index);
                        return (entry.role == erPerimeter || entry.role == erExternalPerimeter) &&
                            (bead.section == BeadSection::RectangularEnvelope || bead.section == BeadSection::RoundedRectangle) &&
                            move.deposition_source_id != 0 &&
                            native(plane) >= native(bead.start_m.z()) - native(bead.height_m) &&
                            native(plane) < native(bead.start_m.z()) &&
                            native(bead.start_m.z()) == native(bead.end_m.z()) &&
                            native(bead.height_m) == native(bead.end_height_m.value_or(bead.height_m));
                    };
                    const auto adjacent = [&](std::size_t a, std::size_t b) {
                        const auto& first = original[a];
                        const auto& second = original[b];
                        const auto& ea = timeline.deposition.at(indices[a]);
                        const auto& eb = timeline.deposition.at(indices[b]);
                        const auto& ma = result.moves.at(ea.move_index);
                        const auto& mb = result.moves.at(eb.move_index);
                        return indices[b] == indices[a] + 1 && eligible(b) && ea.role == eb.role &&
                            first.section == second.section &&
                            ma.deposition_source_id == mb.deposition_source_id && ma.extruder_id == mb.extruder_id &&
                            native(first.width_m) == native(second.width_m) &&
                            native(first.height_m) == native(second.height_m) &&
                            native(first.end_m.z()) == native(second.start_m.z()) &&
                            native(first.end_m.x()) == native(second.start_m.x()) &&
                            native(first.end_m.y()) == native(second.start_m.y());
                    };
                    std::vector<std::pair<std::size_t, std::size_t>> runs;
                    std::size_t eligible_commands = 0, grouped_commands = 0;
                    for (std::size_t i = 0; i < original.size();) {
                        std::size_t end = i + 1;
                        if (eligible(i)) {
                            ++eligible_commands;
                            while (end < original.size() && adjacent(end - 1, end)) { ++end; ++eligible_commands; }
                        }
                        if (end - i > 1) grouped_commands += end - i;
                        runs.emplace_back(i, end);
                        i = end;
                    }
                    std::cout << "Grouped wall eligibility material=" << group.first << " plane_m=" << plane
                              << " commands=" << original.size() << " eligible=" << eligible_commands
                              << " grouped_commands=" << grouped_commands << " runs=" << runs.size()
                              << " rejection_nonexclusive_nonwall=" << rejected[0] << " missing_source=" << rejected[1]
                              << " nonhorizontal=" << rejected[2] << " tapered=" << rejected[3]
                              << " unsupported_shape=" << rejected[4] << " outside=" << rejected[5]
                              << " maximum_z_change_m=" << maximum_z_change_m
                              << " maximum_height_change_m=" << maximum_height_change_m
                              << " geometry diagnostic only, no ownership/contact qualification\n" << std::flush;
                    const auto baseline_start = std::chrono::steady_clock::now();
                    ExactOverlayDiagnostics baseline_diagnostics;
                    const auto baseline = section_deposition(original, plane, scenario.bed_z_m, geometry, probe_check,
                        nullptr, &tile, &baseline_diagnostics);
                    const double baseline_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - baseline_start).count();
                    REQUIRE(baseline.regions.size() == 1);
                    const auto tile_exact = make_exact_section(ExPolygons{tile}, 4000000, probe_check);
                    ExactSection combined;
                    std::size_t work = tile_exact.work;
                    const auto charge = [&](std::size_t value) {
                        if (value > 4000000 - work) throw ResourceLimit("Grouped wall comparison work budget exceeded.");
                        work += value;
                    };
                    double maximum_width_shift_m = 0.;
                    const auto grouped_start = std::chrono::steady_clock::now();
                    for (const auto& run : runs) {
                        probe_check();
                        ExactSection cut;
                        if (run.second - run.first > 1) {
                            const auto& first = original[run.first];
                            double radius = .5 * native(first.width_m);
                            if (first.section == BeadSection::RoundedRectangle) {
                                const double depth = double(native(first.end_m.z()) - native(plane));
                                const double height = double(native(first.height_m));
                                radius = .5 * (native(first.width_m) - height) + std::sqrt(std::max(0., depth * (height - depth)));
                            }
                            if (radius <= 0.) continue;
                            const double width = double(std::llround(2. * radius)) * 1e-9;
                            maximum_width_shift_m = std::max(maximum_width_shift_m, std::abs(width - 2. * radius * 1e-9));
                            std::vector<DepositedSegment> proxies;
                            for (std::size_t i = run.first; i < run.second; ++i) {
                                auto bead = original[i];
                                bead.start_m.z() = bead.end_m.z() = .0002;
                                bead.width_m = width;
                                bead.height_m = .0002;
                                bead.end_height_m.reset();
                                bead.section = BeadSection::RectangularEnvelope;
                                proxies.push_back(bead);
                            }
                            auto native_limits = geometry;
                            native_limits.retain_exact_interface_sections = false;
                            native_limits.exact_per_command_composition = false;
                            const auto layers = layered_deposition(proxies, 0., native_limits, probe_check);
                            REQUIRE(layers.slabs.size() == 1);
                            const auto exact = make_exact_section(layers.slabs.front().footprint.native(), 4000000 - work, probe_check);
                            charge(exact.work);
                            const auto clipped = boolean_exact_sections(exact.section, tile_exact.section,
                                ExactSectionOperation::Intersection, 4000000 - work, probe_check);
                            charge(clipped.work);
                            cut = clipped.section;
                        } else {
                            const auto single = section_deposition({original[run.first]}, plane, scenario.bed_z_m,
                                geometry, probe_check, nullptr, &tile);
                            charge(single.ownership_work);
                            if (!single.regions.empty()) cut = single.regions.front().footprint.exact();
                        }
                        const auto joined = boolean_exact_sections(combined, cut, ExactSectionOperation::Union, 4000000 - work, probe_check);
                        charge(joined.work);
                        combined = joined.section;
                    }
                    const double grouped_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - grouped_start).count();
                    const auto& reference = baseline.regions.front().footprint.exact();
                    const auto difference = boolean_exact_sections(reference, combined, ExactSectionOperation::SymmetricDifference,
                        4000000 - work, probe_check);
                    charge(difference.work);
                    const auto area = exact_section_horizontal_area(difference.section,
                        {50000000,110000000,0}, {5000000,5000000,20000000}, native(plane), 4000000 - work, probe_check);
                    charge(area.work);
                    const auto before = exact_section_interior_components(reference, 4000000 - work, probe_check);
                    charge(before.work);
                    const auto after = exact_section_interior_components(combined, 4000000 - work, probe_check);
                    charge(after.work);
                    std::cout << "Grouped wall material=" << group.first << " baseline_s=" << baseline_seconds
                              << " grouped_s=" << grouped_seconds << " baseline_work=" << baseline.ownership_work
                              << " baseline_successful_overlay_debits_not_full_cost=1"
                              << " overlays=" << baseline_diagnostics.overlays << " edges=" << baseline_diagnostics.edges
                              << " candidates=" << baseline_diagnostics.candidates
                              << " boxes_projection=" << baseline_diagnostics.boxes_projection
                              << " sorting=" << baseline_diagnostics.sorting << " sweep=" << baseline_diagnostics.sweep
                              << " active_scans=" << baseline_diagnostics.active_scans
                              << " event_reservation=" << baseline_diagnostics.event_reservation
                              << " diagnostic_overflow=" << baseline_diagnostics.overflow
                              << " grouped_and_comparison_exact_work=" << work << " native_group_work_uninstrumented=1"
                              << " xor_area_m2=" << area.value << " xor_bounds=" << area.bounds[0] << ',' << area.bounds[1]
                              << " exact_interior_components=" << before.sections.size() << ',' << after.sections.size()
                              << " width_quantization_m=" << maximum_width_shift_m
                              << " grouped_clip_is_exact_baseline_clip_native=1 topology_correspondence_unqualified=1\n" << std::flush;
                    continue;
                }
                std::vector<std::size_t> order;
                for (std::size_t i=0;i<original.size();++i) order.push_back(i);
                std::size_t comparisons=0;
                const auto sort_started=std::chrono::steady_clock::now();
                std::stable_sort(order.begin(),order.end(),[&](std::size_t a,std::size_t b) {
                    probe_check();
                    if (++comparisons>131072) throw ResourceLimit("Tile ordering comparison budget exceeded.");
                    const auto ca=.5*(original[a].start_m+original[a].end_m);
                    const auto cb=.5*(original[b].start_m+original[b].end_m);
                    return ca.x()<cb.x() || (ca.x()==cb.x() && ca.y()<cb.y());
                });
                const double sort_seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-sort_started).count();
                std::vector<DepositedSegment> spatial;
                for (const auto index:order) spatial.push_back(original[index]);
                std::array<std::optional<DepositedSection>,2> sampled;
                for (unsigned variant=0;variant<2;++variant) {
                    const auto began=std::chrono::steady_clock::now();
                    try {
                        sampled[variant]=section_deposition(variant ? spatial : original,plane,scenario.bed_z_m,
                            geometry,probe_check,nullptr,&tile);
                        const auto& section=*sampled[variant];
                        std::cout << "Tile order material=" << group.first << " spatial=" << variant
                                  << " commands=" << original.size() << " construction=" << section.sampler_work.exact_construction
                                  << " composition=" << section.sampler_work.exact_composition << " other=" << section.sampler_work.other
                                  << " total=" << section.ownership_work << " vertices=" << section.boundary_points;
                    } catch (const ResourceLimit& error) {
                        std::cout << "Tile order material=" << group.first << " spatial=" << variant << " rejected=" << error.what();
                    }
                    std::cout << " seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count()
                              << " sort_comparisons=" << (variant ? comparisons : 0)
                              << " sort_seconds=" << (variant ? sort_seconds : 0.) << '\n' << std::flush;
                }
                REQUIRE(sampled[0]);
                REQUIRE(sampled[1]);
                REQUIRE(sampled[0]->regions.size()==sampled[1]->regions.size());
                if (!sampled[0]->regions.empty()) {
                    REQUIRE(sampled[0]->regions.size()==1);
                    CHECK(boolean_exact_sections(sampled[0]->regions.front().footprint.exact(),
                        sampled[1]->regions.front().footprint.exact(),ExactSectionOperation::SymmetricDifference,
                        4000000,probe_check).section.empty());
                }
            }
            }
            return;
        }
        const std::vector<NativeInterfaceVariant> variants = interface_variant == "both" ?
                                                                 std::vector<NativeInterfaceVariant>{NativeInterfaceVariant::FavorFirst,
                                                                                                     NativeInterfaceVariant::FavorSecond} :
                                                                 std::vector<NativeInterfaceVariant>{limits.ownership.interface_variant};
        const auto start                                   = std::chrono::steady_clock::now();
        auto reported                                      = start;
        std::cout << "Fixed history admission: full project, " << activation_count << " batched activations, synthetic400K "
                  << (cooling_control ? "heated-bed/ambient cooling, " : "insulated, ")
                  << scenario.stream_vertical_normal_bins << " angular bins, 4nm Boolean band; NOT curl validation\n"
                  << std::flush;
        std::cout << "Native interface sensitivity=" << interface_variant << '\n' << std::flush;
        std::cout << "Retained exact interface sections=" << exact_sections << '\n' << std::flush;
        std::cout << "Single-sided native boundary assumption=" << limits.ownership.allow_single_sided_native_boundary << '\n';
        if (limits.ownership.allow_dual_native_boundary)
            std::cout<<"Conditional dual-boundary material sensitivity enabled: assess BOTH priorities; unchanged deposited union, NOT a physical uncertainty bound\n";
        std::cout << "Geometry stream tile_mm=" << scenario.stream_tile_size_m->x() * 1000. << '\n' << std::flush;
        std::cout << "Analysis partition increment_um=" << limits.streaming_geometry.partition_increment_m*1e6
                  << "; requires independent refinement\n" << std::flush;
        std::cout << "Homogeneous numerical control=" << homogeneous << '\n' << std::flush;
        std::cout << "Free-body diagnostic=" << limits.mechanics.allow_free_components
                  << "; detached components are not bed adhesion\n" << std::flush;
        std::cout << "Synthetic cooling control=" << cooling_control << '\n' << std::flush;
        std::cout << "Thermal quadrature=" << int(limits.assembly.quadrature) << '\n' << std::flush;
        const auto check = [&] {
            const auto now = std::chrono::steady_clock::now();
            if (now - start > std::chrono::minutes(10))
                throw ResourceLimit("Admission diagnostic exceeded ten minutes.");
            if (now - reported > std::chrono::seconds(20)) {
                std::cout << "Fixed history admission active seconds=" << std::chrono::duration<double>(now - start).count() << '\n'
                          << std::flush;
                reported = now;
            }
        };
        for (const auto variant : variants) {
            limits.ownership.interface_variant = variant;
            limits.ownership.source_centerline_boundary_allowance_nm = variant == NativeInterfaceVariant::Strict ? 0 : 4;
            std::cout << "Native interface variant id=" << int(variant) << '\n' << std::flush;
            std::cout << "Analysis source centerline band_nm=" << limits.ownership.source_centerline_boundary_allowance_nm << '\n'
                      << std::flush;
            try {
                const auto field = simulate_section_gcode_history(result, scenario, limits, check);
                std::cout << "Fixed history admission thermal_state=" << int(field.thermal.state)
                          << " mechanical=" << bool(field.mechanical) << " volume_m3=" << field.occupied_volume_m3
                          << " excluded_m3=" << field.excluded_nominal_volume_m3
                          << " collapsed_omitted_count=" << field.collapsed_envelope.count
                          << " collapsed_omitted_m3=" << field.collapsed_envelope.commanded_volume_m3
                          << " accepted_steps=" << field.thermal.steps << " solver_steps=" << field.thermal.solver_steps << '\n';
                if (field.thermal.rejected_temperature) {
                    std::cout << "Rejected thermal node=" << field.thermal.rejected_temperature->node
                              << " kelvin=" << field.thermal.rejected_temperature->kelvin << '\n';
                    if (field.thermal.rejected_temperature->node<cooling_nodes.size()) {
                        const auto& node=cooling_nodes[field.thermal.rejected_temperature->node];
                        std::cout << "Rejected thermal location_m=" << node.position_m.transpose()
                                  << " component=" << node.component << '\n';
                    }
                }
                REQUIRE(field.thermal.state == ThermalHistory::State::Converged);
                REQUIRE(field.mechanical.has_value());
                CHECK_THAT(field.thermal.energy_balance_error_j, Catch::Matchers::WithinAbs(0., 1e-5));
                REQUIRE_FALSE(field.thermal.temperatures_k.empty());
                const auto extrema=std::minmax_element(field.thermal.temperatures_k.begin(),field.thermal.temperatures_k.end());
                std::cout << "Fixed history thermal: steps=" << field.thermal.steps
                          << " minimum_k=" << *extrema.first << " maximum_k=" << *extrema.second
                          << " energy_balance_error_j=" << field.thermal.energy_balance_error_j << '\n';
                std::cout << "Fixed history mechanics: nodes=" << field.mechanical->nodes_m.size()
                          << " solved_dofs=" << field.mechanical->solved_dofs
                          << " factor_entries=" << field.mechanical->factor_entries
                          << " free_components=" << field.mechanical->free_components
                          << " volume_m3=" << field.mechanical->occupied_volume_m3
                          << " bed_area_m2=" << field.mechanical->bed_area_m2
                          << " elastic_energy_j=" << field.mechanical->elastic_energy_j
                          << " residual=" << field.mechanical->relative_residual << '\n';
                if (cooling_control) {
                    const auto reaction=project_section_bed_reactions(*field.mechanical,scenario.origin_m.head<2>(),
                        scenario.spacing_m.head<2>(),check);
                    std::cout << "Synthetic cooling response: boundary_energy_j=" << field.thermal.boundary_energy_j
                              << " tensile_generalized_force_n=" << reaction.tensile_generalized_force_n << '\n';
                    CHECK(field.thermal.boundary_energy_j<0.);
                    CHECK(field.mechanical->elastic_energy_j>0.);
                } else {
                    // This control remains uniformly at its stress-free temperature.
                    CHECK_THAT(field.mechanical->elastic_energy_j, Catch::Matchers::WithinAbs(0., 1e-8));
                }
            } catch (const ResourceLimit& error) {
                std::cout << "Fixed history admission ResourceLimited; " << error.what() << '\n';
                FAIL_CHECK(error.what());
            } catch (const std::invalid_argument& error) {
                std::cout << "Fixed history admission unsupported; " << error.what() << '\n';
                FAIL_CHECK(error.what());
            }
        }
        std::cout << "Fixed history admission elapsed_s=" << std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count()
                  << "; no physical qualification\n"
                  << std::flush;
        return;
    }
    try {
        const auto geometry = Thermoelastic::source_geometry_diagnostics(result);
        if (geometry.state == Thermoelastic::SourceGeometryDiagnostics::State::Complete) {
            // Sources, horizontal moves, sloped moves, unresolved moves, and
            // sources whose parsed nominal height varies. Display-only timing
            // vertices are excluded by the diagnostic, unlike the raw count above.
            std::map<std::pair<unsigned, unsigned>, std::array<std::size_t, 5>> groups;
            for (std::size_t id = 1; id < geometry.sources.size(); ++id) {
                const auto* source = result.deposition_provenance->source(static_cast<std::uint32_t>(id));
                REQUIRE(source != nullptr);
                const auto& measured = geometry.sources[id];
                auto& group = groups[{source->logical_filament, static_cast<unsigned>(source->height)}];
                ++group[0];
                group[1] += measured.horizontal_moves;
                group[2] += measured.sloped_moves;
                group[3] += measured.stationary_or_unresolved_moves;
                group[4] += measured.height_varies() ? 1 : 0;
            }
            for (const auto& [key, group] : groups) {
                const char* height = key.second == static_cast<unsigned>(DepositionSource::Height::NativePlanar) ? "native_planar" :
                    key.second == static_cast<unsigned>(DepositionSource::Height::Transformed) ? "transformed" : "unknown";
                std::cout << "source_geometry logical_filament=" << key.first << " source_height=" << height
                          << " sources=" << group[0] << " horizontal=" << group[1] << " sloped=" << group[2]
                          << " stationary_or_unresolved=" << group[3] << " variable_height_sources=" << group[4] << '\n';
            }
            const auto& unowned = geometry.sources.front();
            std::cout << "source_geometry unowned_all_filaments horizontal=" << unowned.horizontal_moves
                      << " sloped=" << unowned.sloped_moves
                      << " stationary_or_unresolved=" << unowned.stationary_or_unresolved_moves << '\n';
        } else {
            std::cout << "Source geometry diagnostics unavailable: no verified registry.\n";
        }
    } catch (const Thermoelastic::ResourceLimit&) {
        std::cout << "Source geometry diagnostics unavailable: the two-million-move budget was exceeded.\n";
    }
    if (const char *probe = std::getenv("ORCA_SECTION_PROBE_M")) {
        char *end = nullptr;
        const double z = std::strtod(probe,&end);
        REQUIRE(end != probe);
        REQUIRE(*end == '\0');
        REQUIRE(std::isfinite(z));
        REQUIRE(z >= 0.);
        REQUIRE(result.deposition_sections);
        auto path = Thermoelastic::timed_deposition(result,2000000,2000000);
        audit_collapsed_envelopes(result, path);
        REQUIRE(path.state == Thermoelastic::ToolpathState::Complete);
        const auto excluded = [](const auto &entry) {
            return entry.role == erWipeTower || entry.role == erCustom || entry.role == erNone;
        };
        path.deposition.erase(std::remove_if(path.deposition.begin(),path.deposition.end(),excluded),path.deposition.end());
        path.unresolved_deposition.erase(std::remove_if(path.unresolved_deposition.begin(),path.unresolved_deposition.end(),excluded),
            path.unresolved_deposition.end());
        // Explicit small endpoint-extrusion approximation, reported separately
        // from nominal line geometry; no spatial extrusion is silently erased.
        const Thermoelastic::StationaryExtrusionAllowance allowance{1e-12,1e-6,1e-3};
        auto nominal = Thermoelastic::resolve_nominal_deposition(result,path,path.total_s,
            Thermoelastic::BeadSection::RoundedRectangle,allowance);
        // Explicit diagnostic scenario, e.g. 0,0,1 for two ASA slots and TPU.
        // Do not infer equal mechanical properties from names or colors.
        if (const char *mapping = std::getenv("ORCA_SECTION_MATERIAL_MAP")) {
            std::istringstream input(mapping);
            std::vector<std::size_t> materials;
            std::size_t material = 0;
            while (input >> material) {
                REQUIRE(material < 256);
                materials.push_back(material);
                REQUIRE(materials.size() <= 256);
                if (input.peek() == ',') input.get(); else break;
            }
            REQUIRE(input.eof());
            for (auto &segment : nominal.segments) {
                REQUIRE(segment.material < materials.size());
                segment.material = materials[segment.material];
            }
            std::cout << "Section probe explicit material map=" << mapping << '\n';
        }
        std::cout << "Section probe z_m=" << z << " beads=" << nominal.segments.size()
                  << " stationary_volume_m3=" << nominal.stationary_extrusion.commanded_volume_m3 << '\n';
        const auto section = Thermoelastic::section_deposition(nominal.segments,z,0.);
        std::cout << "Section probe area_m2=" << section.area_m2 << " regions=" << section.regions.size() << '\n';
        REQUIRE_FALSE(section.regions.empty());
    }
    if (const char* audit = std::getenv("ORCA_SECTION_LAYER_AUDIT")) {
        using namespace Thermoelastic;
        const double requested_z_mm = std::stod(audit);
        double distance = std::numeric_limits<double>::infinity();
        unsigned layer = 0;
        for (const auto& move : result.moves)
            if (move.type == EMoveType::Extrude && move.extrusion_role != erWipeTower &&
                move.extrusion_role != erCustom && move.extrusion_role != erNone &&
                std::abs(move.position.z() - requested_z_mm) < distance) {
                distance = std::abs(move.position.z() - requested_z_mm);
                layer = move.layer_id;
            }
        REQUIRE(distance < .001);
        const auto isotropic = [](double e, double nu) {
            const double mu = e / (2. * (1. + nu)), lambda = e * nu / ((1. + nu) * (1. - 2. * nu));
            ElasticityMatrix c = ElasticityMatrix::Zero();
            c.topLeftCorner<3,3>().setConstant(lambda);
            c.diagonal().head<3>().array() += 2. * mu;
            c.diagonal().tail<3>().setConstant(mu);
            return c;
        };
        GCodeSectionScenario scenario(BeadSection::RoundedRectangle);
        // Synthetic load for computational coverage, NOT a material calibration.
        scenario.materials = {{isotropic(2e9,.35)}, {isotropic(26e6,.45)}};
        for (auto& material : scenario.materials) material.eigenstrain.head<3>().setConstant(-.001);
        if (const char* value = std::getenv("ORCA_SECTION_STRAIN_Z_GRADIENT_PER_M")) {
            const std::string text(value);
            REQUIRE(text.size() <= 32);
            std::size_t consumed = 0;
            const double gradient = std::stod(text,&consumed);
            REQUIRE(consumed == text.size());
            REQUIRE(std::isfinite(gradient));
            REQUIRE(std::abs(gradient) <= 1.);
            for (auto& material : scenario.materials)
                material.eigenstrain_gradient_per_m.col(2).head<3>().setConstant(gradient);
            std::cout << "Synthetic prescribed normal-strain Z gradient per metre=" << gradient
                      << "; NOT a computed temperature history or material calibration\n";
        }
        scenario.filament_to_material = {0,0,1};
        scenario.initial_basis_spacing_m = {.02,.02,.005};
        const std::array<const char*,3> basis_variables{
            "ORCA_SECTION_BASIS_X_MM","ORCA_SECTION_BASIS_Y_MM","ORCA_SECTION_BASIS_Z_MM"};
        for (unsigned axis = 0; axis < basis_variables.size(); ++axis) {
            if (const char* spacing = std::getenv(basis_variables[axis])) {
                const double millimetres = std::stod(spacing);
                REQUIRE(std::isfinite(millimetres));
                REQUIRE(millimetres >= .001);
                REQUIRE(millimetres <= 20.);
                scenario.initial_basis_spacing_m[axis] = millimetres * .001;
            }
        }
        scenario.study.max_section_height_m = .0002;
        scenario.study.relative_tolerance = .05;
        scenario.study.max_total_sections = 40000;
        scenario.study.max_total_boundary_points = 8000000;
        scenario.study.max_total_segment_visits = 80000000;
        scenario.study.max_total_work = 1000000000;
        if (const char* tile = std::getenv("ORCA_SECTION_LOCAL_TILE_MM")) {
            const double size_m = std::stod(tile) * .001;
            REQUIRE(std::isfinite(size_m));
            REQUIRE(size_m > 0.);
            scenario.study.local_tile_size_m = Eigen::Vector2d(size_m, size_m);
        }
        if (const char* quantum = std::getenv("ORCA_SECTION_Z_QUANTIZATION_UM")) {
            const double increment_m = std::stod(quantum) * 1e-6;
            REQUIRE(std::isfinite(increment_m));
            REQUIRE(increment_m > 0.);
            scenario.study.initial_z_quantization_m = increment_m;
        }
        scenario.use_native_ownership = std::getenv("ORCA_SECTION_NATIVE_OWNERSHIP") != nullptr;
        scenario.study.sections.geometry.half_native_section_planes =
            std::getenv("ORCA_SECTION_HALF_NATIVE_PLANES") != nullptr;
        if (const char* band = std::getenv("ORCA_SECTION_OWNERSHIP_GRID_UNITS"))
            scenario.study.sections.geometry.ownership_boundary_grid_units = unsigned(std::stoul(band));
        scenario.stationary_extrusion_allowance = StationaryExtrusionAllowance{1e-12,1e-6,1e-3};
        if (std::getenv("ORCA_SECTION_SINGLE_MATERIAL_AUDIT")) {
            scenario.filament_to_material = {0,0,0};
            scenario.use_native_ownership = false;
            std::cout << "Single-material scaling diagnostic only, NOT an actual ASA/TPU stress or Auto result\n";
        }
        if (std::getenv("ORCA_SECTION_STREAM_PROBE")) {
            REQUIRE(std::getenv("ORCA_SECTION_SINGLE_MATERIAL_AUDIT") != nullptr);
            REQUIRE(scenario.study.local_tile_size_m.has_value());
            REQUIRE_FALSE(scenario.study.initial_z_quantization_m.has_value());
            auto path = timed_deposition(result, 2000000, 2000000);
            audit_collapsed_envelopes(result, path);
            REQUIRE(path.state == ToolpathState::Complete);
            const auto excluded = [](const auto& entry) {
                return entry.role == erWipeTower || entry.role == erCustom || entry.role == erNone;
            };
            double end = 0.;
            for (const auto& entry : path.deposition)
                if (result.moves.at(entry.move_index).layer_id == layer && !excluded(entry))
                    end = std::max(end, entry.end_s);
            REQUIRE(end > 0.);
            path.deposition.erase(std::remove_if(path.deposition.begin(), path.deposition.end(), excluded),
                path.deposition.end());
            path.unresolved_deposition.erase(
                std::remove_if(path.unresolved_deposition.begin(), path.unresolved_deposition.end(), excluded),
                path.unresolved_deposition.end());
            for (const auto& entry : path.deposition)
                if (entry.start_s < end) REQUIRE(entry.role != erMixed);
            auto nominal = resolve_nominal_deposition(result, path, end, scenario.section,
                scenario.stationary_extrusion_allowance);
            if (std::getenv("ORCA_SECTION_DETACHED_CONTACT_WITNESS")) {
                // Local source candidates only: never a component-membership or
                // physical bonding certificate. Keep whole original beads.
                constexpr double xmin = .184901674, xmax = .185762328;
                constexpr double ymin = .120555190, ymax = .121391814;
                constexpr double zbase = .0098, halo = .001;
                std::map<std::uint32_t, std::vector<Thermoelastic::DepositedSegment>> candidates;
                std::size_t selected = 0;
                for (std::size_t i = 0; i < nominal.segments.size(); ++i) {
                    const auto& bead = nominal.segments[i];
                    const double radius = .5 * bead.width_m + 2e-9;
                    const double lo = std::min(bead.start_m.z() - bead.height_m,
                        bead.end_m.z() - bead.end_height_m.value_or(bead.height_m));
                    const double hi = std::max(bead.start_m.z(), bead.end_m.z());
                    if (hi < zbase - 2e-9 || lo > zbase + 2e-9 ||
                        std::max(bead.start_m.x(), bead.end_m.x()) + radius < xmin - halo ||
                        std::min(bead.start_m.x(), bead.end_m.x()) - radius > xmax + halo ||
                        std::max(bead.start_m.y(), bead.end_m.y()) + radius < ymin - halo ||
                        std::min(bead.start_m.y(), bead.end_m.y()) - radius > ymax + halo) continue;
                    REQUIRE(++selected <= 4096);
                    const auto& move = result.moves.at(path.deposition.at(i).move_index);
                    REQUIRE(move.deposition_source_id != 0);
                    candidates[move.deposition_source_id].push_back(bead);
                    REQUIRE(candidates.size() <= 32);
                }
                struct SourceSection { std::uint32_t source; std::size_t filament; ExPolygons polygons; };
                std::array<std::vector<SourceSection>, 2> sections;
                auto limits = scenario.study.sections.geometry;
                limits.half_native_section_planes = true;
                for (unsigned side = 0; side < 2; ++side)
                    for (const auto& [id, beads] : candidates) {
                        const auto section = section_deposition(beads, zbase + (side ? .5e-9 : -.5e-9), 0., limits);
                        for (const auto& region : section.regions)
                            sections[side].push_back({id, region.material, region.footprint.native()});
                    }
                std::cout << "Near-boundary source contact witness: original nominal beads, Z=9.8mm +/-0.5nm;"
                          << " 1mm XY halo; positive area is NOT physical bonding or component membership\n";
                std::size_t pairs = 0, work = 0, uppers = 0, distance_pairs = 0;
                for (const auto& upper : sections[1])
                    for (const auto& polygon : upper.polygons) {
                        const auto box = get_extents(polygon);
                        if (box.max.x() * 1e-9 < xmin || box.min.x() * 1e-9 > xmax ||
                            box.max.y() * 1e-9 < ymin || box.min.y() * 1e-9 > ymax) continue;
                        REQUIRE(++uppers <= 32);
                        std::cout << " upper_source=" << upper.source << " filament=" << upper.filament
                                  << " bounds_nm=" << box.min.x() << ',' << box.min.y() << ".."
                                  << box.max.x() << ',' << box.max.y() << '\n';
                        for (const auto& lower : sections[0]) {
                            REQUIRE(++pairs <= 256);
                            REQUIRE(work < 20000000);
                            const auto overlap = exact_section_empty({polygon}, lower.polygons,
                                ExactSectionOperation::Intersection, 1024, std::min<std::size_t>(2000000, 20000000 - work));
                            work += overlap.work;
                            double gap_squared = 0.;
                            if (overlap.empty) {
                                gap_squared = std::numeric_limits<double>::infinity();
                                // With no positive-area intersection between
                                // validated polygons, boundary distance reduces
                                // to endpoint/segment distances. This floating
                                // estimate never replaces the exact predicate.
                                const auto upper_edges = polygon.lines();
                                const auto lower_edges = to_lines(lower.polygons);
                                for (const auto& a : upper_edges)
                                    for (const auto& b : lower_edges) {
                                        REQUIRE(++distance_pairs <= 1000000);
                                        gap_squared = std::min({gap_squared, a.distance_to_squared(b.a),
                                            a.distance_to_squared(b.b), b.distance_to_squared(a.a), b.distance_to_squared(a.b)});
                                    }
                            }
                            REQUIRE(std::isfinite(gap_squared));
                            const auto lower_box = get_extents(lower.polygons);
                            std::cout << "  lower_source=" << lower.source << " filament=" << lower.filament
                                      << " positive_area=" << !overlap.empty
                                      << " approximate_xy_gap_mm=" << std::sqrt(gap_squared) * 1e-6
                                      << " bounds_nm=" << lower_box.min.x() << ',' << lower_box.min.y() << ".."
                                      << lower_box.max.x() << ',' << lower_box.max.y() << '\n';
                        }
                    }
                REQUIRE(uppers > 0);
                std::cout << " complete_local_witness selected_beads=" << selected << " upper_candidates=" << uppers
                          << " pairs=" << pairs << " work=" << work << "; outside-halo contacts not assessed\n";
                return;
            }
            if (std::getenv("ORCA_SECTION_DETACHED_SOURCE_WITNESS")) {
                // Nominal admission preserves the ordered timed prefix. These
                // are conservative AABB candidates, never component membership.
                REQUIRE(nominal.segments.size() <= path.deposition.size());
                struct Group {
                    std::size_t count = 0;
                    unsigned first = std::numeric_limits<unsigned>::max(), last = 0;
                };
                using Key = std::tuple<unsigned, std::uint32_t, unsigned, unsigned>;
                std::map<Key, Group> groups;
                std::size_t omitted_assignments = 0;
                for (std::size_t i = 0; i < nominal.segments.size(); ++i) {
                    const auto& bead = nominal.segments[i];
                    const auto& entry = path.deposition[i];
                    const auto& move = result.moves.at(entry.move_index);
                    constexpr double grid_guard = 2e-9;
                    const double radius = .5 * bead.width_m + grid_guard;
                    if (std::max(bead.start_m.x(),bead.end_m.x()) + radius < .184901674 ||
                        std::min(bead.start_m.x(),bead.end_m.x()) - radius > .185762328 ||
                        std::max(bead.start_m.y(),bead.end_m.y()) + radius < .120555190 ||
                        std::min(bead.start_m.y(),bead.end_m.y()) - radius > .121391814) continue;
                    const double lo = std::min(bead.start_m.z() - bead.height_m,
                        bead.end_m.z() - bead.end_height_m.value_or(bead.height_m)) - grid_guard;
                    const double hi = std::max(bead.start_m.z(),bead.end_m.z()) + grid_guard;
                    for (unsigned band = 0; band < 2; ++band) {
                        if (hi < (band ? .0094 : .0098) || lo > (band ? .0098 : .0114)) continue;
                        const Key key{band,move.deposition_source_id,unsigned(move.extruder_id),unsigned(entry.role)};
                        auto found = groups.find(key);
                        if (found == groups.end()) {
                            if (groups.size() == 32) { ++omitted_assignments; continue; }
                            found = groups.emplace(key,Group{}).first;
                        }
                        ++found->second.count;
                        found->second.first = std::min(found->second.first,entry.gcode_id);
                        found->second.last = std::max(found->second.last,entry.gcode_id);
                    }
                }
                std::cout << "Detached source witness ONLY; no mechanical study: width-expanded AABB candidates"
                          << " with 2nm grid guard, NOT exact component membership; same XY target/support Z=9.8..11.4/9.4..9.8mm\n";
                for (const auto& [key, group] : groups) {
                    const auto [band,id,filament,role] = key;
                    std::cout << "  band=" << (band ? "below" : "target") << " source=" << id
                              << " filament=" << filament << " role=" << role << " count=" << group.count
                              << " gcode_ids=" << group.first << ".." << group.last;
                    const auto* source = result.deposition_provenance ? result.deposition_provenance->source(id) : nullptr;
                    if (source && source->logical_filament == filament)
                        std::cout << " object=" << source->object << " instance=" << source->instance
                                  << " layer=" << source->layer << " region=" << source->region
                                  << " height_class=" << unsigned(source->height);
                    else std::cout << " source_metadata=unavailable";
                    std::cout << '\n';
                }
                std::cout << "  omitted_candidate_assignments=" << omitted_assignments
                          << " (output group cap 32; no analysis geometry removed)\n";
                return;
            }
            path = {};
            for (auto& segment : nominal.segments) segment.material = 0;
            std::optional<SectionPartitionEnvelopes> partition;
            double partition_quantum_m = 0.;
            const char* geometry_quantum = std::getenv("ORCA_SECTION_STREAM_GEOMETRY_QUANTUM_UM");
            const char* mechanical_quantum = std::getenv("ORCA_SECTION_STREAM_MECHANICAL_QUANTUM_UM");
            std::vector<double> quanta;
            if (const char* list = std::getenv("ORCA_SECTION_STREAM_MECHANICAL_QUANTA_UM")) {
                REQUIRE(geometry_quantum == nullptr);
                REQUIRE(mechanical_quantum == nullptr);
                REQUIRE(std::getenv("ORCA_SECTION_STREAM_GEOMETRY_PROBE") == nullptr);
                REQUIRE(std::getenv("ORCA_SECTION_POINT_CONTACT_WITNESS") == nullptr);
                quanta = mechanical_quanta(list);
            }
            if (geometry_quantum) {
                REQUIRE(mechanical_quantum == nullptr);
                if (!std::getenv("ORCA_SECTION_POINT_CONTACT_WITNESS"))
                    REQUIRE(std::getenv("ORCA_SECTION_STREAM_GEOMETRY_PROBE") != nullptr);
            }
            if (mechanical_quantum) {
                REQUIRE(std::getenv("ORCA_SECTION_STREAM_GEOMETRY_PROBE") == nullptr);
                REQUIRE(std::getenv("ORCA_SECTION_POINT_CONTACT_WITNESS") == nullptr);
            }
            if (const char* quantum = geometry_quantum ? geometry_quantum : mechanical_quantum) {
                std::size_t consumed = 0;
                const std::string value(quantum);
                partition_quantum_m = std::stod(value, &consumed) * 1e-6;
                REQUIRE(consumed == value.size());
                REQUIRE(std::isfinite(partition_quantum_m));
                REQUIRE(partition_quantum_m >= 2e-9);
                REQUIRE(partition_quantum_m <= .001);
                partition = make_section_partition_envelopes(nominal.segments, scenario.bed_z_m, partition_quantum_m);
            }
            LocalSectionStreamSettings settings;
            settings.origin_m = scenario.origin_m;
            settings.basis_spacing_m = scenario.initial_basis_spacing_m;
            if (const char* planes = std::getenv("ORCA_SECTION_STREAM_Z_PLANES_MM"))
                settings.basis_z_planes_m = mechanical_z_planes(planes);
            settings.tile_size_m = *scenario.study.local_tile_size_m;
            settings.bed_z_m = scenario.bed_z_m;
            settings.maximum_height_m = scenario.study.max_section_height_m;
            settings.geometry.geometry = scenario.study.sections.geometry;
            settings.geometry.max_sections = scenario.study.max_total_sections;
            settings.geometry.max_boundary_points = scenario.study.max_total_boundary_points;
            settings.geometry.max_segment_visits = scenario.study.max_total_segment_visits;
            settings.mechanics = scenario.study.mechanics;
            settings.mechanics.allow_free_components = std::getenv("ORCA_SECTION_STREAM_ALLOW_FREE") != nullptr;
            settings.topology = scenario.study.topology;
            const auto explicit_budget = [](const char* name, std::size_t maximum) -> std::optional<std::size_t> {
                const char* budget = std::getenv(name);
                if (!budget) return std::nullopt;
                const std::string value(budget);
                REQUIRE_FALSE(value.empty());
                REQUIRE(std::all_of(value.begin(), value.end(), [](char c) { return c >= '0' && c <= '9'; }));
                const auto amount = std::stoull(value);
                REQUIRE(amount > 0);
                REQUIRE(amount <= maximum);
                return std::size_t(amount);
            };
            if (const auto entries = explicit_budget("ORCA_SECTION_STREAM_SPARSE_FACTOR",64000000)) {
                settings.mechanics.sparse_factor = SectionFactorLimits{};
                settings.mechanics.sparse_factor->max_entries = *entries;
            }
            if (const auto work = explicit_budget("ORCA_SECTION_STREAM_SPARSE_WORK",256000000)) {
                REQUIRE(settings.mechanics.sparse_factor.has_value());
                settings.mechanics.sparse_factor->max_work = *work;
            }
            const std::size_t integration_ceiling = settings.mechanics.sparse_factor ? 3200000000u : 600000000;
            if (const auto work = explicit_budget("ORCA_SECTION_STREAM_INTEGRATION_WORK", integration_ceiling)) {
                settings.integration = StreamIntegrationLimits{};
                settings.integration->max_work = *work;
            }
            if (const auto groups = explicit_budget("ORCA_SECTION_STREAM_INTEGRATION_GROUPS",
                    settings.mechanics.sparse_factor ? 125000u : 31250u)) {
                if (!settings.integration) settings.integration = StreamIntegrationLimits{};
                settings.integration->max_groups = *groups;
            }
            if (const auto points = explicit_budget("ORCA_SECTION_STREAM_INTEGRATION_POINTS",
                    settings.mechanics.sparse_factor ? 64000000u : 32000000u)) {
                if (!settings.integration) settings.integration = StreamIntegrationLimits{};
                settings.integration->max_total_points = *points;
                settings.geometry.max_boundary_points = *points;
            }
            if (const auto points = explicit_budget("ORCA_SECTION_STREAM_TOPOLOGY_POINTS", 128000000)) {
                settings.streaming_topology = StreamTopologyLimits{};
                settings.streaming_topology->max_processed_points = *points;
            }
            if (const auto work = explicit_budget("ORCA_SECTION_STREAM_TOPOLOGY_WORK", 400000000)) {
                if (!settings.streaming_topology) settings.streaming_topology = StreamTopologyLimits{};
                settings.streaming_topology->max_work = *work;
            }
            if (const auto entries = explicit_budget("ORCA_SECTION_STREAM_ENTRIES", 250000)) {
                if (!settings.integration) settings.integration = StreamIntegrationLimits{};
                if (!settings.streaming_topology) settings.streaming_topology = StreamTopologyLimits{};
                settings.integration->max_pieces = *entries;
                settings.integration->max_provisional_components = *entries;
                settings.streaming_topology->max_islands = *entries;
                settings.streaming_topology->max_slabs = *entries;
            }
            if (const auto pieces = explicit_budget("ORCA_SECTION_STREAM_INTEGRATION_PIECES",
                    settings.mechanics.sparse_factor ? 500000u : 250000u)) {
                if (!settings.integration) settings.integration = StreamIntegrationLimits{};
                settings.integration->max_pieces = *pieces;
            }
            if (const auto work = explicit_budget("ORCA_SECTION_STREAM_INDEX_WORK", 500000000))
                settings.geometry.max_index_work = *work;
            if (const auto planes = explicit_budget("ORCA_SECTION_STREAM_PLANES", 200000))
                settings.geometry.max_sections = *planes;
            const std::size_t grid_multiplier = settings.mechanics.sparse_factor ? 16 : 1;
            Eigen::Vector2d reporting_spacing_m = scenario.initial_basis_spacing_m.head<2>();
            if (const auto reporting_um = explicit_budget("ORCA_SECTION_REPORT_GRID_UM",1000000))
                reporting_spacing_m.setConstant(double(*reporting_um)*1e-6);
            if (const auto cells = explicit_budget("ORCA_SECTION_STREAM_CELLS", 4000 * grid_multiplier))
                settings.mechanics.max_cells = *cells;
            if (const auto nodes = explicit_budget("ORCA_SECTION_STREAM_NODES", 6000 * grid_multiplier))
                settings.mechanics.max_nodes = *nodes;
            if (const auto dofs = explicit_budget("ORCA_SECTION_STREAM_DOFS", 12000 * grid_multiplier))
                settings.mechanics.max_free_dofs = *dofs;
            if (const auto work = explicit_budget("ORCA_SECTION_STREAM_ASSEMBLY_WORK", 2000000 * grid_multiplier))
                settings.mechanics.max_work = *work;
            if (const auto triplets = explicit_budget("ORCA_SECTION_STREAM_TRIPLETS", 1200000 * grid_multiplier))
                settings.mechanics.max_triplets = *triplets;
            if (const auto points = explicit_budget("ORCA_SECTION_STREAM_GEOMETRY_POINTS", 128000000)) {
                REQUIRE(std::getenv("ORCA_SECTION_STREAM_GEOMETRY_PROBE") != nullptr);
                settings.geometry.max_boundary_points = *points;
                settings.geometry.max_sections = 200000;
            }
            std::cout << "Visitor-wrapper SINGLE fixed-basis UNQUALIFIED synthetic single-material diagnostic"
                      << " layer=" << layer << " end_s=" << end << " segments=" << nominal.segments.size()
                      << " half_native=" << settings.geometry.geometry.half_native_section_planes
                      << " original_bead_coordinates_unchanged=1 partition_quantum_m=" << partition_quantum_m << '\n';
            if (partition)
                std::cout << "Explicit quadrature approximation only, NOT geometry/convergence qualification: maximum_partition_shift_m="
                          << partition->maximum_shift_m << " collapsed_partition_intervals=" << partition->collapsed_intervals
                          << " expanded_zero_partition_endpoints=" << partition->expanded_zero_endpoints << '\n';
            const auto started = std::chrono::steady_clock::now();
            std::size_t checkpoints = 0;
            StreamIntegrationStatistics statistics;
            StreamTopologyStatistics topology_statistics;
            LocalSectionStatistics geometry_statistics;
            if (std::getenv("ORCA_SECTION_STREAM_GEOMETRY_PROBE")) {
                std::size_t polygons = 0, bands = 0;
                LocalSectionPacketVisitor visitor;
                visitor.begin_tile = [&](const LocalSectionTile&) {
                    if (geometry_statistics.completed_tiles % 32 == 0)
                        std::cout << "Geometry workload progress: " << geometry_statistics.completed_tiles
                                  << '/' << geometry_statistics.indexed_tiles << " tiles" << std::endl;
                };
                visitor.band = [&](const LocalSectionTile& packet) {
                    ++bands;
                    for (const auto& slab : packet.slabs) polygons += slab.footprint.native().size();
                };
                visitor.end_tile = [] {};
                try {
                    const auto geometry = visit_local_deposition_packets(nominal.segments,
                        settings.origin_m.head<2>(), settings.tile_size_m, settings.bed_z_m,
                        settings.maximum_height_m, visitor, settings.geometry, {}, {},
                        settings.interval_subdivisions, &geometry_statistics, partition ? &partition->segments : nullptr);
                    REQUIRE(geometry_statistics.completed_tiles == geometry_statistics.indexed_tiles);
                    std::cout << "Geometry-only workload, NOT connectivity or mechanics: tiles="
                              << geometry_statistics.completed_tiles << " bands=" << bands
                              << " polygons=" << polygons << " sections=" << geometry.sections
                              << " boundary_points=" << geometry.boundary_points
                              << " segment_visits=" << geometry.segment_visits
                              << " sampled_volume_m3=" << geometry.volume_m3
                              << " elapsed_s=" << std::chrono::duration<double>(
                                  std::chrono::steady_clock::now() - started).count() << '\n';
                } catch (...) {
                    std::cout << "Incomplete geometry workload: completed_tiles=" << geometry_statistics.completed_tiles
                              << " elapsed_s=" << std::chrono::duration<double>(
                                  std::chrono::steady_clock::now() - started).count() << '\n';
                    throw;
                }
                return;
            }
            if (const char* witness = std::getenv("ORCA_SECTION_POINT_CONTACT_WITNESS")) {
                const bool hole_contact = std::string(witness) == "hole";
                struct Captured {};
                LocalSectionPacketVisitor visitor;
                visitor.begin_tile = [](const LocalSectionTile&) {};
                visitor.end_tile = [] {};
                visitor.band = [hole_contact](const LocalSectionTile& packet) {
                    if (packet.index != (hole_contact ? std::array<int, 2>{6, 3} : std::array<int, 2>{4, 12})) return;
                    for (const auto& slab : packet.slabs) {
                        if (std::llround(slab.bottom_m * 1e9) != (hole_contact ? 440000 : 563094) ||
                            std::llround(slab.top_m * 1e9) != (hole_contact ? 450000 : 564000)) continue;
                        for (const auto& polygon : slab.footprint.native()) {
                            if (hole_contact) {
                                if (polygon.holes.size() < 3 || polygon.holes[2].points.size() != 18) continue;
                                std::cout << "Point-contact witness target hole2; total_holes=" << polygon.holes.size()
                                          << " outer_points=" << polygon.contour.points.size() << '\n';
                                for (const auto& point : polygon.holes[2].points)
                                    std::cout << point.x() << ',' << point.y() << ' ';
                                std::cout << '\n';
                                throw Captured{};
                            } else if (polygon.contour.points.size() != 620) continue;
                            REQUIRE(polygon.holes.size() <= 32);
                            std::size_t points = polygon.contour.points.size();
                            for (const auto& hole : polygon.holes) points += hole.points.size();
                            REQUIRE(points <= 4096);
                            const auto print_ring = [](const Polygon& ring) {
                                std::size_t column = 0;
                                for (const auto& point : ring.points) {
                                    std::cout << point.x() << ',' << point.y() << ' ';
                                    if (++column % 32 == 0) std::cout << '\n';
                                }
                                std::cout << '\n';
                            };
                            std::cout << "Point-contact witness outer:\n";
                            print_ring(polygon.contour);
                            for (const auto& hole : polygon.holes) {
                                std::cout << "Point-contact witness hole:\n";
                                print_ring(hole);
                            }
                            throw Captured{};
                        }
                    }
                };
                bool captured = false;
                try {
                    visit_local_deposition_packets(nominal.segments, settings.origin_m.head<2>(),
                        settings.tile_size_m, settings.bed_z_m, settings.maximum_height_m,
                        visitor, settings.geometry, {}, {}, settings.interval_subdivisions,
                        nullptr, partition ? &partition->segments : nullptr);
                } catch (const Captured&) { captured = true; }
                REQUIRE(captured);
                return;
            }
            const char* basis_mode = std::getenv("ORCA_SECTION_STREAM_BASIS_PROBE");
            const char* release_box = std::getenv("ORCA_SECTION_RELEASE_BOX_MM");
            const bool cell_contact = std::getenv("ORCA_SECTION_CELL_CONTACT") != nullptr;
            if (cell_contact) REQUIRE(release_box != nullptr);
            const bool release_refine_xy = std::getenv("ORCA_SECTION_RELEASE_REFINE_XY") != nullptr;
            if (release_refine_xy) REQUIRE(release_box != nullptr);
            const bool basis_probe = basis_mode != nullptr;
            std::vector<Eigen::Vector3d> bases;
            if (basis_probe) {
                // Graded Z is fixed across lateral probes. A Z-spacing change
                // would be ineffective and must instead insert explicit planes.
                const std::string mode(basis_mode);
                const bool valid_grid = settings.basis_z_planes_m.empty() || mode == "x" || mode == "y" || mode == "xy";
                REQUIRE(valid_grid);
                REQUIRE(quanta.empty());
                REQUIRE(mechanical_quantum != nullptr);
                quanta.assign(4, partition_quantum_m);
                bases = mechanical_bases(scenario.initial_basis_spacing_m,basis_mode);
            }
            std::vector<SectionBedReleaseRegion> release_regions;
            if (release_box) {
                REQUIRE_FALSE(basis_probe);
                REQUIRE(quanta.empty());
                REQUIRE(mechanical_quantum != nullptr);
                release_regions.push_back(mechanical_release_region(release_box));
                quanta.assign(release_refine_xy ? 4 : 2,partition_quantum_m);
                if (release_refine_xy) {
                    bases.assign(4,scenario.initial_basis_spacing_m);
                    bases[2].head<2>() *= .5;
                    bases[3].head<2>() *= .5;
                }
                std::cout << "Prescribed bed-release pair: fixed rectangle_mm=" << release_box
                          << "; NOT peel propagation or interface calibration\n";
            }
            const bool shared_prefix_comparison = !quanta.empty();
            if (quanta.empty()) quanta.push_back(partition_quantum_m);
            bool all_fields_completed = true;
            std::size_t completed_fields = 0;
            std::optional<double> bonded_energy;
            if (shared_prefix_comparison)
                std::cout << (basis_probe ? "Shared-prefix basis comparison: identical nominal segments=" :
                    "Shared-prefix fixed-basis quadrature comparison: identical nominal segments=")
                          << nominal.segments.size() << " runs=" << quanta.size()
                          << "; no re-slicing and NO convergence verdict\n";
            for (std::size_t run = 0; run < quanta.size(); ++run) {
                const auto run_started = std::chrono::steady_clock::now();
                auto last_report = run_started;
                statistics = {};
                topology_statistics = {};
                geometry_statistics = {};
                checkpoints = 0;
                partition_quantum_m = quanta[run];
                settings.basis_spacing_m = scenario.initial_basis_spacing_m;
                settings.released_bed_regions = release_box && run % 2 == 1 ? release_regions : std::vector<SectionBedReleaseRegion>{};
                settings.cell_contact.reset();
                if (cell_contact && run % 2 == 1) settings.cell_contact = SectionCellContact{};
                if (!bases.empty()) settings.basis_spacing_m = bases[run];
                std::cout << "Mechanical diagnostic row=" << run << " partition_quantum_m=" << partition_quantum_m
                          << " basis_m=" << settings.basis_spacing_m.transpose()
                          << " graded_z_planes=" << settings.basis_z_planes_m.size() << std::endl;
                if (!settings.basis_z_planes_m.empty()) {
                    std::cout << "Graded basis Z coordinates (m):";
                    for (double z : settings.basis_z_planes_m) std::cout << ' ' << z;
                    std::cout << '\n';
                }
                bool completed_mechanical_field = false;
                try {
                    if (shared_prefix_comparison && !basis_probe && (run == 0 || quanta[run] != quanta[run-1])) {
                        partition.reset();
                        partition = make_section_partition_envelopes(nominal.segments, scenario.bed_z_m, partition_quantum_m);
                        std::cout << "Shared-prefix partition: maximum_shift_m=" << partition->maximum_shift_m
                                  << " collapsed_intervals=" << partition->collapsed_intervals
                                  << " expanded_zero_endpoints=" << partition->expanded_zero_endpoints << '\n';
                    }
                    const auto solved = solve_local_deposition_sections(nominal.segments,
                        {scenario.materials.front()}, settings, {}, [&] {
                            if ((++checkpoints & 0xfffff) != 0) return;
                            const auto now = std::chrono::steady_clock::now();
                            if (now - last_report < std::chrono::seconds(10)) return;
                            last_report = now;
                            std::cout << "Mechanical workload progress: " << geometry_statistics.completed_tiles
                                      << '/' << geometry_statistics.indexed_tiles << " tiles, "
                                      << statistics.piece_count << " pieces" << std::endl;
                        }, partition ? &partition->segments : nullptr, {}, &statistics,
                        &topology_statistics, &geometry_statistics);
                    double tensile_generalized_force = 0.;
                    for (const auto& reaction : solved.field.generalized_bed_reaction_n)
                        tensile_generalized_force += std::max(0., -reaction.z());
                    REQUIRE(std::isfinite(tensile_generalized_force));
                    const auto reported = project_section_bed_reactions(solved.field,settings.origin_m.head<2>(),
                        reporting_spacing_m);
                    completed_mechanical_field = true;
                    if (release_box) {
                        if (run % 2 == 0) bonded_energy = solved.field.elastic_energy_j;
                        else if (bonded_energy) {
                            REQUIRE(solved.field.released_bed_area_m2 > 0.);
                            CHECK(solved.field.elastic_energy_j <= *bonded_energy*(1.+1e-8));
                            std::cout << "Prescribed release finite-area demand: released_area_m2=" << solved.field.released_bed_area_m2
                                      << " energy_drop_j=" << *bonded_energy-solved.field.elastic_energy_j
                                      << " average_release_j_m2=" << (*bonded_energy-solved.field.elastic_energy_j)/solved.field.released_bed_area_m2
                                      << " bed_lift_upper_bound_m=" << solved.field.bed_lift_upper_bound_m
                                      << " maximum_footprint_lift_m=" << solved.field.maximum_footprint_lift_m
                                      << " maximum_footprint_penetration_m=" << solved.field.maximum_footprint_penetration_m
                                      << " conservative_cell_contact=" << solved.field.conservative_cell_contact
                                      << " contact_iterations=" << solved.field.contact_iterations << '\n';
                        }
                    }
                    std::cout << "Stream diagnostic completed (NOT convergence or physical qualification): sections="
                              << solved.geometry.sections << " points=" << solved.geometry.boundary_points
                              << " partition_quantum_m=" << partition_quantum_m
                              << " sampled_volume_m3=" << solved.field.occupied_volume_m3
                              << " energy_j=" << solved.field.elastic_energy_j
                              << " tensile_generalized_bed_force_n=" << tensile_generalized_force
                              << " fixed_grid_generalized_force_n=" << reported.tensile_generalized_force_n
                              << " reporting_xy_m=" << reporting_spacing_m.transpose()
                              << " bed_area_m2=" << solved.field.bed_area_m2
                              << " relative_residual=" << solved.field.relative_residual
                              << " nodes=" << solved.field.nodes_m.size()
                              << " free_components=" << solved.field.free_components
                              << " affine_free_components=" << solved.field.affine_free_components
                              << " quadratic_free_components=" << solved.field.quadratic_free_components
                              << " backward_error=" << solved.field.componentwise_relative_residual
                              << " linear_refinement_steps=" << solved.field.linear_refinement_steps
                              << " reference_dofs=" << solved.field.reference_dofs
                              << " total_displacement_dofs=" << 3 * solved.field.nodes_m.size()
                              << " solved_dofs=" << solved.field.solved_dofs
                              << " factor_entries=" << solved.field.factor_entries
                              << " assembly_triplets=" << solved.field.assembly_triplets
                              << " matrix_entries=" << solved.field.matrix_entries
                              << " occupied_axis_cells=" << solved.field.occupied_axis_cells[0] << ','
                              << solved.field.occupied_axis_cells[1] << ',' << solved.field.occupied_axis_cells[2]
                              << " final_components=not_exposed free_dofs=not_exposed"
                              << "; generalized reaction is NOT local bed traction or adhesion capacity\n";
                } catch (const ResourceLimit& error) {
                    std::cout << "Stream diagnostic ResourceLimited; no field: " << error.what() << '\n';
                } catch (const std::invalid_argument& error) {
                    std::cout << "Stream diagnostic unsupported; no field: " << error.what() << '\n';
                } catch (const std::runtime_error& error) {
                    std::cout << "Stream diagnostic unresolved; no field: " << error.what() << '\n';
                }
                std::cout << "Stream diagnostic elapsed_s="
                          << std::chrono::duration<double>(std::chrono::steady_clock::now() - run_started).count()
                          << " checkpoints=" << checkpoints
                          << " processed_points=" << statistics.processed_points
                          << " integration_work=" << statistics.integration_work
                          << " assembly_work=" << statistics.assembly_work
                          << " peak_groups=" << statistics.peak_groups
                          << " peak_memberships=" << statistics.peak_memberships
                          << " peak_input_points=" << statistics.peak_input_points
                          << " peak_clip_points=" << statistics.peak_clip_points
                          << " piece_count=" << statistics.piece_count
                          << " provisional_mapping_size=" << statistics.provisional_mapping_size
                          << " topology_input_points=" << topology_statistics.input_points
                          << " topology_processed_points=" << topology_statistics.processed_points
                          << " topology_work=" << topology_statistics.work
                          << " topology_slabs=" << topology_statistics.slabs
                          << " topology_islands=" << topology_statistics.islands
                          << " topology_processed_islands=" << topology_statistics.processed_islands
                          << " topology_labels=" << topology_statistics.labels
                          << " topology_face_records=" << topology_statistics.face_records
                          << " topology_peak_face_records=" << topology_statistics.peak_face_records
                          << " topology_trace_intervals=" << topology_statistics.trace_intervals
                          << " topology_peak_live_points=" << topology_statistics.peak_live_points
                          << " completed_tiles=" << geometry_statistics.completed_tiles
                          << " indexed_tiles=" << geometry_statistics.indexed_tiles
                          << "; failed traversals intentionally publish no partial geometry metrics\n";
                all_fields_completed = all_fields_completed && completed_mechanical_field;
                completed_fields += completed_mechanical_field ? 1 : 0;
            }
            std::cout << "Mechanical diagnostic fields_completed=" << completed_fields << '/' << quanta.size()
                      << "; no convergence verdict; failed rows have no field\n";
            // An admission probe may intentionally stop at a bound, but a
            // completion check must not pass merely because it caught that stop.
            if (std::getenv("ORCA_SECTION_STREAM_REQUIRE_FIELD"))
                REQUIRE(all_fields_completed);
            return;
        }
        if (std::getenv("ORCA_OWNERSHIP_WITNESS")) {
            // Read-only offline investigation, with an independent bounded
            // budget. No diagnostic work is charged to or enabled in Auto.
            struct Captured : std::runtime_error { Captured() : std::runtime_error("captured ownership") {} };
            auto path = timed_deposition(result, 2000000, 2000000);
            audit_collapsed_envelopes(result, path);
            REQUIRE(path.state == ToolpathState::Complete);
            double end = 0.;
            for (const auto& entry : path.deposition)
                if (result.moves.at(entry.move_index).layer_id == layer && entry.role != erWipeTower &&
                    entry.role != erCustom && entry.role != erNone) end = std::max(end, entry.end_s);
            path = {};
            std::vector<DepositedSection::Region> native;
            ExPolygons first_claim, second_claim, first_native, second_native;
            std::size_t remaining_work = 20000000;
            bool captured = false;
            ConflictOwnershipLimits limits;
            limits.geometry = scenario.study.sections.geometry;
            limits.observe_masks = [&](double z, bool clipped, const auto& regions) {
                if (std::abs(z - .0024945) > 1e-12) return;
                if (!clipped) {
                    std::size_t points = 0;
                    for (const auto& region : regions) points += count_points(region.footprint.native());
                    REQUIRE(points <= 4500);
                    native = regions;
                    return;
                }
                for (std::size_t i = 0; i < regions.size(); ++i)
                    for (std::size_t j = 0; j < i; ++j) {
                        const auto& a = regions[i];
                        const auto& b = regions[j];
                        if (a.material == b.material || a.footprint.native().empty() || b.footprint.native().empty()) continue;
                        if (count_points(a.footprint.native()) + count_points(b.footprint.native()) > 1024) continue;
                        const auto tested = exact_section_empty(a.footprint.native(), b.footprint.native(), ExactSectionOperation::Intersection,
                            1024, std::min<std::size_t>(remaining_work, 4000000));
                        remaining_work -= tested.work;
                        if (tested.empty) continue;
                        first_claim = a.footprint.native();
                        second_claim = b.footprint.native();
                        for (const auto& region : native) {
                            if (region.material == a.material)
                                first_native.insert(first_native.end(), region.footprint.native().begin(), region.footprint.native().end());
                            if (region.material == b.material)
                                second_native.insert(second_native.end(), region.footprint.native().begin(), region.footprint.native().end());
                        }
                        captured = true;
                        throw Captured();
                    }
            };
            SampledConflictOwnershipResolver resolver;
            prepare_native_model_ownership(result, scenario.section, scenario.bed_z_m, scenario.filament_to_material,
                end, scenario.stationary_extrusion_allowance, limits, {}, false, &resolver);
            scenario.use_native_ownership = false;
            const auto ownership = [&](double z) -> std::optional<SectionOwnership> {
                SectionOwnership request{z, {}};
                request.resolve_sampled_conflict = [&](const SampledSectionConflict& conflict) { return resolver(conflict).regions; };
                return request;
            };
            try {
                const auto rejected = preview_section_layer(result, layer, scenario, 2000000, 2000000, 200000, ownership);
                REQUIRE(rejected.samples.empty());
                std::cout << "Offline capture stop: " << rejected.detail << '\n';
            }
            catch (const Captured&) { captured = true; }
            REQUIRE(captured);
            const auto box_a = get_extents(first_claim), box_b = get_extents(second_claim);
            const Point domain_min(std::max(box_a.min.x(), box_b.min.x()), std::max(box_a.min.y(), box_b.min.y()));
            const Point domain_max(std::min(box_a.max.x(), box_b.max.x()), std::min(box_a.max.y(), box_b.max.y()));
            const auto relevant = [&](const Polygon& ring) {
                const auto box = get_extents(ring);
                return box.max.x() >= domain_min.x() && box.min.x() <= domain_max.x() &&
                    box.max.y() >= domain_min.y() && box.min.y() <= domain_max.y();
            };
            const auto restrict_domain = [&](const ExPolygons& native) {
                ExPolygons filtered;
                for (const auto& polygon : native) {
                    if (!relevant(polygon.contour)) continue;
                    filtered.emplace_back(polygon.contour);
                    for (const auto& hole : polygon.holes)
                        if (relevant(hole)) filtered.back().holes.push_back(hole);
                }
                return filtered;
            };
            const auto first_local = restrict_domain(first_native), second_local = restrict_domain(second_native);
            std::cout << "Offline original-coordinate domain pruning native_vertices=" << count_points(first_local)
                      << ',' << count_points(second_local) << " claim_vertices=" << count_points(first_claim)
                      << ',' << count_points(second_claim) << '\n';
            const auto certificate = exact_intersection_ownership(first_claim, second_claim, first_local, second_local);
            REQUIRE(certificate.work <= remaining_work);
            remaining_work -= certificate.work;
            std::cout << "Offline whole-overlap certificate nonempty=" << certificate.nonempty
                      << " unique_native=" << certificate.first_unique << ',' << certificate.second_unique
                      << " native_partition=" << certificate.native_partition
                      << " native_covers=" << certificate.native_covers << " native_disjoint=" << certificate.native_disjoint
                      << " total_work=" << 20000000 - remaining_work << "; no integer correction applied\n";
            const auto variants = partition_uncertain_interface(first_claim, second_claim, first_local, second_local);
            REQUIRE(variants.work <= remaining_work);
            remaining_work -= variants.work;
            std::cout << "Offline unqualified interface variants work=" << variants.work
                      << " allowance_nm=" << variants.interface_allowance_nm
                      << " representation_allowance_nm=" << variants.representation_allowance_nm
                      << " total_work=" << 20000000 - remaining_work
                      << "; full-region topology/contact validation still required, no field accepted\n";
            return;
        }
        const auto started = std::chrono::steady_clock::now();
        const auto report = preview_section_layer(result, layer, scenario);
        std::cout << "Static layer computational audit z_mm=" << requested_z_mm << " layer=" << layer
                  << " state=" << int(report.state) << " samples=" << report.samples.size()
                  << " tiles=" << report.geometry_tiles << " slabs=" << report.geometry_slabs
                  << " planes=" << report.section_count << " visits=" << report.segment_visits
                  << " points=" << report.boundary_points
                  << " seconds=" << std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()
                  << " reason=" << report.detail << '\n';
        for (const auto& level : report.basis_history)
            std::cout << "  basis spacing_m=" << level.spacing_m.transpose() << " nodes=" << level.nodes
                      << " energy_J=" << level.elastic_energy_j << " generalized_tensile_N=" << level.tensile_generalized_bed_force_n << '\n';
        for (const auto& level : report.geometry_history)
            std::cout << "  geometry z_quantization_m=" << level.z_quantization_m
                      << " maximum_partition_shift_m=" << level.maximum_partition_shift_m
                      << " maximum_boundary_uncertainty_m=" << level.maximum_boundary_uncertainty_m
                      << " added_representation_volume_bound_m3=" << level.added_representation_volume_bound_m3
                      << " removed_representation_volume_bound_m3=" << level.removed_representation_volume_bound_m3
                      << " coalesced_endpoint_intervals=" << level.collapsed_endpoint_intervals
                      << " expanded_zero_endpoints=" << level.expanded_zero_endpoints
                      << " volume_m3=" << level.volume_m3 << " components=" << level.components << '\n';
        if (report.state != SectionLayerPreview::State::Ready) CHECK(report.samples.empty());
        std::cout << "Synthetic prescribed-strain audit is not a physical adhesion or thermal-history validation.\n";
    }
    std::cout << "This audits source attribution and optional nominal sections, not mechanical validity.\n";
}
