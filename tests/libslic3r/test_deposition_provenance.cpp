#include <catch2/catch_test_macros.hpp>

#include "libslic3r/GCode/DepositionProvenance.hpp"

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Slic3r;

TEST_CASE("Deposition sources retain explicit wiping assignment without inferring purity", "[DepositionProvenance]")
{
    DepositionProvenance registry;
    DepositionSource source;
    const auto unknown=registry.add(source);
    source.assignment=DepositionSource::Assignment::NominalRole;
    const auto normal=registry.add(source);
    source.assignment=DepositionSource::Assignment::WipingOverride;
    const auto wiping=registry.add(source);
    REQUIRE(registry.source(unknown));
    REQUIRE(registry.source(normal));
    REQUIRE(registry.source(wiping));
    CHECK(registry.source(unknown)->assignment==DepositionSource::Assignment::Unknown);
    CHECK(registry.source(normal)->assignment==DepositionSource::Assignment::NominalRole);
    CHECK(registry.source(wiping)->assignment==DepositionSource::Assignment::WipingOverride);
    source.assignment=static_cast<DepositionSource::Assignment>(255);
    CHECK_THROWS_AS(registry.add(source),std::invalid_argument);
}

TEST_CASE("Output provenance capture waits for a sealed producer registry", "[DepositionProvenance]")
{
    auto registry = std::make_shared<DepositionProvenance>();
    DepositionProvenanceCapture capture(registry->generation());
    const auto id = registry->add({});
    CHECK(capture.consume(registry->begin(id)));
    CHECK(capture.source_id() == id);
    CHECK_FALSE(capture.consume("G1 X1 E1\n"));
    CHECK(capture.consume(registry->clear()));
    CHECK(capture.source_id() == 0);
    registry->seal();
    CHECK(capture.finish(registry, {id, 0}));
    CHECK(capture.source_id() == 0);
    CHECK_FALSE(capture.finish(registry, {id, 0}));
}

TEST_CASE("Output provenance capture fails closed on missing unknown or over-budget events", "[DepositionProvenance]")
{
    auto registry = std::make_shared<DepositionProvenance>();
    const auto id = registry->add({});
    registry->seal();
    DepositionProvenanceCapture omitted(registry->generation());
    CHECK_FALSE(omitted.finish(registry, {id, 0}));
    DepositionProvenanceCapture unknown(registry->generation());
    CHECK(unknown.consume(";_ORCA_INTERNAL_SOURCE:" + std::to_string(registry->generation()) + ":999\n"));
    CHECK(unknown.consume(registry->clear()));
    CHECK_FALSE(unknown.finish(registry, {}));
    DepositionProvenanceCapture bounded(registry->generation(), 1);
    CHECK(bounded.consume(registry->begin(id)));
    CHECK(bounded.consume(registry->clear()));
    CHECK_FALSE(bounded.valid());
    CHECK(bounded.source_id() == 0);
    CHECK_FALSE(bounded.finish(registry, {id, 0}));
    DepositionProvenanceCapture nested(registry->generation());
    nested.consume(registry->begin(id));
    nested.consume(registry->begin(id));
    CHECK_FALSE(nested.valid());
    CHECK(nested.source_id() == 0);
    DepositionProvenanceCapture malformed(registry->generation());
    CHECK(malformed.consume(";_ORCA_INTERNAL_SOURCE malformed\n"));
    CHECK_FALSE(malformed.valid());
    CHECK(malformed.consume(registry->begin(id)));
    CHECK(malformed.source_id() == 0);
    auto stale = std::make_shared<DepositionProvenance>();
    DepositionProvenanceCapture wrong_generation(stale->generation());
    CHECK(wrong_generation.consume(registry->begin(id)));
    CHECK_FALSE(wrong_generation.valid());
    CHECK_THROWS_AS(DepositionProvenanceCapture(0), std::invalid_argument);
    CHECK_THROWS_AS(DepositionProvenanceCapture(registry->generation(), 0), std::invalid_argument);
    CHECK_THROWS_AS(DepositionProvenanceCapture(registry->generation(), DepositionProvenanceDecoder::maximum_events + 1),
                    std::invalid_argument);
}

TEST_CASE("Expected deposition transcripts reject missing duplicated and reordered scopes", "[DepositionProvenance]")
{
    auto registry     = std::make_shared<DepositionProvenance>();
    const auto first  = registry->add({});
    const auto second = registry->add({1});
    registry->seal();
    const std::vector<std::uint32_t> expected{first, 0, second, 0};
    const auto run = [&](const std::vector<std::uint32_t>& events) {
        DepositionProvenanceDecoder decoder(registry, expected);
        CHECK_FALSE(decoder.transcript_verified());
        for (const auto id : events)
            CHECK(decoder.consume(id == 0 ? registry->clear() : registry->begin(id)));
        decoder.finish();
        return decoder.transcript_verified();
    };
    CHECK(run(expected));
    CHECK_FALSE(run({first, 0}));
    CHECK_FALSE(run({second, 0}));
    CHECK_FALSE(run({}));
    CHECK_FALSE(run({first, 0, first, 0, second, 0}));
    CHECK_FALSE(run({second, 0, first, 0}));
}

TEST_CASE("Expected deposition transcripts are bounded immutable and explicitly finalized", "[DepositionProvenance]")
{
    auto registry = std::make_shared<DepositionProvenance>();
    const auto id = registry->add({});
    registry->seal();
    CHECK_THROWS_AS((DepositionProvenanceDecoder(registry, {id})), std::invalid_argument);
    CHECK_THROWS_AS((DepositionProvenanceDecoder(registry, {id, id, 0})), std::invalid_argument);
    CHECK_THROWS_AS((DepositionProvenanceDecoder(registry, {id + 1, 0})), std::invalid_argument);
    const std::vector<std::uint32_t> too_many(DepositionProvenanceDecoder::maximum_events + 1, 0);
    CHECK_THROWS_AS((DepositionProvenanceDecoder(registry, too_many)), std::length_error);
    std::vector<std::uint32_t> expected{id, 0};
    DepositionProvenanceDecoder decoder(registry, expected);
    expected.clear();
    decoder.consume(registry->begin(id));
    decoder.consume(registry->clear());
    CHECK_FALSE(decoder.transcript_verified());
    decoder.finish();
    CHECK(decoder.transcript_verified());
    decoder.consume(registry->begin(id));
    CHECK_FALSE(decoder.valid());
    CHECK_FALSE(decoder.transcript_verified());
    DepositionProvenanceDecoder unchecked(registry);
    unchecked.finish();
    CHECK(unchecked.valid());
    CHECK_FALSE(unchecked.transcript_verified());
    DepositionProvenanceDecoder empty(registry, {});
    empty.finish();
    CHECK(empty.transcript_verified());
}

TEST_CASE("Local deposition sources retain identity without changing ordinary output", "[DepositionProvenance]")
{
    auto registry = std::make_shared<DepositionProvenance>();
    const auto id = registry->add({7, 2, 11, 3, 4, DepositionSource::Height::Transformed});
    registry->seal();
    DepositionProvenanceDecoder decoder(registry);
    CHECK_FALSE(decoder.source_id(4).has_value());
    REQUIRE(decoder.consume(registry->begin(id)));
    REQUIRE(decoder.source_id(4).has_value());
    CHECK(*decoder.source_id(4) == id);
    CHECK(registry->source(id)->object == 7);
    CHECK(registry->source(id)->instance == 2);
    CHECK(registry->source(id)->layer == 11);
    CHECK(registry->source(id)->region == 3);
    CHECK(registry->source(id)->height == DepositionSource::Height::Transformed);
    const std::string ordinary = "G1 X3 E1 ; unmodified\n";
    CHECK_FALSE(decoder.consume(ordinary));
    CHECK(ordinary == "G1 X3 E1 ; unmodified\n");
    CHECK(decoder.consume(registry->clear()));
    CHECK_FALSE(decoder.source_id(4).has_value());
    decoder.finish();
    CHECK(decoder.valid());
}

TEST_CASE("Deposition provenance rejects stale malformed and nested contexts", "[DepositionProvenance]")
{
    auto registry = std::make_shared<DepositionProvenance>();
    const auto id = registry->add({});
    registry->seal();
    auto other          = std::make_shared<DepositionProvenance>();
    const auto other_id = other->add({});
    other->seal();
    const std::string prefix = ";_ORCA_INTERNAL_SOURCE:" + std::to_string(registry->generation()) + ":";
    for (const std::string& invalid : {other->begin(other_id), prefix + "999\n", prefix + "-1\n", prefix + "4294967296\n",
                                       prefix + "1 trailing\n", std::string(";_ORCA_INTERNAL_SOURCE:bad\n")}) {
        DepositionProvenanceDecoder decoder(registry);
        CHECK(decoder.consume(invalid));
        CHECK_FALSE(decoder.valid());
        CHECK(decoder.consume(registry->begin(id)));
        CHECK_FALSE(decoder.source_id(0).has_value());
    }
    DepositionProvenanceDecoder nested(registry);
    nested.consume(registry->begin(id));
    nested.consume(registry->begin(id));
    CHECK_FALSE(nested.valid());
}

TEST_CASE("Deposition context cannot outlive its span or change logical filament", "[DepositionProvenance]")
{
    auto registry = std::make_shared<DepositionProvenance>();
    const auto id = registry->add({});
    registry->seal();
    DepositionProvenanceDecoder unfinished(registry);
    unfinished.consume(registry->begin(id));
    unfinished.finish();
    CHECK_FALSE(unfinished.valid());
    CHECK_FALSE(unfinished.source_id(0).has_value());
    DepositionProvenanceDecoder mismatch(registry);
    mismatch.consume(registry->begin(id));
    CHECK_FALSE(mismatch.source_id(1).has_value());
    CHECK_FALSE(mismatch.valid());
}

TEST_CASE("Deposition registry is bounded and immutable before decoding", "[DepositionProvenance]")
{
    CHECK_THROWS_AS(DepositionProvenance(0), std::invalid_argument);
    CHECK_THROWS_AS(DepositionProvenance(DepositionProvenance::maximum_sources + 1), std::invalid_argument);
    auto registry = std::make_shared<DepositionProvenance>(1);
    CHECK_THROWS_AS(DepositionProvenanceDecoder{registry}, std::invalid_argument);
    registry->add({});
    CHECK_THROWS_AS(registry->add({}), std::length_error);
    CHECK_THROWS_AS(registry->begin(0), std::invalid_argument);
    registry->seal();
    CHECK_THROWS_AS(registry->add({}), std::logic_error);
    CHECK(registry->source(0) == nullptr);
    CHECK(registry->source(2) == nullptr);
}
