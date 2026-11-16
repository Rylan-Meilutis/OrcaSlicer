#pragma once

#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/Polygon.hpp"

namespace Slic3r {
class ModelObject;

enum class PaintedSupportStyle { Default, Normal, Tree };

struct PaintedSupportContacts {
    ExPolygons normal;
    ExPolygons tree;
    ExPolygons automatic;
};

bool has_typed_support_paint(const ModelObject &object);
bool has_normal_support_paint(const ModelObject &object);

// Blockers take precedence; coincident projections of different painted faces
// assign the overlap to normal columns. All outputs are disjoint demand subsets.
PaintedSupportContacts partition_painted_support_contacts(const ExPolygons &demand,
    const Polygons &normal, const Polygons &tree, const Polygons &blockers);
}
