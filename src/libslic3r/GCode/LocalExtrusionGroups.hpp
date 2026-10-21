#pragma once

#include "../BoundingBox.hpp"
#include <algorithm>
#include <cstddef>
#include <numeric>
#include <vector>

namespace Slic3r {

// Conservative spatial components: overlapping path bounds must stay together.
// Inputs are whole ordered extrusion entities, never individual wall segments.
// This deliberately merges some concave/disjoint footprints rather than splitting
// a wall stack or moving its fill to a different visit.
inline std::vector<std::vector<size_t>> local_extrusion_groups(const std::vector<BoundingBox> &boxes)
{
    std::vector<size_t> parent(boxes.size()), order(boxes.size());
    std::iota(parent.begin(), parent.end(), 0);
    std::iota(order.begin(), order.end(), 0);
    auto root = [&](size_t i) {
        while (parent[i] != i) {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        return i;
    };
    std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        return boxes[a].min.x() < boxes[b].min.x();
    });
    for (size_t i = 0; i < order.size(); ++i)
        for (size_t j = i + 1; j < order.size() && boxes[order[j]].min.x() <= boxes[order[i]].max.x(); ++j)
            if (boxes[order[i]].overlap(boxes[order[j]]))
                parent[root(order[j])] = root(order[i]);

    std::vector<std::vector<size_t>> groups;
    std::vector<size_t> group_ids(boxes.size(), boxes.size());
    for (size_t i = 0; i < boxes.size(); ++i) {
        const size_t r = root(i);
        if (group_ids[r] == boxes.size()) {
            group_ids[r] = groups.size();
            groups.emplace_back();
        }
        groups[group_ids[r]].push_back(i);
    }
    return groups;
}

} // namespace Slic3r
