#pragma once

#include <algorithm>
#include <cstddef>
#include <string>

namespace Slic3r {
// TriangleSelector writes a depth-first nibble tree in reverse string order.
// Preserve subdivisions while mapping typed support leaves to legacy ENFORCER.
inline bool legacy_support_paint(const std::string &typed, std::string &legacy)
{
    legacy.clear();
    size_t remaining = typed.size(), pending = 1;
    const auto nibble = [&]() -> int {
        if (remaining == 0) return -1;
        const char c = typed[--remaining];
        return c >= '0' && c <= '9' ? c - '0' : c >= 'A' && c <= 'F' ? c - 'A' + 10 : -1;
    };
    while (pending != 0) {
        const int code = nibble();
        if (code < 0) return false;
        --pending;
        const int split = code & 3;
        if (split != 0) {
            if (split != 3 && (code >> 2) == 3) return false;
            pending += size_t(split + 1);
            legacy += "0123456789ABCDEF"[code];
        } else if (code == 12) {
            // Extended leaf codes zero and one represent support states 3 and 4.
            const int extended = nibble();
            if (extended != 0 && extended != 1) return false;
            legacy += '4';
        } else {
            legacy += "0123456789ABCDEF"[code];
        }
    }
    if (remaining != 0) return false; // Exactly one complete triangle tree.
    std::reverse(legacy.begin(), legacy.end());
    return true;
}

inline std::string select_support_paint(const std::string &legacy, const std::string &versioned)
{
    std::string downgraded;
    if (versioned.compare(0, 2, "1:") == 0 && legacy_support_paint(versioned.substr(2), downgraded))
        return versioned.substr(2);
    return legacy;
}

} // namespace Slic3r
