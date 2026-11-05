"""Check planar arc-fill crossings in emitted G-code, including fitted G2/G3.

This checks centerline intersections, not full nozzle clearance or bead coverage.
Endpoint tolerance applies to complete extrusion paths, never tessellation chords.
"""
import argparse
from collections import defaultdict
import json
import math

from audit_gcode_moves import WORDS, linearized_gcode


def audit_crossings(stream, chord_error=0.0005, endpoint_tolerance=0.005):
    if not math.isfinite(endpoint_tolerance) or endpoint_tolerance < 0:
        raise ValueError("endpoint tolerance must be finite and nonnegative")
    position = dict.fromkeys("XYZ", 0.)
    role = ""
    active = None
    paths, segments, candidates = [], [], []
    grid = defaultdict(list)
    arc_roles = {"Arc bridge", "Arc overhang"}
    for source, line in linearized_gcode(stream, chord_error, exclude_cleaner_extrusion=True):
        if line.startswith(";TYPE:"):
            new_role = line.strip()[6:]
            if new_role != role:
                active = None
            role = new_role
        opcode = line.split(maxsplit=1)[0] if line.strip() else ""
        if opcode not in ("G0", "G1", "G92"):
            continue
        words = {key: float(value) for key, value in WORDS.findall(line.split(";")[0])}
        a = tuple(position.values())
        position.update((key, words[key]) for key in "XYZ" if key in words)
        b = tuple(position.values())
        if opcode == "G92":
            if any(key in words for key in "XYZ"):
                active = None
            continue
        if words.get("E", 0) <= 0 or math.dist(a, b) < 1e-9 or abs(a[2] - b[2]) > 1e-6:
            active = None
            continue
        if active is None:
            active = len(paths)
            paths.append([a[:2], b[:2]])
        else:
            paths[active][1] = b[:2]
        cells = [(round(b[2], 6), x, y)
                 for x in range(math.floor(min(a[0], b[0])), math.floor(max(a[0], b[0])) + 1)
                 for y in range(math.floor(min(a[1], b[1])), math.floor(max(a[1], b[1])) + 1)]
        seen = set()
        r = (b[0] - a[0], b[1] - a[1])
        for cell in cells:
            for idx in grid[cell]:
                if idx in seen:
                    continue
                seen.add(idx)
                c, d, other_role, other_source, other_path = segments[idx]
                if role not in arc_roles and other_role not in arc_roles:
                    continue
                if active == other_path and idx == len(segments) - 1:
                    continue
                s = (d[0] - c[0], d[1] - c[1])
                determinant = r[0] * s[1] - r[1] * s[0]
                if abs(determinant) < 1e-12:
                    continue  # Collinear bead overlap is a separate check.
                q = (c[0] - a[0], c[1] - a[1])
                t = (q[0] * s[1] - q[1] * s[0]) / determinant
                u = (q[0] * r[1] - q[1] * r[0]) / determinant
                if 0 <= t <= 1 and 0 <= u <= 1:
                    point = (a[0] + t * r[0], a[1] + t * r[1])
                    candidates.append((active, other_path, point, source, other_source, b[2]))
        for cell in cells:
            grid[cell].append(len(segments))
        segments.append((a, b, role, source, active))

    crossings, seen = [], set()
    for first, second, point, source, other_source, z in candidates:
        distances = [min(math.dist(point, end) for end in paths[idx]) for idx in (first, second)]
        if min(distances) <= endpoint_tolerance:
            continue
        key = (min(first, second), max(first, second), *(round(value, 6) for value in point))
        if key in seen:
            continue
        seen.add(key)
        crossings.append({"lines": [source, other_source], "xy": point, "z": z,
                          "endpoint_distances": distances})
    return {"planar_segments": len(segments), "crossings": crossings}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gcode")
    parser.add_argument("--chord-error", type=float, default=0.0005)
    parser.add_argument("--endpoint-tolerance", type=float, default=0.005)
    args = parser.parse_args()
    with open(args.gcode) as stream:
        result = audit_crossings(stream, args.chord_error, args.endpoint_tolerance)
    print(json.dumps(result, indent=2))
    raise SystemExit(bool(result["crossings"]))
